classdef OrbitPropagator < Propagator
    %ORBITPROPAGATOR Parent module for a class of spacecraft orbit
    %propagators. Notable subclasses include LunarPropagator and
    %EarthPropagator.
    
    properties
        % struct containing info about the primary body
        pri     (1,1)   struct
        % array of structs containing info about additional planets to consider
        sec     (1,:)   struct
        % starting time of sim, seconds past J2000
        t0
        % starting states of satellites @ t0 in J2000 frame
        x0
        % dimension of state
        dim             = 6
        % max degree/order of gravity model to use
        ord     (1,1)   {mustBeInteger,mustBePositive} = 1
        % number of satellites
        nsats   (1,1)   {mustBePositive,mustBeInteger} = 1
        % ODE45 options
        opts    (1,1)   struct
        % info stored from latest run (if frame == '', no run yet)
        ts
        xs
        frame   (1,:)   char = ''

        % solar radiation pressure variables %
        % coefficient of reflectivity
        Cr      (1,1)   double {mustBeNonnegative} = 1
        % SRP area-to-mass ratio
        Am      (1,1)   double {mustBeNonnegative} = 0

        % preallocation variables %
        % are transformations preallocated?
        pre     (1,1)   = 0
        % time span covering preallocation period
        t_pre   (1,:)   double
        % inertial to body-fixed transformations for t_pre (used in
        % spherical harmonics)
        q_I2F   (1,:)   quaternion
    end
    properties(Access=private)
        % precomputed factorial for harmonics calculations %
        fact (:,:) = 0
    end
    properties (Constant, Access=private)
        % solar radiation pressure constants %
        % kN/m^2 (Pa), solar radiation pressure at Earth, given by the solar
        % constant at Earth (W/m^2) divided by the speed of light (m/s)
        Ps = 1361/299792458 * 1e-3
        % km, 1 astronomical unit (average Earth-sun distance)
        AU = 149597870.700
        % km, radius of the sun
        Rs = 695700
        % km, radius of the Earth
        Re = 6378.1366
        % km, radius of the moon
        Rm = 1737.4
    end
    
    methods
        function obj = OrbitPropagator(t0,ord,varargin)
            %ORBITPROPAGATOR Construct an OrbitPropagator instance.
            %   Inputs:
            %    - t0; character string, 'DD-MMM-YYYY XX:XX:XX'
            %    - ord; maximum degree and order of gravity model to use
            %    - opts; optional name-value arg, ODE45 integration tolerances
            %    - Cr; optional name-value arg, coefficient of reflectivity
            %       for computing solar radiation pressure
            %    - A/m; optional name-value arg, area-to-mass ratio for
            %       computing solar radiation pressure
            %    - preallocate; optional name-value arg, used to
            %       preallocate transformation matrices so SPICE isn't called
            %       during propagation. Can be used to speed up (maybe?)
            %       propagation, but main intent is to be able to parallelize
            %       calls. In this case, the argument is then time steps of
            %       period in seconds past J2000. This range should
            %       encompass any propagation intervals called in the
            %       future, otherwise errors will occur.
            arguments
                t0      (1,:)
                ord     (1,1)   {mustBeInteger,mustBePositive}
            end
            arguments (Repeating)
                varargin
            end

            % manage basic arguments
            obj.ord = ord;
            obj.opts = odeset("RelTol", 1e-9, "AbsTol", 1e-11);
            obj.t_pre = [];
            default = 2;
            try
                varargin = varargin{1};
            catch
                varargin = {};
            end

            if nargin > default
                for i=1:2:length(varargin)
                    if strcmp(varargin{i}, "opts")
                        obj.opts = varargin{i+1};
                    elseif strcmp(varargin{i}, "Cr")
                        obj.Cr = varargin{i+1};
                    elseif strcmp(varargin{i}, "A/m")
                        obj.Am = varargin{i+1};
                    elseif strcmp(varargin{i}, "preallocate")
                        obj.pre = 1;
                        obj.t_pre = varargin{i+1};
                    elseif ~isempty(varargin)
                        error("OrbitPropagator:invalidArgument", ...
                            "%s is not a valid optional argument.", ...
                            convertCharsToStrings(varargin{i}));
                    end
                end
            elseif nargin < default
                error("OrbitPropagator:nargin", "Too few arguments.");
            end

            % parse start time
            if isa(t0, 'double')
                obj.t0 = t0;
            else
                obj.t0 = cspice_str2et(t0);
            end

            % precompute factorial for harmonics
            n = 1:ord+1; m = 1:ord+1;
            m = m(2:end);
            obj.fact = factorial(abs(n'-m+2))./factorial(abs(n'-m));
        end
        
        function [ts,xs,fail] = run(obj,tf,n,frame,assign)
            %RUN Propagate the stored states for tf seconds (n steps
            %between). Data returned in provided frame.
            %   Input:
            %    - tf; final time, seconds past t0
            %    - n; number of time steps
            %    - frame; reference frame to return data in
            %    - assign; optional argument (default true), should
            %         propagation data be assigned to object properties
            arguments
                obj     (1,1)   OrbitPropagator
                tf      (1,1)   double {mustBePositive}
                n       (1,1)   {mustBeInteger,mustBePositive}
                frame   (1,:)   char
                assign  (1,1) = true
            end

            ts = linspace(0, tf, n);
            [ts,xs,fail] = obj.runat(ts,frame,assign);
        end

        function [ts,xs,fail] = runat(obj,ts,frame,assign)
            %RUNAT Propagate the stored states over the provided time steps. 
            %Data returned in provided frame.
            %   Input:
            %    - ts; eval time steps, seconds past t0
            %    - frame; reference frame to return data in
            %    - assign; optional argument (default true), should
            %         propagation data be assigned to object properties
            arguments
                obj     (1,1)   OrbitPropagator
                ts      (1,:)   double
                frame   (1,:)   char
                assign  (1,1) = true
            end

            n = length(ts);
            ts = obj.t0 + ts;
            xs = zeros(6,n,obj.nsats);
            fail = 0;       % flag if propagation failed
            
            for i=1:obj.nsats
                [~,X] = ode89(@obj.dynamics, [obj.t0 ts], obj.x0(:,i), obj.opts);
                X = X(2:end,:)';
                % catch and bounce if propagation failed
                if size(X,2) ~= n, fail = 1; return; end

                if ~strcmp(frame, 'J2000')  % transform if necessary
                    for j=1:length(ts)
                        X(:,j) = cspice_sxform('J2000', frame, ts(j)) * X(:,j);
                    end
                end

                xs(:,:,i) = X;              % assign to output
            end

            if assign       % only overwrite object properties if requested
                obj.ts = ts;
                obj.xs = xs;
                obj.frame = frame;
            end
        end

        function [ts,xs,fail] = run_prealloc(obj,x0_,ts)
            %RUN_PREALLOC Propagate the input states over the provided time
            %steps. Designed for use in parallelized setups, so read
            %extended description for details.
            %   CAUTION: this is designed pretty much exclusively for
            %   parallelized work, causing several changes from
            %   OrbitPropagator.runat:
            %    - starting state is passed in to avoid assignment
            %    - ODE45 is used instead of ODE89 for speed
            %    - results are not assigned back to the object
            %    - data assumed given/returned in J2000
            %    - ts are seconds past J2000 rather than past t0
            %    - generally a bunch of features are stripped out for
            %      efficiency
            %
            %   Input:
            %    - x0_; starting states of satellites
            %    - ts; eval time steps, seconds past J2000
            %    - frame; reference frame to return data in
            %    - assign; optional argument (default true), should
            %         propagation data be assigned to object properties
            arguments
                obj     (1,1)   OrbitPropagator
                x0_     (6,:,:) double
                ts      (1,:)   double
            end

            n = length(ts);
            xs = zeros(6,n,obj.nsats);
            fail = 0;       % flag if propagation failed
            
            for i=1:obj.nsats
                [~,X] = ode45(@obj.dynamics, ts, x0_(:,i), obj.opts);
                % catch and bounce if propagation failed
                if size(X,1) ~= n, fail = 1; return; end

                xs(:,:,i) = X';     % assign to output
            end
        end

        function dxdt = dynamics(obj,t,x)
            %DYNAMICS Inertial dynamics for a satellite orbiting a primary
            %body with gravitational influence from a secondary (and tertiary);
            %includes nonspherical gravity.
            %   Input:
            %    - t; simulation time in seconds past J2000
            %    - x; state [pos (km); vel (km/s)] of spacecraft
            %    - pri; struct, {GM: gravitational parameter, x: @(t) position [km]}
            %           for primary body
            %    - N; maximum degree and order of harmonics to compute
            %    - sec; array of struct, {GM: gravitational parameter, x: @(t) position [km]}
            %           for secondary bodies
            arguments
                obj     (1,1)   OrbitPropagator
                t       (1,1)   double
                x       (6,1)   double
            end

            dxdt = zeros(6,1);
            x_1s = x(1:3);
            r_1s = norm(x_1s);
            
            % get transformation from inertial to body-fixed
            if obj.pre
                T = obj.rotate_prealloc(t, obj.q_I2F);
            else
                T = cspice_pxform('J2000', obj.pri.frame, t);
            end

            % get nonspherical gravity effects
            x_me = T * x_1s;
            f_ns = obj.fast_harmonics(x_me);
            f_ns = T' * f_ns;
            
            f_sec = zeros(3,1);
            % get influence of secondary bodies
            for i=1:length(obj.sec)
                x_13 = obj.sec(i).x(t);
                r_13 = norm(x_13);
                x_s3 = x_13 - x_1s;
                r_s3 = norm(x_s3);
                % Formulation according to Geyling and Westerman for numerical
                % improvements
                % Fundamentals of Astrodynamics, Vallado, pp.575 8-36
                % f_sec = f_sec - obj.sec(i).GM * (x_1s - 3*x_13*(x_1s'*x_13)/r_13^2 - 15/2*((x_1s'*x_13)/r_13^2)^2*x_13)/r_13^3;
                % Standard formulation
                f_sec = f_sec + obj.sec(i).GM * (x_s3/r_s3^3 - x_13/r_13^3);
                % Formulation according to Roy for numerical improvements
                % Q = (r_1s^2 + 2*x_1s'*x_s3) * (r_13^2 + r_13*r_s3 + r_s3^2) / (r_13^3 * r_s3^3 * (r_13 + r_s3));
                % f_sec = f_sec + obj.sec(i).GM * (x_s3 * Q - x_1s/r_s3^3);
            end 
            
            % compute drag
            f_d = 0;
            % compute SRP
            f_SRP = 0;
            if obj.Am > 0
                x_1sun = cspice_spkpos('SUN', t, 'J2000', 'NONE', obj.pri.name);
                x_ssun = x_1sun - x_1s;
                r_ssun = norm(x_ssun);
                x_1e = cspice_spkpos('EARTH', t, 'J2000', 'NONE', obj.pri.name);
                x_se = x_1e - x_1s;
                r_se = norm(x_se);
                x_1m = cspice_spkpos('MOON', t, 'J2000', 'NONE', obj.pri.name);
                x_sm = x_1m - x_1s;
                r_sm = norm(x_sm);

                % determine shadown conditions
                a = asin(obj.Rs/r_ssun);        % apparent radius of sun
                b1 = asin(obj.Re/r_se);         % apparent radius of earth
                b2 = asin(obj.Rm/r_sm);         % apparent radius of moon
                % apparent separation of body centers
                c1 = acos(x_se' * x_ssun / (r_se * r_ssun));
                c2 = acos(x_sm' * x_ssun / (r_sm * r_ssun));

                % handle the earth
                if a + b1 <= c1
                    v1 = 1;
                elseif c1 < abs(b1 - a)
                    v1 = 0;
                else
                    % intermediates
                    x1 = (c1^2 + a^2 - b1^2) / (2*c1);
                    y1 = sqrt(a^2 - x1^2);
                    % area
                    A1 = a^2 * acos(x1/a) + b1^2 * acos((c1-x1)/b1) - c1*y1;
                    % fraction of sunlight left
                    v1 = 1 - A1 / (pi*a^2);
                end
                % handle the moon
                if a + b2 <= c2
                    v2 = 1;
                elseif c2 < abs(b2 - a)
                    v2 = 0;
                else
                    % intermediates
                    x2 = (c2^2 + a^2 - b2^2) / (2*c2);
                    y2 = sqrt(a^2 - x2^2);
                    % area
                    A2 = a^2 * acos(x2/a) + b2^2 * acos((c2-x2)/b2) - c2*y2;
                    % fraction of sunlight left
                    v2 = 1 - A2 / (pi*a^2);
                end

                v = min(v1,v2);     % find maximum occultation

                % assume surface normal points in the direction of the sun
                % (reasonable if solar panels point towards it, which they
                % normally do)
                f_SRP = -v*obj.Ps*obj.Cr*obj.Am*obj.AU^2 * x_ssun/r_ssun^3;
            end
            
            dxdt(1:3) = x(4:6);
            dxdt(4:6) = f_ns + f_sec + f_d + f_SRP;
        end

        function A = partials(obj,t,x)
            %PARTIALS Returns the Jacobian of dynamics w.r.t. x (spherical 
            %harmonics approximated to J2).
            %   Input:
            %    - t; time, seconds past J2000
            %    - x_; state [pos (km); vel (km/s)] of spacecraft
            arguments
                obj     (1,1)   OrbitPropagator
                t       (1,1)   double
                x       (6,1)   double
            end
            
            x_1s = x(1:3);
            r_1s = norm(x_1s);
            J2 = -obj.pri.C(3,1);
            A = zeros(6,6);
            
            A(1:3,4:6) = eye(3);
            % point mass gravity gradient
            A(4:6,1:3) = obj.pri.GM * (3*(x_1s*x_1s')/r_1s^5 - 1/r_1s^3 * eye(3));
            
            % third-body perturbations
            for i=1:length(obj.sec)
                x_si = obj.sec(i).x(t) - x_1s;
                r_si = norm(x_si);
                A(4:6,1:3) = A(4:6,1:3) + obj.sec(i).GM * ...
                             (3*(x_si*x_si')/r_si^5 - 1/r_si^3 * eye(3));
            end
            
            % J2 gravity gradient
            % x is in J2000, transform to Moon ME to apply this partial
            T = cspice_pxform('MOON_ME', 'J2000', t);
            x_1s = T' * x_1s;
            z = x_1s(3);
            S = diag([1 1 3]);
            Q = [1 1 3; 1 1 3; 3 3 3];
            V1 = (5*z^2-r_1s^2)/r_1s^7 * S - 35*z^2*(x_1s*x_1s')/r_1s^9 + ...
                 Q .* (5*(x_1s*x_1s')/r_1s^7);
            % transform the partial back to J2000
            A(4:6,1:3) = A(4:6,1:3) + T * 3*obj.pri.GM*obj.pri.R^2*J2/2 * V1 * T';
        end

        function P = proplyapunov(obj,ts,P0,sv)
            %PROPLYAPUNOV Propagates Lyapunov equations (obj.lyapunov) from
            %given to next time and provides covariance matrices.
            %   Input
            %    - ts; propagation times, seconds past J2000
            %    - P0; covariance of state x0
            %    - sv; (optional) index of space vehicle to propagate for,
            %          if nsats > 1
            arguments
                obj (1,1) OrbitPropagator
                ts  (1,:) double
                P0  (6,6) double
                sv  (1,1) {mustBePositive,mustBeInteger} = 1
            end

            % catch invalid sv input
            if sv ~= 1 && sv > obj.nsats
                error("proplyapunov:invalidInput", ...
                    "sv must be 1 < sv < obj.nsats.");
            end

            ts = ts + obj.t0;

            n = length(ts);
            P = zeros(obj.dim, obj.dim, n);
            
            % reshape starting P to correct format
            P0 = reshape(P0, obj.dim*obj.dim, 1);
            % joint dynamics, since lyapunov requires obj.partials, which
            % requires the current state.
            jointdyn = @(t,x) [obj.dynamics(t, x(1:6)); ...
                obj.lyapunov(x(7:end), obj.partials(t,x(1:6)), 0)];
            [~,Y] = ode89(jointdyn, ts, [obj.x0(:,sv); P0], obj.opts);

            % store covariance matrices in appropriate structure
            for i=1:length(ts)
                P(:,:,i) = reshape(Y(i,7:end)', obj.dim, obj.dim);
            end
        end

        function [fx,C] = modelfit(obj,type,dt,N)
            %MODELFIT Fits given surrogate model type to propagated data.
            %   Input:
            %    - type; "Kepler" or "polynomial", type of model fit -- fit
            %            to error from solving Kepler's problem, or entire
            %            trajectory
            %    - dt; time after t0 to propagate to
            %    - N; number of interpolation points
            %   Output:
            %    - fx; function handle @(t), takes seconds past t0 and returns 
            %          s/c state in J2000
            %    - C; model coefficients
            arguments
                obj     (1,1)   OrbitPropagator
                type    (1,:)   {mustBeText}
                dt      (1,1)   double {mustBePositive}
                N       (1,1)   {mustBeNonnegative,mustBeInteger}
            end
            
            x_init = obj.x0;

            % use chebichev nodes for interpolation
            span = chebichev(N);
            t_interp = (span + 1) * dt / 2;

            % Parse different model types
            if strcmpi(type, "Kepler")
                x_base = obj.keplertool(t_interp, x_init);
                f_base = @(tau) obj.keplertool(tau, x_init);
            elseif strcmpi(type, "polynomial")
                x_base = zeros(6,N+1);
                f_base =  @(tau) zeros(6,length(tau));
            else
                error("modelfit:invalidType", ...
                    "Model type must be either 'Kepler' or 'polynomial'.");
            end

            if N ~= 0
                % propagate state at interpolation points
                [~,x] = obj.runat(t_interp, 'J2000');
                dx = x - x_base;
    
                % compute coefficients for basis and generate model function
                phi = (span').^(0:N);
                B = pinv(phi);
                C = B * dx(1:3,:)';
                basis = @(tau) [tau'.^(0:N) ...
                    (0:N).*(tau'.^([0 0:N - 1])) * 2/dt];
                D = [C zeros(size(C)); zeros(size(C)) C];
    
                % final model function
                fx = @(tau) (basis(2*(tau)/dt - 1) * D)' ...
                            + f_base(tau);
            else
                fx = @(tau) f_base(tau);
            end
        end

        function handles = statetotrajectory(obj)
            %STATETOTRAJECTORY Converts the last propagation to a
            %spline-interpolated trajectory handle.
            arguments
                obj     (1,1)   OrbitPropagator
            end

            data = obj.xs;

            % catch a lack of run happening
            if isempty(obj.frame)
                error("plotlastorbits:noData", ...
                    "No data has been generated yet!");
            end

            % convert trajectory into splines
            handles = cell(obj.nsats, 1);
            for k=1:obj.nsats
                pp = spline(obj.ts, data(:,:,k));
                handles{k} = @(tau,frame) cspice_sxform(obj.frame,frame,tau) ...
                             * ppval(pp, tau);
            end
        end

        function R = rotate_prealloc(obj,t,qs)
            %ROTATE_PREALLOC provides a rotation matrix at time t, given a
            %list of quaternions associated with times t_pre.
            %   This function is used to replace SPICE cspice_*xform calls
            %   from within a propagation so that they can be parallelized.
            %   qs are provided as an extra parameter so that this
            %   function can be used for multiple different
            %   transformations.
            %
            %   Input:
            %    - t; time to get the rotation at (must be w/in bounds of
            %         t_pre)
            %    - qs; list of quaternions corresponding to t_pre
            
            % get surrounding time steps to t
            low = find(obj.t_pre < t, 1, 'last');
            high = find(obj.t_pre >= t, 1);
            % catch if at end or beginning of index
            if isempty(low), low = 1; high = high+1; end
            if isempty(high), high = low; low = low-1; end

            % scale t onto [0,1] and call slerp, then return rotation matrix
            tau = (t - obj.t_pre(low))/(obj.t_pre(high) - obj.t_pre(low));
            q = slerp(qs(low), qs(high), tau);
            R = quat2rotm(q);
        end
    end

    methods (Access = private)
        function f = fast_harmonics(obj,x)
            %FAST_HARMONICS compute gravitational acceleration from spherical
            %harmonics of body
            %   Source: Stellite Orbits - Models, Methods, and
            %   Applications; Montenbruck and Gill; pp. 66-68
            %
            %   Input:
            %    - x; position vector of point (3,) [km]
        
            nn = obj.ord;       % max degree (and order m) of harmonics
            mu = obj.pri.GM;    % km^3/s^2, gravitational parameter of body
            R = obj.pri.R;      % km, reference radius of body
            C = obj.pri.C;      % non-normalized C coefficients of spherical harmonics (n,n)
            S = obj.pri.S;      % non-normalized S coefficients of spherical harmonics (n,n)

            [V, W] = obj.VandW(x, nn+1);
        
            n = 1:nn+1; m = 1:nn+1;
            A = -mu/R^2;
            f= A * [C(1:nn+1,1)' * V(2:nn+2,2)
                    C(1:nn+1,1)' * W(2:nn+2,2)
                    trace((n'-m+1) * (C(n,m).*V(n+1,m) + S(n,m).*W(n+1,m))')];
        
            m = m(2:end);
            f(1) = f(1) + A/2 * trace((C(n,m)*V(n+1,m+1)' + S(n,m)*W(n+1,m+1)') ...
                   - obj.fact * (C(n,m).*V(n+1,m-1) + S(n,m).*W(n+1,m-1))');
            f(2) = f(2) + A/2 * trace((C(n,m)*W(n+1,m+1)' - S(n,m)*V(n+1,m+1)') ...
                   + obj.fact * (C(n,m).*W(n+1,m-1) - S(n,m).*V(n+1,m-1))');
        end

        function [V,W] = VandW(obj,x,nn)
            %VANDW Computes normalized recursive Legendre polynomials for
            %spherical harmonics of body.
            %   Source: Stellite Orbits - Models, Methods, and
            %   Applications; Montenbruck and Gill; pp. 66-68
            %  
            %   Inputs:
            %    - x; position vector of point (3,) [km]
            %    - nn; max degree and order of harmonics

            R = obj.pri.R;
        
            r = norm(x);
            V = zeros(nn+1,nn+1); W = V;
            V(1,1) = R / r; W(1,1) = 0;         % initial conditions
            A = R / r^2;
        
            for m=0:nn          % iterate over m cols (add 1 for indexing)
                if m ~= 0       % if not V_00, W_00 (already defined)
                    V(m+1,m+1) = (2*m-1)*A*(x(1)*V(m,m) - x(2)*W(m,m));
                    W(m+1,m+1) = (2*m-1)*A*(x(1)*W(m,m) + x(2)*V(m,m));
                end
        
                for n=m+1:nn    % iterate over n rows (add 1 for indexing)  
                    V(n+1,m+1) = (2*n-1)/(n-m)*x(3)*A*V(n,m+1);
                    W(n+1,m+1) = (2*n-1)/(n-m)*x(3)*A*W(n,m+1);
                    if n ~= m+1
                        V(n+1,m+1) = V(n+1,m+1) - (n+m-1)/(n-m)*A*R * V(n-1,m+1);
                        W(n+1,m+1) = W(n+1,m+1) - (n+m-1)/(n-m)*A*R * W(n-1,m+1);
                    end
                end
            end

            % handle numbers that are too big (i.e. eliminate degree/order
            % terms that are too small to matter).
            V(or(isinf(V), isnan(V))) = 0;
            W(or(isinf(W), isnan(W))) = 0;
        end

        function x = keplertool(obj,ts,x0)
            %KEPLERTOOL Wrapper for Kepler_universal to handle multiple
            %time requests.
            %   Input:
            %    - ts; states to propagate to by solving Kepler's problem
            arguments
                obj (1,1)   OrbitPropagator
                ts  (1,:)   double
                x0  (6,:)   double
            end

            x = zeros(6,length(ts));
            r0 = x0(1:3);
            v0 = x0(4:6);

            for i=1:length(ts)
                try
                    [rf,vf] = Kepler_universal(r0, v0, ts(i), obj.pri.GM, 1e-10);
                catch
                    0;
                end
                x(:,i) = [rf; vf];
            end
        end
    end
end

