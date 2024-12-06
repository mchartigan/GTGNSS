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
    end
    
    methods
        function obj = OrbitPropagator(t0,ord,varargin)
            %ORBITPROPAGATOR Construct an OrbitPropagator instance.
            %   Inputs:
            %    - t0; character string, 'DD-MMM-YYYY XX:XX:XX'
            %    - ord; maximum degree and order of gravity model to use
            %    - opts; optional name-value arg, ODE45 integration tolerances
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
            default = 2;

            if nargin > default
                for i=1:2:nargin-default
                    if strcmp(varargin{i}, "opts")
                        obj.opts = varargin{i+1};
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
        end
        
        function [ts,xs] = run(obj,tf,n,frame)
            %RUN Propagate the input states for tf seconds (n steps
            %between). Data returned in provided frame.
            %   Input:
            %    - tf; final time, seconds past t0
            %    - n; number of time steps
            %    - frame; reference frame to return data in
            arguments
                obj     (1,1)   OrbitPropagator
                tf      (1,1)   double {mustBePositive}
                n       (1,1)   {mustBeInteger,mustBePositive}
                frame   (1,:)   char
            end

            ts = linspace(0, tf, n);
            [ts,xs] = obj.runat(ts,frame);
        end

        function [ts,xs] = runat(obj,ts,frame)
            %RUNAT Propagate the input states over the provided time steps. 
            %Data returned in provided frame.
            %   Input:
            %    - ts; eval time steps, seconds past t0
            %    - frame; reference frame to return data in
            arguments
                obj     (1,1)   OrbitPropagator
                ts      (1,:)   double
                frame   (1,:)   char
            end

            n = length(ts);
            ts = obj.t0 + ts;
            xs = zeros(6,n,obj.nsats);
            
            for i=1:obj.nsats
                [~,X] = ode45(@obj.dynamics, [obj.t0 ts], obj.x0(:,i), obj.opts);
                X = X(2:end,:)';
                for j=1:length(ts)
                    X(:,j) = cspice_sxform('J2000', frame, ts(j)) * X(:,j);
                end
                xs(:,:,i) = X;
            end

            obj.ts = ts;
            obj.xs = xs;
            obj.frame = frame;
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
            
            % get lunar nonspherical gravity effects
            T = cspice_pxform('J2000', obj.pri.frame, t);
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
            %given to next time and provides covariance matrices. This
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

            n = length(ts);
            P = zeros(obj.dim, obj.dim, n);
            
            % reshape starting P to correct format
            P0 = reshape(P0, obj.dim*obj.dim, 1);
            % joint dynamics, since lyapunov requires obj.partials, which
            % requires the current state.
            jointdyn = @(t,x) [obj.dynamics(t, x(1:6)); ...
                obj.lyapunov(x(7:end), obj.partials(t,x(1:6)), 0)];
            [~,Y] = ode45(jointdyn, ts, [obj.x0(:,sv); P0], obj.opts);

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
            %   Input:
            %    - frame; reference frame for trajectories
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
            arguments
                obj (1,1)   OrbitPropagator
                x   (3,1)   double
            end
        
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
            fact = factorial(abs(n'-m+2))./factorial(abs(n'-m));
            f(1) = f(1) + A/2 * trace((C(n,m)*V(n+1,m+1)' + S(n,m)*W(n+1,m+1)') ...
                   - fact * (C(n,m).*V(n+1,m-1) + S(n,m).*W(n+1,m-1))');
            f(2) = f(2) + A/2 * trace((C(n,m)*W(n+1,m+1)' - S(n,m)*V(n+1,m+1)') ...
                   + fact * (C(n,m).*W(n+1,m-1) - S(n,m).*V(n+1,m-1))');
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

