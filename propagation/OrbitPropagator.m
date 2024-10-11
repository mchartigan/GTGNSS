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
                ts      (1,:)   double {mustBeNonnegative}
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
            %DYNAMICS Invokes the orbital dynamics, with relevant settings,
            %for this propagator instance.
            %   Input:
            %    - t; simulation time in seconds past J2000
            %    - x; satellite state
            arguments
                obj     (1,1)   OrbitPropagator
                t       (1,1)   double
                x       (6,1)   double
            end

            dxdt = orbitaldynamics(t, x, obj.pri, obj.ord, obj.sec);
        end

        function A = partials(obj,t,x)
            %PARTIALS Invokes the partials of dynamics, with relevant
            %settings, for this propagator instance.
            %   Input:
            %    - t; simulation time in seconds past J2000
            %    - x; satellite state
            arguments
                obj     (1,1)   OrbitPropagator
                t       (1,1)   double
                x       (6,1)   double
            end

            A = orbitalpartials(t, x, obj.pri, obj.sec);
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
            %    - fx; function handle @(t), returns s/c state in J2000
            %    - C; model coefficients
            arguments
                obj     (1,1)   OrbitPropagator
                type    (1,:)   {mustBeText}
                dt      (1,1)   double {mustBePositive}
                N       (1,1)   {mustBeNonnegative,mustBeInteger}
            end

            % basis = @(tau) cell2mat(arrayfun(@(x) chebyshevT(0:N_APPX,x), tau, 'UniformOutput', false));

            % use chebichev nodes for interpolation
            span = chebichev(N);
            t_interp = (span + 1) * dt / 2;

            % Parse different model types
            if strcmpi(type, "Kepler")
                x_base = obj.keplertool(t_interp);
                f_base = @(tau) obj.keplertool(tau - obj.t0);
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
                fx = @(tau) (basis(2*(tau - obj.t0)/dt - 1) * D)' ...
                            + f_base(tau);
            else
                fx = @(tau) f_base(tau);
            end
        end

        function x = keplertool(obj,ts)
            %KEPLERTOOL Wrapper for Kepler_universal to handle multiple
            %time requests.
            %   Input:
            %    - ts; states to propagate to by solving Kepler's problem
            arguments
                obj (1,1)   OrbitPropagator
                ts  (1,:)   double {mustBeNonnegative}
            end

            x = zeros(6,length(ts));
            r0 = obj.x0(1:3);
            v0 = obj.x0(4:6);

            for i=1:length(ts)
                try
                    [rf,vf] = Kepler_universal(r0, v0, ts(i), obj.pri.GM, 1e-10);
                catch
                    0;
                end
                x(:,i) = [rf; vf];
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
end

