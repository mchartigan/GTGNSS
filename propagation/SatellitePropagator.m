classdef SatellitePropagator < Propagator
    %SATELLITEPROPAGATOR Propagator for a satellite, which tracks its
    %position, velocity, and time (PVT) information.

    properties
        orbit   (1,1)   OrbitPropagator = OrbitPropagator(1)
        clock   (1,1)   Clock = Clock("none", zeros(4,1))
        % options for IMU and fixed-step integration ('flight' mode)
        imu     (1,1)   IMU
        flight  (1,1)   {mustBeNonnegative,mustBeInteger} = 0
        step    (1,1)   double {mustBePositive} = 0.025
        % dimension of state
        dim     = 9
        % optional parameter to set default frame of dynamics
        frame   (1,:)   {mustBeText} = 'J2000'
        % optional scaling of process noise
        scale   (1,1)   double = 1
        % measurement bias and noise to add to propagator
        bias    (1,:)   Propagator = RandomRun.empty
        nbias   (1,1)   {mustBeInteger,mustBeNonnegative} = 0
    end

    methods
        function obj = SatellitePropagator(orbit,clock,imu,options)
            %SATELLITEPROPAGATOR Construct a Satellite instance.
            %   Input:
            %    - orbit; orbit propagator with appropriate fidelity
            %    - clock; clock propagator representative of what's onboard
            %    - imu; IMU instance, if satellite is inertially aided
            arguments
                orbit           (1,1)   OrbitPropagator
                clock           (1,1)   Clock
                imu             (1,1)   IMU = IMU()
                options.flight  (1,1)   = 0
                options.tol     (1,1)   double {mustBePositive} = 1e-5
                options.bias    (1,:)   Propagator = RandomRun.empty
            end
            
            obj.orbit = orbit;
            obj.clock = clock;
            obj.imu = imu;
            obj.flight = options.flight;

            if obj.flight
                % compute RK4 step size. Since alg error is O(h^5), solve for h
                % and then round down to nearest size that evenly divides 0.1.
                % This effectively makes the minimum tolerance 1e-5.
                h = options.tol^(1/5);
                obj.step = 0.1 / ceil(0.1 / h);
                obj.imu.step = obj.step;
            end

            if ~isempty(options.bias)
                % store bias models
                obj.bias = options.bias;
                obj.nbias = length(obj.bias);

                % add new dimensions
                for i=1:obj.nbias
                    obj.dim = obj.dim + obj.bias(i).dim;
                end
            end
        end

        function [ts,xs] = run(obj,ts,x0,n,frame,noise)
            %RUN Propagate the provided state(s) for ts(2)-ts(1) seconds
            %(n steps between). Data returned in provided frame.
            %   Input:
            %    - ts; [intial time, final time], seconds past J2000
            %    - x0 (6,:) double; starting state(s) in frame
            %    - n; number of time steps
            %    - frame; reference frame of x0; data is also returned in
            %       this frame
            arguments
                obj     (1,1)   SatellitePropagator
                ts      (1,2)   double
                x0      (:,:)   double
                n       (1,1)   {mustBeInteger,mustBePositive}
                frame   (1,:)   char = obj.frame
                noise   (1,1)   = false
            end

            ts = linspace(ts(1), ts(2), n);
            xs = obj.runat(ts,x0,frame,noise);
        end

        function xs = runat(obj,ts,x0,frame,noise)
            %RUNAT Propagate the provided states over the provided time steps. 
            %Data returned in indicated frame.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - x0; starting states in frame
            %    - frame; reference frame of x0; data is also returned in
            %       this frame
            arguments
                obj     (1,1)   SatellitePropagator
                ts      (1,:)   double
                x0      (:,:)   double
                frame   (1,:)   char = obj.frame
                noise   (1,1)   = false
            end

            xs = zeros(obj.dim, length(ts));
            if ~obj.flight
                xs(1:6,:) = obj.orbit.runat(ts,x0(1:6),frame);
            else
                % fixed-step integration for IMU purposes
                n = length(ts);
                xs(1:6,1) = x0(1:6);
                for i=2:n
                    xs(1:6,i) = obj.RK4(ts(i-1:i), xs(1:6,i-1));
                end
            end
            xs(7:9,:) = obj.clock.runat(ts,x0(7:9),noise);

            k = 9;
            for i=1:obj.nbias
                ind = k+1:k+obj.bias(i).dim;
                xs(ind,:) = obj.bias(i).runat(ts,x0(ind),noise);
                k = k + obj.bias(i).dim;
            end
        end
        
        function xf = RK4(obj,ts,x0)
            %RK4 Fixed-step numerical integrator to provide xs, given ts
            %and x0. Functions best when step sizes ts(k)-ts(k-1) are
            %multiples of obj.h (or 0.1s).
            %   Input:
            %    - ts; eval time steps [t0 tf], seconds past J2000
            %    - x0; starting states in frame

            teval = ts(1):obj.step:ts(end);
            if teval(end) ~= ts(end), teval = [teval ts(end)]; end

            % get dynamics around equilibrium (t0). Using 1st order Taylor
            % expansion to minimize obj.dynamics calls.
            g = obj.orbit.dynamics(ts(1), x0);
            G = obj.orbit.numpart(ts(1), x0);
            xf = x0;
            for i=1:length(teval)-1
                % time step
                ti = teval(i);
                dt = teval(i+1) - ti;

                % compute coefficients w/ Taylor expansion
                C = g - G*x0;
                A1 = [zeros(3,1); obj.imu.read(ti+dt/2)];
                F0 = [zeros(3,1); obj.imu.read(ti)] + G*xf + C;
                F1 = A1 + G*(xf + dt/2*F0) + C;
                F2 = A1 + G*(xf + dt/2*F1) + C;
                F3 = [zeros(3,1); obj.imu.read(ti+dt)] + G*(xf + dt*F2) + C;

                % compute Runge-Kutta update
                xf = xf + (dt/6) * (F0 + 2*F1 + 2*F2 + F3);
            end
        end

        function dxdt = dynamics(obj,t,x)
            %DYNAMICS Inertial dynamics for an orbiting satellite carrying
            %a clock.
            %   Input:
            %    - t; time, seconds past J2000
            %    - x (9,1) double; state

            dxdt = zeros(obj.dim,1);
            dxdt(1:6) = obj.orbit.dynamics(t,x(1:6));
            dxdt(4:6) = dxdt(4:6) + obj.imu.read(t);
            dxdt(7:9) = obj.clock.dynamics(t,x(7:9));
            
            k = 9;
            for i=1:obj.nbias
                ind = k+1:k+obj.bias(i).dim;
                dxdt(ind) = obj.bias(i).dynamics(t,x(ind));
                k = k + obj.bias(i).dim;
            end
        end

        function A = partials(obj,t,x)
            %PARTIALS Jacobian of satellite and clock dynamics w.r.t. the state.
            A = zeros(obj.dim, obj.dim);
            A(1:6,1:6) = obj.orbit.numpart(t,x(1:6));
            A(7:9,7:9) = obj.clock.partials(t,x(7:9));
            
            k = 9;
            for i=1:obj.nbias
                ind = k+1:k+obj.bias(i).dim;
                A(ind,ind) = obj.bias(i).partials();
                k = k + obj.bias(i).dim;
            end
        end

        function P = proplyapunov(obj,ts,x0,P0)
            %PROPLYAPUNOV Propagates Lyapunov equations between times and
            %provides covariance matrices.
            %
            %   Input:
            %    - ts (1,:) double; propagation times, seconds past J2000
            %    - x0 (6,1) double; starting state in J2000
            %    - P0 (6,6) double; covariance of state x0
            arguments
                obj (1,1)   SatellitePropagator
                ts  (1,:)   double
                x0  (9,1)   double
                P0  (9,9)   double
            end

            n = length(ts);
            P = zeros(obj.dim, obj.dim,n);
            P(1:6,1:6,:) = obj.orbit.proplyapunov(ts, x0(1:6), P0(1:6,1:6));
            P(7:9,7:9,:) = obj.clock.proplyapunov(ts, P0(7:9,7:9));
        end

        function Q = noise(obj,dt,x)
            %PNC Returns the dicrete-time process noise covariance of the
            %propagator over a time step dt.
            %   Input:
            %    - dt; time step, in s
            %    - x; state at starting time

            Q = zeros(obj.dim,obj.dim);
            Q(1:6,1:6) = obj.orbit.noise(dt,x(1:6)) + obj.imu.noise(dt);
            Q(7:9,7:9) = obj.clock.noise(dt);

            k = 9;
            for i=1:obj.nbias
                ind = k+1:k+obj.bias(i).dim;
                Q(ind,ind) = obj.bias(i).noise(dt);
                k = k + obj.bias(i).dim;
            end

            Q = Q * obj.scale;
        end

        function [fx,C] = modelfit(obj,trajs,type,N)
            %MODELFIT Fits chosen surrogate model type to trajectories of
            %satellite and clock.
            %   Input:
            %    - traj (2,1) Trajectory; satellite states and clock
            %       states, sharing t0
            %    - type; "Kepler" or "polynomial", type of model fit -- fit
            %       to error from solving Kepler's problem, or entire
            %       trajectory
            %    - N; number of interpolation points
            %   Output:
            %    - fx; function handle @(t), takes seconds past t0 and returns 
            %       s/c state in J2000 and clock state
            %    - C; model coefficients

            [fx1,C1] = obj.orbit.modelfit(trajs(1), type, N);
            [fx2,C2] = obj.clock.modelfit(trajs(2));

            fx = @(dt) [fx1(dt); fx2(dt)];
            C = [C1; C2];
        end
    end
end