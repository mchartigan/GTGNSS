classdef IMUPropagator < Propagator
    %IMUPROPAGATOR Propagates body dynamics with inertial aiding through an
    %Inertial Measurement Unit (IMU).
    %   Methods adapted from "Navigation Filtering Best Practices" and the
    %   GEONS Math Spec.

    properties
        % Inertial (unforced) dynamics of body
        prop    (1,1)   Propagator = OrbitPropagator(1)
        opts    (1,1)   struct
        step    (1,1)   double {mustBePositive} = 0.025
        % dimension of state
        dim     = 6
        % Trajectory instance of s/c acceleration over time
        a_m     (1,1)   Trajectory
    end

    methods
        function obj = IMUPropagator(prop,traj,options)
            %IMUPROPAGATOR Construct an IMUPropagator instance.
            %   Input:
            %    - prop; Propagator instance that defines the gravitational
            %       acceleration on the body.
            %    - opts; optional name-value arg, options for ODE45
            %       integrator.
            arguments
                prop            (1,1)   Propagator
                traj            (1,1)   Trajectory
                options.tol     (1,1)   double {mustBePositive} = 1e-5
                options.opts    (1,1)   struct = odeset(RelTol=1e-8, AbsTol=1e-9)
            end

            obj.prop = prop;
            obj.dim = prop.dim;
            obj.opts = options.opts;
            obj.a_m = traj;

            % compute RK4 step size. Since alg error is O(h^5), solve for h
            % and then round down to nearest size that evenly divides 0.1.
            % This effectively makes the minimum tolerance 1e-5.
            h = options.tol^(1/5);
            obj.step = 0.1 / ceil(0.1 / h);
        end

        function [ts,xs] = run(obj,ts,x0,n,frame)
            %RUN Propagate the provided states over the provided time steps. 
            %Data returned in indicated frame.
            %   Input:
            %    - ts; [intial time, final time], seconds past J2000
            %    - x0 (6,:) double; starting state(s) in frame
            %    - n; number of time steps
            %    - frame; reference frame of x0; data is also returned in
            %       this frame
            
            ts = linspace(ts(1), ts(2), n);
            xs = obj.runat(ts,x0,frame);
        end

        function xs = runat(obj,ts,x0,frame)
            %RUNAT Propagate the provided states over the provided time steps. 
            %Data returned in indicated frame.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - x0; starting states in frame
            %    - frame; reference frame of x0; data is also returned in
            %       this frame

            if ts(1) < obj.a_m.ts(1) || ts(end) > obj.a_m.ts(end)
                error("IMUPropagator:timeOutOfBounds", ...
                    "ts exceeds provided IMU measurements.");
            end

            n = length(ts);

            % switch to inertial frame, since that's where we'll be operating
            % x0 = cspice_sxform(frame, 'J2000', ts(1)) * x0;
            xs = zeros(obj.dim, n);
            xs(:,1) = cspice_sxform(frame, 'J2000', ts(1)) * x0;

            % % integrate
            % [~,X] = ode45(@obj.dynamics, ts, x0, obj.opts);
            % if n == 1
            %     xs = X(end,:)';
            % elseif n == 2
            %     xs = [X(1,:)' X(end,:)'];
            % else
            %     xs = X';
            % end
            % fixed-step integration for IMU purposes
            for i=2:n
                xs(:,i) = obj.RK4(ts(i-1:i), xs(:,i-1));
            end

            if ~strcmp(frame, 'J2000')  % transform back if necessary
                for j=1:length(ts)
                    xs(:,j) = cspice_sxform('J2000', frame, ts(j)) * xs(:,j);
                end
            end
        end

        function a = readIMU(obj,t)
            %READIMU Returns sensed acceleration (in ICRF frame) at time t.
            %   Input:
            %    - t; time in s past J2000

            a = obj.a_m.get(t);
        end

        function dxdt = dynamics(obj,t,x)
            %DYNAMICS Inertial dynamics for body, based on IMU and provided
            %gravity model.
            %   Input:
            %    - t (1,1) double; simulation time in seconds past J2000
            %    - x (6,1) double; state [pos (km); vel (km/s)] of
            %        body in J2000

            dxdt = obj.prop.dynamics(t,x);
            % dxdt(4:6) = dxdt(4:6) + obj.readIMU(t);
        end

        function out = modelfit(obj,traj)
            out = 0;
        end

        function A = partials(obj,t,x)
            %NUMPART Returns the CT dynamics matrix at time t, computed
            %using central differences.
            %   Input:
            %    - t (1,1) double; simulation time in seconds past J2000
            %    - x (6,1) double; state [pos (km); vel (km/s)] of
            %        body in J2000

            A = obj.prop.numpart(t,x);
        end

        function P = proplyapunov(obj,ts,P0)
            P = zeros(obj.dim, obj.dim);
        end

        function Q = noise(obj,dt,x)
            Q = zeros(obj.dim, obj.dim);
        end
    end
end