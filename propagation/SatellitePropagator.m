classdef SatellitePropagator < Propagator
    %SATELLITEPROPAGATOR Propagator for a satellite, which tracks its
    %position, velocity, and time (PVT) information.

    properties
        orbit   (1,1)   OrbitPropagator = OrbitPropagator(1)
        clock   (1,1)   Clock = Clock("none", zeros(4,1))
        dim     = 9         % dimension of state
    end

    methods
        function obj = SatellitePropagator(orbit,clock)
            %SATELLITEPROPAGATOR Construct a Satellite instance.
            %   Input:
            %    - orbit; orbit propagator with appropriate fidelity
            %    - clock; clock propagator representative of what's onboard

            if nargin ~= 0
                obj.orbit = orbit;
                obj.clock = clock;
            end
        end

        function [ts,xs] = run(obj,ts,x0,n,frame)
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
                ts      (1,:)   double
                x0      (9,:)   double
                n       (1,1)   {mustBeInteger,mustBePositive}
                frame   (1,:)   char
            end

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
            arguments
                obj     (1,1)   SatellitePropagator
                ts      (1,:)   double
                x0      (9,:)   double
                frame   (1,:)   char
            end

            xs = zeros(obj.dim, length(ts));
            xs(1:6,:) = obj.orbit.runat(ts,x0(1:6),frame);
            xs(7:9,:) = obj.clock.runat(ts,x0(7:9));
        end

        function dxdt = dynamics(obj,t,x)
            %DYNAMICS Inertial dynamics for an orbiting satellite carrying
            %a clock.
            %   Input:
            %    - t; time, seconds past J2000
            %    - x (9,1) double; state

            dxdt = [obj.orbit.dynamics(t,x(1:6)); obj.clock.dynamics(t,x(7:9))];
        end

        function A = partials(obj,t,x)
            %PARTIALS Jacobian of satellite and clock dynamics w.r.t. the state.
            A = zeros(obj.dim, obj.dim);
            A(1:6,1:6) = obj.orbit.partials(t,x(1:6));
            A(7:9,7:9) = obj.clock.partials(t,x(7:9));
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