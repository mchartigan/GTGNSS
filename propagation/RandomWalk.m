classdef RandomWalk < Propagator
    %RANDOMWALK Random walk model, where the state is a random walk
    %process.

    properties
        % continuous-time variance of rate diffusion coefficient
        var (1,1)   double {mustBeNonnegative} = 0
        % dimension of state
        dim = 1
    end

    methods
        function obj = RandomWalk(var)
            %RANDOMWALK Construct a RandomWalk instance
            %   Input:
            %    - variance of random walk process

            obj.var = var;
        end

        function [ts,xs] = run(obj,ts,x0,n,noise)
            %RUN Propagate the input states for tf seconds (n steps
            %between).
            %   Input:
            %    - ts; [intial time, final time], seconds past J2000
            %    - x0 (3,1) double; starting state
            %    - n; number of time steps
            %    - noise; optional boolean, include random noise
            %   Output:
            %    - ts; times in seconds past J2000
            %    - xs; states at ts
            arguments
                obj     (1,1)   RandomWalk
                ts      (1,2)   double {mustBePositive}
                x0      (:,1)   double
                n       (1,1)   {mustBeInteger,mustBePositive}
                noise   (1,1)   = false
            end

            % initialize variables
            ts = linspace(ts(1),ts(end),n);
            xs = obj.runat(ts,x0,noise);
        end

        function xs = runat(obj,ts,x0,noise)
            %RUNAT Propagate the input states over the provided time steps.
            %   Input:
            %    - ts; eval time steps, seconds past t0
            %   Output:
            %    - ts; times in seconds past J2000
            %    - xs; clock states at ts
            %    - vs; clock covariance at ts
            arguments
                obj     (1,1)   RandomWalk
                ts      (1,:)   double {mustBeNonnegative}
                x0      (:,1)   double
                noise   (1,1)   = false
            end

            % initialize variables
            n = length(ts);
            xs = zeros(obj.dim,n);
            xs(:,1) = x0;

            for i=2:n
                dt = ts(i) - ts(i-1);
                stm = obj.STM(dt);
            
                xs(:,i) = stm * xs(:,i-1);
                if noise
                    % innovation vector, J ~ N(0,Q)
                    J = mvnrnd(zeros(1,obj.dim), obj.noise(dt), 1)';
                    xs(:,i) = xs(:,i) + J;
                end
            end
        end

        function Q = noise(obj,dt)
            %NOISE Returns the discrete-time process noise covariance over
            %the interval dt.
            %   Input:
            %    - dt; time interval, in s
            Q = obj.var * dt;
        end
    end

    methods (Static)
        function dxdt = dynamics(~,~)
            %DYNAMICS Returns the derivative of the state
            %   Input:
            %    - t; simulation time
            %    - x; state at time t
            dxdt = 0;
        end

        function A = partials()
            %PARTIALS Returns the jacobian of the state
            A = 0;
        end

        function Phi = STM(~)
            %STM Returns the state transition matrix over a time step dt
            Phi = 1;
        end
    end
end