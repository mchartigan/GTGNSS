classdef RandomRun < Propagator
    %RANDOMRUN Random run model, where the rate of change is a random walk
    %process.

    properties
        % continuous-time variance of rate diffusion coefficient
        var (1,1)   double {mustBeNonnegative} = 0
        % dimension of state
        dim = 2
        % scale factor
        a   (1,1)   double = 1
    end

    methods
        function obj = RandomRun(var)
            %RANDOMRUN Construct a RandomRun instance
            %   Input:
            %    - variance of random run process

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
                obj     (1,1)   RandomRun
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
                obj     (1,1)   RandomRun
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

        function dxdt = dynamics(obj,~,x)
            %DYNAMICS Returns the derivative of the state
            %   Input:
            %    - t; simulation time
            %    - x; state at time t
            dxdt = obj.partials() * x;
        end

        function Q = noise(obj,dt)
            %NOISE Returns the discrete-time process noise covariance over
            %the interval dt.
            %   Input:
            %    - dt; time interval, in s
            Q = obj.var * [obj.a^2*dt^3/3 obj.a*dt^2/2; obj.a*dt^2/2 dt];
        end

        function A = partials(obj)
            %PARTIALS Returns the jacobian of the state
            A = [0 obj.a; 0 0];
        end

        function Phi = STM(obj,dt)
            %STM Returns the state transition matrix over a time step dt
            Phi = [1 obj.a*dt; 0 1];
        end
    end
end