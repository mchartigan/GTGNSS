classdef Clock < Propagator
    %CLOCK Generic propagator for an oscillator/clock drift over time.
    %Provides various default oscillator options and implements all common
    %properties/methods of the abstract parent class, Propagator (also
    %shared by OrbitPropagator).
    %
    %   For additional reading, see Zucca and Tavella, 2005: 
    %   https://www.doi.org/10.1109/TUFFC.2005.1406554
    
    properties
        t0                              % sim start time, seconds past J2000
        x0                              % sim start state
        var     (3,3)   double          % variance of clock noise
        dim             = 3             % dimension of state
        % info stored from latest run (if ts = [], no run yet)
        ts
        xs
    end
    
    methods
        function obj = Clock(t0,x0,type)
            %CLOCK Construct an oscillator/clock instance.
            %   Inputs:
            %    - t0; starting time in seconds past J2000
            %    - x0; starting state, [bias (s); drift (s/s); aging (1/s)]
            %    - type; options include "CSAC", "RAFS", or "USO"
            arguments
                t0      (1,1)   double
                x0      (3,1)   double
                type    (1,:)
            end
            
            obj.t0 = t0;
            obj.x0 = x0;
            % parse clock type
            if isa(type, 'string')
                if strcmp(type, "USO")          % ultra-stable oscillator
                    [s1,s2,s3] = DiffCoeffUSO();
                elseif strcmp(type, "CSAC")     % chip-scale atomic clock
                    [s1,s2,s3] = DiffCoeffCSAC();
                elseif strcmp(type, "RAFS")     % Rubidium atomic frequency standard
                    [s1,s2,s3] = DiffCoeffRAFS();
                elseif strcmpi(type, "none")    % no clock (default initialization)
                    return
                else
                    error("Clock:invalidType", ...
                        "Invalid clock type of %s. See documentation.", type);
                end

                obj.var = diag([s1 s2 s3].^2);

            elseif all(size(type) == [1 3]) && isa(type, 'double')
                obj.var = diag(type .^ 2);
            else
                error("Clock:invalidType", ...
                    "Clock type must either be string or list of coefficients.");
            end
        end

        function [ts,xs,vs] = run(obj,tf,n)
            %RUN Propagate the input states for tf seconds (n steps
            %between).
            %   Input:
            %    - tf; final time, seconds past t0
            %    - n; number of time steps
            %   Output:
            %    - ts; times in seconds past J2000
            %    - xs; clock states at ts
            %    - vs; clock covariance at ts
            arguments
                obj (1,1)   Clock
                tf  (1,1)   double {mustBePositive}
                n   (1,1)   {mustBeInteger,mustBePositive}
            end

            % initialize variables
            ts = linspace(0,tf,n);
            [ts,xs,vs] = obj.runat(ts);
        end

        function [ts,xs,vs] = runat(obj,ts)
            %RUNAT Propagate the input states over the provided time steps.
            %   Input:
            %    - ts; eval time steps, seconds past t0
            %   Output:
            %    - ts; times in seconds past J2000
            %    - xs; clock states at ts
            %    - vs; clock covariance at ts
            arguments
                obj (1,1)   Clock
                ts  (1,:)   double {mustBeNonnegative}
            end

            % initialize variables
            n = length(ts);
            xs = zeros(3,n);
            xs(:,1) = obj.x0;
            vs = zeros(3,3,n);

            for i=2:n
                dt = ts(i) - ts(i-1);
                stm = Clock.STM(dt);
            
                % innovation vector, J ~ N(0,Q)
                J = mvnrnd([0 0 0], obj.dtcovariance(dt), 1)';
                % covariance of x at ts(i)
                vs(:,:,i) = obj.dtcovariance(ts(i)); 
                xs(:,i) = stm * xs(:,i-1) + J;
            end

            obj.ts = ts + obj.t0;
        end

        function P = proplyapunov(obj,ts,P0)
            %PROPLYAPUNOV Propagates Lyapunov equations (obj.lyapunov) from
            %given to next time and provides covariance matrices. This
            %   Input
            %    - ts; propagation times, seconds past J2000
            %    - P0; covariance of state x0
            arguments
                obj (1,1) Clock
                ts  (1,:) double
                P0  (3,3) double
            end

            n = length(ts);
            P = zeros(obj.dim, obj.dim, n);
            
            % reshape starting P to correct format
            P0 = reshape(P0, obj.dim*obj.dim, 1);
            [~,Y] = ode45(@(t,p) obj.lyapunov(p, Clock.partials(t,obj.x0), obj.var), ...
                          ts, P0, odeset("RelTol", 1e-9, "AbsTol", 1e-11));

            % store covariance matrices in appropriate structure
            for i=1:length(ts)
                P(:,:,i) = reshape(Y(i,:)', obj.dim, obj.dim);
            end
        end

        function Q = dtcovariance(obj,dt)
            %DTCOVARIANCE Returns the discrete-time covariance associated
            %w/ the Wiener processes (obj.var) is the continuous-time
            %coviarance.
            %   Input:
            %    - dt; time step
            arguments
                obj (1,1)   Clock
                dt  (1,1)   double {mustBePositive}
            end

            v1 = obj.var(1,1);
            v2 = obj.var(2,2);
            v3 = obj.var(3,3);

            Q = [v1*dt + v2/3*dt^3 + v3/20*dt^5 v2/2*dt^2 + v3/8*dt^4 v3/6*dt^3;
                 v2/2*dt^2 + v3/8*dt^4          v2*dt + v3/3*dt^3     v3/2*dt^2;
                 v3/6*dt^3                      v3/2*dt^2             v3*dt   ];
        end

        function [fx,C] = modelfit(obj)
            %MODELFIT Returns a second-order polynomial model for the clock
            %state over time. Starting epoch is the current t0, x0
            %   Output:
            %    - fx; @(t) function handle, input seconds past J2000 and
            %          it returns clock state
            %    - C; current bias, drift, aging used
            arguments
                obj (1,1)   Clock
            end

            C = obj.x0;
            fx = @(t) [obj.x0(1) + (t-obj.t0)*obj.x0(2) + (t-obj.t0).^2/2*obj.x0(3); ...
                       obj.x0(2) + (t-obj.t0)*obj.x0(3); ...
                       ones(size(t))*obj.x0(3)];
        end
    end

    methods (Static)
        function stm = STM(dt)
            %STM Returns the DT state transition matrix based on the dynamics
            %defined in Zucca and Tavella.
            %   Input:
            %    - dt; time step
            arguments
                dt  (1,1)   double {mustBePositive}
            end

            stm = [1 dt dt^2/2; 0 1 dt; 0 0 1];
        end

        function dxdt = dynamics(~,x)
            %DYNAMICS Invokes the clock dynamics based on the Zucca and
            %Tavella paper.
            %   Input:
            %    - t; simulation time, not used
            %    - x; current clock state
            arguments
                ~
                x   (3,1)   double
            end

            dxdt = [0 1 0; 0 0 1; 0 0 0] * x;
        end

        function A = partials(~,~)
            %PARTIALS Invokes the partials of dynamics.
            %   Input:
            %    - t; simulation time, not used
            %    - x; current clock state, also not used lol

            A = [0 1 0; 0 0 1; 0 0 0];
        end
    end
end

