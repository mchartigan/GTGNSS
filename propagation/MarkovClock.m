classdef MarkovClock < Propagator
    %MARKOVCLOCK Markov process-based clock model for propagation and
    %filtering.
    %   NOTE: Everything is normalized (multiplied) by c by default
    
    properties
        % starting time of simulation, seconds past J2000
        t0
        % starting state of simulation
        x0
        % size of state
        dim
        % no. of Markov processes
        m           (1,1)   {mustBePositive,mustBeInteger} = 1
        % propagation seed so noise is consistent between runs
        seed        (1,1)   {mustBeInteger} = 0
        % info stored from latest run (if ts = [], no run yet)
        ts
        xs

        % FIT UNCERTAINTY PARAMETERS %
        % white frequency noise
        sigma_w     (1,1)   double {mustBeNonnegative} = 0
        % random walk frequency noise
        sigma_rw    (1,1)   double {mustBeNonnegative} = 0
        % frequency drift (aging rate) uncertainty
        sigma_a     (1,1)   double {mustBeNonnegative} = 0
        % Markov process noise
        sigma_m     (:,1)   double {mustBeNonnegative} = 0
        % Markov process time constants
        R           (:,1)   double {mustBePositive} = 1

        % CLOCK INFO FROM JSON %
        % Hz, oscillator frequency
        f           (1,1)   double {mustBePositive} = 1
        % s/s/s, aging rate
        a           (1,1)   double = 0
        % s, stability intervals
        t_stab      (:,1)   double {mustBePositive} = []
        % s/s, stability (standard deviations)
        s_stab      (:,1)   double {mustBePositive} = []
        % Hz, phase noise frequency offsets
        f_noise     (:,1)   double {mustBePositive} = []
        % dBc/Hz, phase noise
        n_noise     (:,1)   double = []
    end
    properties (Constant)
        % m/s, speed of light
        c = 299792458;
    end
    
    methods
        function obj = MarkovClock(t0,name,minout)
            %MARKOVCLOCK Construct an instance of MarkovClock
            %   Inputs:
            %    - t0; starting time in seconds past J2000
            %    - name; filename (no extension) of clock info .json
            %         available on path
            %    - minout; output struct from MarkovOpt that dictates the
            %         stability

            % NOTE: Uncertainties (standard deviations) and states are
            % normalized by c throughout.
            
            % assign uncertainties
            obj.sigma_w  = minout.x(1) * obj.c;
            obj.sigma_rw = minout.x(2) * obj.c;
            % if output is odd, it means aging uncertainty was fit
            fit_a = mod(length(minout.x), 2);
            obj.sigma_a = minout.x(3) * fit_a * obj.c;
            % assign set number of Markov processes
            obj.m = (length(minout.x) - 2 - fit_a)/2;
            obj.sigma_m = zeros(obj.m, 1);
            obj.R       = ones(obj.m, 1);
            for i=1:2:2*obj.m
                obj.sigma_m(i) = minout.x(2+fit_a+i) * obj.c;
                obj.R(i)       = minout.x(3+fit_a+i);
            end

            % set starting time
            obj.t0 = t0;
            % set state size
            obj.dim = 3 + obj.m;
            % set static seed
            obj.seed = randi([1 1e9]);
            % get info from datasheet .json
            obj.assignclockdata(name);
            obj.a = obj.a * obj.c;
            obj.x0 = zeros(obj.dim,1);
            obj.x0(3) = obj.a;
        end

        function [ts,xs] = run(obj,x0,tf,n)
            %RUN Propagate the input states for tf seconds (n steps
            %between).
            %   Input:
            %    - x0; starting clock state (phase offset, frequency offset,
            %         frequency drift)
            %    - tf; final time, seconds past t0
            %    - n; number of time steps
            %   Output:
            %    - ts; times in seconds past J2000
            %    - xs; clock states at ts
            arguments
                obj (1,1)   MarkovClock
                x0  (3,1)   double
                tf  (1,1)   double {mustBePositive}
                n   (1,1)   {mustBeInteger,mustBePositive}
            end

            % initialize variables
            ts = linspace(0,tf,n);
            [ts,xs] = obj.runat(x0,ts);
        end
        
        function [ts,xs] = runat(obj,x0,ts)
            %RUNAT Propagate the input states over the provided time steps.
            %   Input:
            %    - x0; starting clock state (phase offset, frequency offset,
            %         frequency drift)
            %    - ts; eval time steps, seconds past t0
            %   Output:
            %    - ts; times in seconds past J2000
            %    - xs; clock states at ts
            arguments
                obj (1,1)   MarkovClock
                x0  (3,1)   double
                ts  (1,:)   double {mustBeNonnegative}
            end
            
            % rng(obj.seed)       % initialize rng for consistency

            % initialize variables
            n = length(ts);
            xs = zeros(3,n);
            xs(1:3,1) = x0;
            % set starting state of Markov processes as RV with mean 0 and
            % variance U = sigma_m^2/(2*R)
            for i=1:obj.m
                xs(3+i,1) = mvnrnd(0, obj.sigma_m(i)^2/(2*obj.R(i)));
            end
            % assign starting state for obj.modelfit()
            obj.x0 = xs(:,1);

            for i=2:n
                dt = ts(i) - ts(i-1);
                phi = obj.stm(dt);
            
                % innovation vector, J ~ N(0,Q)
                J = mvnrnd(zeros(1,obj.dim), obj.processnoise(dt), 1)';

                xs(:,i) = phi * xs(:,i-1) + J;
            end

            obj.ts = ts + obj.t0;
            obj.xs = xs;
        end

        function phi = stm(obj,tau)
            %STM Returns the state transition matrix for the time step tau.
            %   Input:
            %    - tau; time interval, in seconds
            arguments
                obj (1,1)   MarkovClock
                tau (1,1)   double {mustBeNonnegative}
            end

            phi = zeros(obj.dim, obj.dim);
            % traditional 3-state model
            phi(1:3,1:3) = [1 tau tau^2/2; 0 1 tau; 0 0 1];
            % contribution from Markov processes
            for i=1:obj.m
                phi(1,3+i) = (1 - exp(-obj.R(i)*tau))/obj.R(i);
                phi(3+i,3+i) = exp(-obj.R(i)*tau);
            end
        end

        function Q = processnoise(obj,tau)
            %PROCESSNOISE Returns the discrete-time clock process noise for
            %the time step tau.
            %   Input:
            %    - tau; time interval, in seconds
            arguments
                obj (1,1)   MarkovClock
                tau (1,1)   double {mustBeNonnegative}
            end
            
            Q = zeros(obj.dim, obj.dim);
            % traditional white + random walk model
            Q(1:2,1:2) = ...
                [obj.sigma_w^2*tau + obj.sigma_rw^2/3*tau^3 obj.sigma_rw^2/2*tau^2
                                     obj.sigma_rw^2/2*tau^2     obj.sigma_rw^2*tau];
            % contributions from Markov processes
            for i=1:obj.m
                Q(1,1) = Q(1,1) + obj.sigma_m(i)^2 * ...
                    (-3/2 + obj.R(i)*tau + 2*exp(-obj.R(i)*tau) - exp(-2*obj.R(i)*tau)/2) / obj.R(i)^3;
                Q(1,3+i) = obj.sigma_m(i)^2 * ...
                    (1/2 - exp(-obj.R(i)*tau) + exp(-2*obj.R(i)*tau)/2) / obj.R(i)^2;
                Q(3+i,1) = Q(1,3+i);
                Q(3+i,3+i) = obj.sigma_m(i)^2 * ...
                    (1 - exp(-2*obj.R(i)*tau)) / (2*obj.R(i));
            end
        end

        function assignclockdata(obj,name)
            %ASSIGNCLOCKDATA Sets a number of object properties based on
            %the name of a given oscillator, accessing data from .json.
            %   Input:
            %    - name; string name of oscillator (must match filename exactly)
            arguments
                obj     (1,1)   MarkovClock
                name    (1,:)   {mustBeText}
            end

            fname = name + ".json";

            try
                data = jsondecode(fileread(fname));
            catch
                error("assignclockdata:fileNotFound", ...
                    "File %s could not be found.", fname);
            end

            % assign data from JSON file
            obj.f = data.frequency;
            obj.a = data.aging / 86400;     % convert s/s/day -> s/s/s
            obj.t_stab  = data.stability.int;
            obj.s_stab  = data.stability.dev;
            obj.f_noise = data.phase_noise.freq;
            obj.n_noise = data.phase_noise.noise;
        end

        function [fx,C] = modelfit(obj)
            %MODELFIT Returns a second-order polynomial model for the clock
            %state over time. Starting epoch is the current t0, x0
            %   Output:
            %    - fx; @(t) function handle, input seconds past J2000 and
            %          it returns clock state
            %    - C; current bias, drift, aging used
            arguments
                obj (1,1)   MarkovClock
            end

            C = obj.x0;
            fx = @(tau) [C(1) + (tau)*C(2) + (tau).^2/2*C(3); ...
                       C(2) + (tau)*C(3); ...
                       ones(size(tau))*C(3)];
        end

        function dxdt = dynamics(obj)
            %DYNAMICS Placeholder
        end

        function A = partials(obj)
            %PARTIALS Placeholder
        end

        function P = proplyapunov(obj,P0,ts)
            %PROPLYAPUNOV Propagates discrete-time Lyapunov equations over
            %given time steps. If t0 not included, it is prepended.
            %   Inputs:
            %    - P0; starting covariance of (x(t), y(t), z(t))
            %    - ts; seconds past t0
            arguments
                obj (1,1)   MarkovClock
                P0  (3,3)   double {mustBeNonnegative}
                ts  (:,1)   double
            end

            if ts(1) ~= 0, ts = [obj.t0 ts]; end
            n = length(ts);
            P = zeros(obj.dim,obj.dim,n);
            P(1:3,1:3,1) = P0;

            for i=2:n
                tau = ts(i) - ts(i-1);
                phi = obj.stm(tau);
                P(:,:,i) = phi * P(:,:,i-1) * phi' + obj.processnoise(tau);
            end
        end
    end
end

