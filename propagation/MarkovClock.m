classdef MarkovClock < Propagator
    %MARKOVCLOCK Markov process-based clock model for propagation and
    %filtering.
    %   NOTE: Everything is normalized (multiplied) by c by default
    
    properties
        % size of state
        dim
        % no. of Markov processes
        m           (1,1)   {mustBeNonnegative,mustBeInteger} = 1
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
        % random run frequency noise
        sigma_rr    (1,1)   double {mustBeNonnegative} = 0
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
        % s/s, stability (standard deviations) from Hadamard deviations
        s_had       (:,1)   double {mustBePositive} = []
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
        function obj = MarkovClock(name,minout)
            %MARKOVCLOCK Construct an instance of MarkovClock
            %   Inputs:
            %    - t0; starting time in seconds past J2000
            %    - name; filename (no extension) of clock info .json
            %         available on path
            %    - minout;  

            % NOTE: Uncertainties (standard deviations) and states are
            % normalized by c throughout.
            
            % assign uncertainties
            obj.sigma_w  = minout.x(1) * obj.c;
            obj.sigma_rw = minout.x(2) * obj.c;
            % if output is odd, it means random run FM was fit
            fit_a = mod(length(minout.x), 2);
            obj.sigma_rr = minout.x(3) * fit_a * obj.c;
            % assign set number of Markov processes
            obj.m = (length(minout.x) - 2 - fit_a)/2;
            obj.sigma_m = zeros(obj.m, 1);
            obj.R       = ones(obj.m, 1);
            for i=1:2:2*obj.m
                obj.sigma_m(i) = minout.x(2+fit_a+i) * obj.c;
                obj.R(i)       = minout.x(3+fit_a+i);
            end

            % set state size
            obj.dim = 3 + obj.m;
            % set static seed
            obj.seed = randi([1 1e9]);
            % get info from datasheet .json
            obj.assignclockdata(name);
            obj.a = obj.a * obj.c;
        end

        function [ts,xs] = run(obj,ts,x0,n)
            %RUN Propagate the input states for tf seconds (n steps
            %between).
            %   Input:
            %    - ts; [intial time, final time], seconds past J2000
            %    - x0 (3,1) double; starting state
            %    - n; number of time steps
            %   Output:
            %    - ts; times in seconds past J2000
            %    - xs; clock states at ts
            arguments
                obj (1,1)   MarkovClock
                ts  (1,1)   double {mustBePositive}
                x0  (3,1)   double
                n   (1,1)   {mustBeInteger,mustBePositive}
            end

            % initialize variables
            ts = linspace(ts(1),ts(end),n);
            xs = obj.runat(ts,x0);
        end
        
        function [ts,xs] = runat(obj,ts,x0)
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
                ts  (1,:)   double {mustBeNonnegative}
                x0  (3,1)   double
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

            for i=2:n
                dt = ts(i) - ts(i-1);
                phi = obj.stm(dt);
            
                % innovation vector, J ~ N(0,Q)
                J = mvnrnd(zeros(1,obj.dim), obj.processnoise(dt), 1)';

                xs(:,i) = phi * xs(:,i-1) + J;
            end
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

            s1 = obj.sigma_w;
            s2 = obj.sigma_rw;
            s3 = obj.sigma_rr;
            
            Q = zeros(obj.dim, obj.dim);
            % traditional white + random walk model
            Q(1:3,1:3) = ...
                [s1^2*tau + s2^2/3*tau^3 + s3^2/20*tau^5, s2^2/2*tau^2 + s3^2/8*tau^4, s3^2/6*tau^3
                            s2^2/2*tau^2 + s3^2/8 *tau^4, s2^2  *tau   + s3^2/3*tau^3, s3^2/2*tau^2
                                           s3^2/6 *tau^3,                s3^2/2*tau^2, s3^2*tau    ];
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
            obj.s_had   = data.stability.hadamard;
            obj.f_noise = data.phase_noise.freq;
            obj.n_noise = data.phase_noise.noise;
        end

        function [err,var] = getjitter(obj,fc,Bn)
            %GETJITTER Returns the jitter noise of a clock at a specific noise
            %bandwidth, based on the phase noise statistics provided in the 
            %datasheets.
            %   Inputs:
            %    - fc; carrier frequency (to determine multiplication of
            %          clock frequency needed)
            %    - Bn; carrier loop noise bandwidth
            %   Outputs:
            %    - err; sample error, in rad
            %    - var; variance of clock jitter, rad^2
            %
            %   Ref: Zucca, C. and Tavella, P.; doi.org/10.1109/TUFFC.2005.1406554
            arguments
                obj (1,1)   MarkovClock
                fc  (1,1)   double {mustBePositive}
                Bn  (1,1)   double {mustBePositive}
            end

            % noise bandwidth presumed two-sided, so get one side
            Bn = Bn / 2;        
            N = fc / obj.f;
            noise = 10.^((obj.n_noise + 20*log10(N))/10);
            n_Bn = interp1(obj.f_noise, noise, Bn);
            ii = find(obj.f_noise < Bn);
            f_int = [0; obj.f_noise(ii); Bn];
            n_int = [noise(1); noise(ii); n_Bn];
            A = trapz(f_int, n_int);
            
            var = 2*A;
            err = mvnrnd(0, var);
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

