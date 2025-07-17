classdef Clock < Propagator
    %CLOCK Generic propagator for an oscillator/clock drift over time.
    %Provides various default oscillator options and implements all common
    %properties/methods of the abstract parent class, Propagator (also
    %shared by OrbitPropagator).
    %
    %   For additional reading, see Zucca and Tavella, 2005: 
    %   https://www.doi.org/10.1109/TUFFC.2005.1406554
    
    properties
        Q       (:,:)   double              % process noise of clock
        % simulation seed, generated on initialization and reused on runat() for
        % consistency across runs
        seed    (1,1)   {mustBeInteger} = 0;
        DEBUG   (1,1)   {mustBeNumericOrLogical} = false        % should debug output be used?
        norm    (1,1)   double = 1                              % normalization coefficient
        markov  (1,1)   {mustBeNumericOrLogical} = false        % are markov processes included?
        dim     = 3                                             % dimension of state
        m       (1,1)   {mustBeInteger,mustBeNonnegative} = 0   % number of markov processes

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
        c   = 299792458;    % m/s, speed of light
    end
    
    methods
        function obj = Clock(name,minout,options)
            %CLOCK Construct an oscillator/clock instance.
            %   Inputs:
            %    - name; filename (no extension) of clock info .json
            %         available on path
            %    - minout; output struct from ClockOpt that dictates the
            %       stability
            %    - debug; optional name-value pair (default false), print debug
            %       output
            %    - normalize; optional name-value pair (default false),
            %       multiply state by the speed of light to improve
            %       numerical accuracy
            arguments
                name                (1,:)   {mustBeText}
                minout              (1,1)   struct
                options.debug       (1,1)   {mustBeNumericOrLogical} = false
                options.normalize   (1,1)   {mustBeNumericOrLogical} = false
            end
            
            % ASSIGN FROM INPUT OPTIONS %
            obj.assignclockdata(name);
            obj.DEBUG = options.debug;
            if options.normalize, obj.norm = obj.c; end

            % ASSIGN FIT VARIANCES %
            obj.sigma_w  = minout.x(2) * obj.norm;
            obj.sigma_rw = minout.x(3) * obj.norm;
            obj.sigma_rr = minout.x(5) * obj.norm;
            % Markov process info
            obj.m = (length(minout.x) - 5)/2;
            obj.dim = 5 + obj.m;
            obj.sigma_m  = zeros(obj.m, 1);
            obj.R        = ones(obj.m, 1);
            for i=1:2:2*obj.m
                obj.sigma_m(i) = minout.x(3+i) * obj.norm;
                obj.R(i)       = minout.x(4+i);
            end

            obj.seed = randi([1 1e9]);
            obj.a = obj.a * obj.norm;
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
            %    - vs; clock covariance at ts
            arguments
                obj (1,1)   Clock
                ts  (1,1)   double {mustBePositive}
                x0  (3,1)   double
                n   (1,1)   {mustBeInteger,mustBePositive}
            end

            % initialize variables
            ts = linspace(ts(1),ts(end),n);
            xs = obj.runat(ts,x0);
        end

        function xs = runat(obj,ts,x0)
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
                x0  (3,1)   double
            end
            
            rng(obj.seed)       % initialize rng for consistency

            % initialize variables
            n = length(ts);
            xs = zeros(obj.dim,n);
            xs(1:3,1) = x0;
            % set starting state of Markov processes as RV with mean 0 and
            % variance U = sigma_m^2/(2*R)
            for i=1:obj.m
                xs(3+i,1) = mvnrnd(0, obj.sigma_m(i)^2/(2*obj.R(i)));
            end

            for i=2:n
                dt = ts(i) - ts(i-1);
                stm = obj.STM(dt);
            
                % innovation vector, J ~ N(0,Q)
                J = mvnrnd(zeros(1,obj.dim), obj.processnoise(dt), 1)';
                xs(:,i) = stm * xs(:,i-1) + J;
            end
        end

        function assignclockdata(obj,name)
            %ASSIGNCLOCKDATA Sets a number of object properties based on
            %the name of a given oscillator, accessing data from .json.
            %   Input:
            %    - name; string name of oscillator (must match filename exactly)
            arguments
                obj     (1,1)   Clock
                name    (1,:)   {mustBeText}
            end

            % handle empty case
            if strcmpi(name, "none")
                obj.f = 1;
                obj.a = 0;
                obj.t_stab = 1;
                obj.s_stab = 0;
                obj.s_had  = 0;
                obj.f_noise = 1;
                obj.n_noise = -Inf;
                return;
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
            P(:,:,1) = P0;

            % store covariance matrices in appropriate structure
            for i=2:length(ts)
                Phi = obj.STM(ts(i)-ts(i-1));
                Qi = obj.processnoise(ts(i)-ts(i-1));
                P(:,:,i) = Phi*P(:,:,i-1)*Phi' + Qi;
            end
        end

        function s = stability(obj,dt)
            %STABILITY Returns the short-term stability of the oscillator
            %at the given measurement interval. Hadamard variance(! not
            %deviation), in (s/s)^2 or (m/s)^2, depending on normalization.
            %   Input:
            %    - dt; measurement interval in s
            arguments
                obj (1,1)   Clock
                dt  (1,:)   double {mustBePositive}
            end

            s1 = obj.sigma_w;
            s2 = obj.sigma_rw;
            s3 = obj.sigma_rr;
            s = s1^2./dt + s2^2*dt/6 + 11/120*s3^2*dt.^3;
        end

        function Q = processnoise(obj,dt)
            %PROCESSNOISE Returns the discrete-time covariance associated
            %w/ the Wiener processes.
            %   Input:
            %    - tau; time step

            s1 = obj.sigma_w;
            s2 = obj.sigma_rw;
            s3 = obj.sigma_rr;
            
            Q = zeros(obj.dim, obj.dim);
            % traditional white + random walk model
            Q(1:3,1:3) = ...
                [s1^2*dt + s2^2/3*dt^3 + s3^2/20*dt^5, s2^2/2*dt^2 + s3^2/8*dt^4, s3^2/6*dt^3
                           s2^2/2*dt^2 + s3^2/8 *dt^4, s2^2  *dt   + s3^2/3*dt^3, s3^2/2*dt^2
                                         s3^2/6 *dt^3,               s3^2/2*dt^2, s3^2*dt    ];
            % contributions from Markov processes
            for i=1:obj.m
                Q(1,1) = Q(1,1) + obj.sigma_m(i)^2 * ...
                    (-3/2 + obj.R(i)*dt + 2*exp(-obj.R(i)*dt) - exp(-2*obj.R(i)*dt)/2) / obj.R(i)^3;
                Q(1,3+i) = obj.sigma_m(i)^2 * ...
                    (1/2 - exp(-obj.R(i)*dt) + exp(-2*obj.R(i)*dt)/2) / obj.R(i)^2;
                Q(3+i,1) = Q(1,3+i);
                Q(3+i,3+i) = obj.sigma_m(i)^2 * ...
                    (1 - exp(-2*obj.R(i)*dt)) / (2*obj.R(i));
            end
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
                obj (1,1)   Clock
                fc  (1,1)   double {mustBePositive}
                Bn  (1,1)   double {mustBePositive}
            end

            % noise bandwidth presumed two-sided, so get one side
            Bn = Bn / 2;        
            N = fc / obj.f;
            noise = 10.^((obj.n_noise + 20*log10(N))/10);
            n_Bn = interp1(obj.f_noise, noise, Bn);
            ii = find(obj.f_noise < Bn);
            f_int = [obj.f_noise(ii) Bn];
            n_int = [noise(ii) n_Bn];
            A = trapz(f_int, n_int);
            
            var = 2*A;
            err = mvnrnd(0, var);
        end

        function stm = STM(obj,dt)
            %STM Returns the DT state transition matrix based on the dynamics
            %defined in Zucca and Tavella.
            %   Input:
            %    - dt; time step

            stm = [1 dt dt^2/2; 0 1 dt; 0 0 1];
            % contribution from Markov processes
            for i=1:obj.m
                stm(1,3+i) = (1 - exp(-obj.R(i)*dt))/obj.R(i);
                stm(3+i,3+i) = exp(-obj.R(i)*dt);
            end
        end

        function plot(obj,traj)
            %PLOT Plots the phase, freq. offset, and freq. drift for the
            %provided trajectory.
            %   Input:
            %    - traj; Trajectory object for Clock output
            arguments
                obj     (1,1)   Clock
                traj    (1,1)   Trajectory
            end

            ts = traj.ts;
            xs = traj.xs;

            if obj.norm == 1
                units = "ns";
                xs = xs * 1e9;
            else
                units = "m";
            end

            dt = ts - ts(1);

            if dt(end) > 86400
                time = "days";
                tplot = dt / 86400;
            elseif dt(end) > 3600
                time = "hrs";
                tplot = dt / 3600;
            else
                time = "s";
            end

            
            figure();
            plotformat("APA", 0.9);
            tiledlayout(3,1);
            
            nexttile;
            plot(tplot, xs(1,:));
            grid on;
            ylabel(sprintf("Phase offset (%s)", units));
            
            nexttile;
            plot(tplot, xs(2,:));
            grid on;
            ylabel(sprintf("Freq. offset (%s/s)", units));
            
            nexttile;
            plot(tplot, xs(3,:));
            grid on;
            ylabel(sprintf("Freq. drift (%s/s^2)", units));
            xlabel(sprintf("Time (%s)", time));

            sgtitle("Clock trajectory");
        end
    end

    methods (Static)
        function dxdt = dynamics(~,x)
            %DYNAMICS Invokes the clock dynamics based on the Zucca and
            %Tavella paper.
            %   Input:
            %    - t; simulation time, not used
            %    - x; current clock state

            dxdt = [0 1 0; 0 0 1; 0 0 0] * x;
        end

        function A = partials(~,~)
            %PARTIALS Invokes the partials of dynamics.
            %   Input:
            %    - t; simulation time, not used
            %    - x; current clock state, also not used lol

            A = [0 1 0; 0 0 1; 0 0 0];
        end

        function [fx,C] = modelfit(traj,t0)
            %MODELFIT Returns a second-order polynomial model for the clock
            %state over time. Starting epoch is the current t0, x0
            %   Input:
            %    - traj; Trajectory object of clock states
            %    - t0; starting epoch (s past J2000)
            %   Output:
            %    - fx; @(t) function handle, input seconds past J2000 and
            %          it returns clock state
            %    - C; current bias, drift, aging used
            arguments
                traj    (1,1)   Trajectory
                t0      (1,1)   double
            end

            C = traj.get(t0);
            fx = @(tau) [C(1) + (tau)*C(2) + (tau).^2/2*C(3); ...
                       C(2) + (tau)*C(3); ...
                       ones(size(tau))*C(3)];
        end
    end
end

