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
        % override noise param (for filtering)
        random  = true

        % FIT UNCERTAINTY PARAMETERS %
        % white phase noise
        sigma_wp    (1,1)   double {mustBeNonnegative} = 0
        % white frequency noise
        sigma_wf    (1,1)   double {mustBeNonnegative} = 0
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
        s_stab      (:,1)   double {mustBeNonnegative} = []
        % s/s, stability (standard deviations) from Hadamard deviations
        s_had       (:,1)   double {mustBeNonnegative} = []
        % Hz, phase noise frequency offsets
        f_noise     (1,:)   double {mustBePositive} = []
        % dBc/Hz, phase noise
        n_noise     (1,:)   double = []
    end
    properties (Constant)
        c   = 299792458     % m/s, speed of light
    end
    
    methods
        function obj = Clock(name,x,options)
            %CLOCK Construct an oscillator/clock instance.
            %   Inputs:
            %    - name; filename (no extension) of clock info .json
            %         available on path
            %    - minout; output from ClockOpt that dictates the
            %       stability, in form
            %          [sigma_WP sigma_WF sigma_RW sigma_RR sigma_m1 R1 ... sigma_mn Rn]'
            %       all entries after sigma_RR are optional. units in s 
            %    - debug; optional name-value pair (default false), print debug
            %       output
            %    - normalize; optional name-value pair (default false),
            %       multiply state by the speed of light to improve
            %       numerical accuracy
            arguments
                name                (1,:)   {mustBeText}
                x                   (:,1)   double {mustBeNonnegative}
                options.debug       (1,1)   {mustBeNumericOrLogical} = false
                options.normalize   (1,1)   {mustBePositive} = 1
            end
            
            % ASSIGN FROM INPUT OPTIONS %
            obj.assignclockdata(name);
            obj.DEBUG = options.debug;
            obj.norm = options.normalize;

            % ASSIGN FIT VARIANCES %
            obj.sigma_wp = x(1) * obj.norm;
            obj.sigma_wf = x(2) * obj.norm;
            obj.sigma_rw = x(3) * obj.norm;
            obj.sigma_rr = x(4) * obj.norm;
            % Markov process info
            obj.m = (length(x) - 4)/2;
            obj.dim = 3 + obj.m;
            obj.sigma_m  = zeros(obj.m, 1);
            obj.R        = ones(obj.m, 1);
            for i=1:2:2*obj.m
                obj.sigma_m(i) = x(4+i) * obj.norm;
                obj.R(i)       = x(5+i);
            end

            obj.seed = randi([1 1e9]);
            obj.a = obj.a * obj.norm;
        end

        function [ts,xs] = run(obj,ts,x0,n,noise)
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
                obj     (1,1)   Clock
                ts      (1,2)   double {mustBePositive}
                x0      (:,1)   double
                n       (1,1)   {mustBeInteger,mustBePositive}
                noise   (1,1)   = true
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
                obj     (1,1)   Clock
                ts      (1,:)   double {mustBeNonnegative}
                x0      (:,1)   double
                noise   (1,1)   = true
            end
            
            % rng(obj.seed)       % initialize rng for consistency

            % initialize variables
            n = length(ts);
            xs = zeros(obj.dim,n);
            xs(:,1) = x0;
            % set starting state of Markov processes as RV with mean 0 and
            % variance U = sigma_m^2/(2*R)
            for i=1:obj.m
                xs(3+i,1) = mvnrnd(0, obj.sigma_m(i)^2/(2*obj.R(i)));
            end

            for i=2:n
                dt = ts(i) - ts(i-1);
                stm = obj.STM(dt);
            
                xs(:,i) = stm * xs(:,i-1);
                if noise && obj.random
                    % innovation vector, J ~ N(0,Q)
                    J = mvnrnd(zeros(1,obj.dim), obj.noise(dt), 1)';
                    xs(:,i) = xs(:,i) + J;
                end
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
            if isfield(data.stability, "hadamard")
                obj.s_had   = data.stability.hadamard;
            end
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
                P0  (:,:) double
            end

            n = length(ts);
            P = zeros(obj.dim, obj.dim, n);
            
            % reshape starting P to correct format
            P(:,:,1) = P0;

            % store covariance matrices in appropriate structure
            for i=2:length(ts)
                Phi = obj.STM(ts(i)-ts(i-1));
                Qi = obj.noise(ts(i)-ts(i-1));
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

            s1 = obj.sigma_wf;
            s2 = obj.sigma_rw;
            s3 = obj.sigma_rr;

            % contribution of white frequency modulation
            part_WFM  = s1^2 ./ dt;
            % ... random walk frequency modulation
            part_RWFM = s2^2/6 * dt;
            % ... random run frequency modulation
            part_RRFM = 11/120*s3^2 * dt.^3;
            % ... stationary Markov (Wiener) processes
            part_M = zeros(size(part_WFM));
            for i=1:obj.m
                sm = obj.sigma_m(i);
                Ri = obj.R(i);
                part_M = part_M + ...
                    sm^2*(Ri*dt - 5/3 + 5/2*exp(-Ri*dt) - exp(-2*Ri*dt) + exp(-3*Ri*dt)/6) ./ ...
                    (Ri^3*dt.^2);
            end

            s = part_WFM + part_RWFM + part_RRFM + part_M;
        end

        function y = adev(obj,dt)
            %ADEV Returns the Allan deviation at time intervals ts. Clock
            %behavior is governed by x.
            %   Input:
            %    - x (obj.dim,1) double {mustBeNonnegative}; variances and
            %       time constants for the 3-state clock model with obj.m
            %       Markov processes.
            %    - ts (1,:) double {mustBeNonnegative}; evaluation times,
            %       in seconds

            s1 = obj.sigma_wf;
            s2 = obj.sigma_rw;

            y = s1^2./dt + s2^2*dt/3 + obj.a^2*dt.^2/2;

            % for j=6:2:obj.dim
            %     sm = x(j);
            %     Rm = x(j+1);
            %     y = y + sm^2 * (-3/2 + Rm*ts + 2*exp(-Rm*ts) - ...
            %                     exp(-2*Rm*ts)/2) ./ (Rm^3 * ts.^2);
            % end

            y = sqrt(y);
        end

        function Q = noise(obj,dt,~)
            %PNC Returns the discrete-time process noise covariance.
            %   Input:
            %    - tau; time step

            s1 = obj.sigma_wf;
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

            % catch if clock is none
            if isinf(obj.n_noise)
                err = 0;
                var = 0;
                return;
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

        function [pos,vel] = getsisecontrib(obj,tf,P0)
            %GETSISECONTRIB Returns the signal-in-space error contribution
            %of the clock.
            %   Input:
            %    - tf; end time projection
            %    - P0; starting clock covariance
            %   Output:
            %    - pos; position variance Trajectory (m)^2
            %    - vel; velocity variance Trajectory (mm/s)^2

            ts = 0:1:tf;
            T = 10;
            fc = 2492.028e6;
            Bn = 20;
            P = obj.proplyapunov(ts,P0);
            % white phase measurement noise
            [~,v_WP] = obj.getjitter(fc, Bn);
            % convert to m
            v_WP = v_WP * (2*pi*fc)^(-2) * obj.c^2;

            pos = reshape(P(1,1,:), 1, []) + v_WP;
            vel = reshape(P(2,2,:), 1, []);
            S = obj.noise(T);
            vel = vel(1:end-5) + 2*v_WP/T^2 + S(1,1)/T^2;

            pos = Trajectory(ts, pos);
            vel = Trajectory(ts(6:end), 1e6*vel);
        end

        function stm = STM(obj,dt)
            %STM Returns the DT state transition matrix based on the dynamics
            %defined in Zucca and Tavella.
            %   Input:
            %    - dt; time step

            stm = zeros(obj.dim,obj.dim);
            stm(1:3,1:3) = [1 dt dt^2/2; 0 1 dt; 0 0 1];
            % contribution from Markov processes
            for i=1:obj.m
                stm(1,3+i) = (1 - exp(-obj.R(i)*dt))/obj.R(i);
                stm(3+i,3+i) = exp(-obj.R(i)*dt);
            end
        end

        function [ax,tplot] = plot(obj,traj,options)
            %PLOT Plots the phase, freq. offset, and freq. drift for the
            %provided trajectory.
            %   Input:
            %    - traj; Trajectory object for Clock output
            arguments
                obj     (1,1)   Clock
                traj    (1,1)   Trajectory
                options.scale  (1,:) double = []
                options.labels (1,:) {mustBeText} = strings(0)
                options.axes   (1,:) matlab.graphics.axis.Axes = []
            end

            ts = traj.ts;
            xs = traj.xs;

            % error out if wrong units are provided
            if ~all(size(options.scale) == size(options.labels)) || ~ismember(length(options.scale), [0 3])
                error("plot:invalidArg", ...
                    "If custom units are provided, scale and labels must be the same size.");
            end

            if isempty(options.scale)
                if obj.norm == 1
                    units = ["ns","ns/s","ns/s^2"];
                    xs = xs * 1e9;
                elseif obj.norm == 1e9
                    units = ["ns","ns/s","ns/s^2"];
                elseif obj.norm == obj.c
                    units = ["m","m/s","m/s^2"];
                else
                    error("Clock:plot", "Unsupported normalization scheme.");
                end
            else
                xs = xs / obj.norm;

                units = options.labels;
                xs(1,:) = xs(1,:) * options.scale(1);
                xs(2,:) = xs(2,:) * options.scale(2);
                xs(3,:) = xs(3,:) * options.scale(3);

                for i=4:size(xs,1)
                    xs(i,:) = xs(i,:) * options.scale(2);
                end
            end

            dt = ts - ts(1);

            if dt(end) > 86400
                time = "days";
                tplot = dt / 86400;
            elseif dt(end) > 3600
                time = "hrs";
                tplot = dt / 3600;
            elseif dt(end) > 120
                time = "min";
                tplot = dt / 60;
            else
                time = "s";
                tplot = dt;
            end


            if length(options.axes) < 3
                figure();
                plotformat("APA", 0.3*size(xs,1));
                tiledlayout(size(xs,1),1);
                
                a1 = nexttile;
                plot(tplot, xs(1,:), Color="k", LineStyle=":");
                grid on;
                ylabel(sprintf("Phase offset (%s)", units(1)));
                
                a2 = nexttile;
                plot(tplot, xs(2,:), Color="k", LineStyle=":");
                grid on;
                ylabel(sprintf("Freq. offset (%s)", units(2)));

                ax = [a1 a2];

                for i=4:size(xs,1)
                    ai = nexttile;
                    ax = [ax ai];
                    plot(tplot, xs(i,:), Color="k", LineStyle=":");
                    grid on;
                    ylabel(sprintf("Freq. drift (%s)", units(2)));
                end
                
                a3 = nexttile;
                plot(tplot, xs(3,:), Color="k", LineStyle=":");
                grid on;
                ylabel(sprintf("Freq. drift (%s)", units(3)));
                xlabel(sprintf("Time (%s)", time));

                ax = [ax a3];
                
                % sgtitle("Clock trajectory");
            else
                hold(options.axes(1), "on");
                plot(options.axes(1), tplot, xs(1,:), Color="k", ...
                    LineStyle=":", HandleVisibility="off");
                hold(options.axes(2), "on");
                plot(options.axes(2), tplot, xs(2,:), Color="k", ...
                    LineStyle=":", HandleVisibility="off");
                hold(options.axes(3), "on");
                plot(options.axes(3), tplot, xs(3,:), Color="k", ...
                    LineStyle=":", HandleVisibility="off");

                ax = options.axes;
            end
        end

        function dxdt = dynamics(obj,~,x)
            %DYNAMICS Invokes the clock dynamics based on the Zucca and
            %Tavella paper.
            %   Input:
            %    - t; simulation time, not used
            %    - x; current clock state

            dxdt = obj.partials() * x;
        end

        function A = partials(obj,~,~)
            %PARTIALS Invokes the partials of dynamics.
            %   Input:
            %    - t; simulation time, not used
            %    - x; current clock state, used for dimension
            A = zeros(3+obj.m,3+obj.m);

            A(1:3,1:3) = [0 1 0; 0 0 1; 0 0 0];

            if obj.m > 0
                A(1,4:end) = 1;
                A(4:end,4:end) = diag(-obj.R);
            end
        end
    end

    methods (Static)
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

