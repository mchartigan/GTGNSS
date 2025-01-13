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
        Q       (3,3)   double          % process noise of clock
        dim             = 3             % dimension of state
        % info stored from latest run (if ts = [], no run yet)
        ts
        xs
        % simulation seed, generated on initialization and reused on runat() for
        % consistency across runs
        seed    (1,1)   {mustBeInteger} = 0;
        % should debug output be used?
        DEBUG   (1,1)   {mustBeNonnegative,mustBeInteger} = 0
        
    end
    properties (SetAccess = private)
        % fundamental oscillator-specific info
        type    (1,:)   {mustBeText} = "Allan"          % fit type
        R2      (1,1)   double                          % R-squared value
        f       (1,1)   double {mustBePositive} = 1     % oscillator frequency, Hz
        a       (1,1)   double = 0                      % aging rate, s/s/s
        t_allan (:,1)   double {mustBePositive} = []    % Allan deviation intervals, s
        s_allan (:,1)   double {mustBePositive} = []    % Allan deviations, s/s
        f_noise (1,:)   double {mustBePositive} = []    % phase noise frequency offsets, Hz
        n_noise (1,:)   double = []                     % phase noise, dBc/Hz
    end
    properties (Constant)
        c = 299792458;                  % m/s, speed of light
    end
    
    methods
        function obj = Clock(t0,x0,name,varargin)
            %CLOCK Construct an oscillator/clock instance.
            %   Inputs:
            %    - t0; starting time in seconds past J2000
            %    - x0; starting state, [bias (s); drift (s/s); aging (1/s)]
            %    - name; string, see Clock.assignclockdata() for options
            %    - debug; optional name-value pair (default false), print debug
            %             output
            %    - fit; optional name-value pair (default "Allan"), model
            %           fit for short-term stability data
            arguments
                t0      (1,1)   double
                x0      (3,1)   double
                name    (1,:)   {mustBeText}
            end
            arguments (Repeating)
                varargin
            end
            
            obj.t0 = t0;
            obj.x0 = x0;
            obj.DEBUG = false;
            obj.type = "Allan";

            if nargin > 3
                for i=1:2:nargin-3
                    if strcmp(varargin{i}, "debug")
                        obj.DEBUG = varargin{i+1};
                    elseif strcmp(varargin{i}, "fit")
                        obj.type = varargin{i+1};
                    end
                end
            elseif nargin < 3
                error("Clock:nargin", "Clock.Clock() accepts >= 3 args.");
            end

            % parse clock type
            if ~isa(name, 'string')
                error("Clock:invalidType", "Clock type must either be string.");
            end

            if strcmpi(name, "none")    % no clock (default initialization)
                return
            else                        % get from supported clocks
                obj.assignclockdata(name);
            end

            obj.seed = randi([1 1e9]);
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
            
            rng(obj.seed)       % initialize rng for consistency

            % initialize variables
            n = length(ts);
            xs = zeros(3,n);
            % scale by c^2 to improve integration
            xs(:,1) = obj.x0 * obj.c^2;
            vs = zeros(3,3,n);

            for i=2:n
                dt = ts(i) - ts(i-1);
                stm = Clock.STM(dt);
            
                % innovation vector, J ~ N(0,Q)
                J = mvnrnd([0 0 0], obj.dtcovariance(dt), 1)';
                % covariance of x at ts(i)
                vs(:,:,i) = obj.dtcovariance(ts(i)); 
                % scale J by c^2 to improve integration
                xs(:,i) = stm * xs(:,i-1) + J * obj.c^2;
            end

            obj.ts = ts + obj.t0;
            % rescale xs back by c^2
            xs = xs ./ obj.c^2;
            obj.xs = xs;
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
            P0 = reshape(P0, obj.dim*obj.dim, 1) * obj.c^4;
            [~,Y] = ode45(@(t,p) obj.lyapunov(p, Clock.partials(t,obj.x0), obj.Q * obj.c^4), ...
                          ts, P0, odeset("RelTol", 1e-9, "AbsTol", 1e-11));

            % store covariance matrices in appropriate structure
            for i=1:length(ts)
                P(:,:,i) = reshape(Y(i,:)', obj.dim, obj.dim);
            end

            P = P ./ obj.c^4;
        end

        function s = stability(obj,dt)
            %STABILITY Returns the short-term stability of the oscillator
            %at the given measurement interval. Either Allan or Hadamard
            %variance(! not deviation), in (s/s)^2.
            %   Input:
            %    - dt; measurement interval in s
            arguments
                obj (1,1)   Clock
                dt  (1,:)   double {mustBePositive}
            end

            if strcmpi(obj.type, "Allan")
                s = obj.Q(1,1)./dt + obj.Q(2,2)*dt/3 + dt.^2/2 * obj.a^2;
            elseif strcmpi(obj.type, "Hadamard")
                s = obj.Q(1,1)./dt + obj.Q(2,2)*dt/6 + 11/120*obj.Q(3,3)*dt.^3;
            elseif strcmpi(obj.type, "constant")
                s = obj.s_allan^2;
            else
                error("stability:invalidFit", "Not sure how you got here really...");
            end
        end

        function Q = dtcovariance(obj,dt)
            %DTCOVARIANCE Returns the discrete-time covariance associated
            %w/ the Wiener processes (obj.Q) is the continuous-time
            %coviarance.
            %   Input:
            %    - dt; time step
            arguments
                obj (1,1)   Clock
                dt  (1,1)   double {mustBePositive}
            end

            v1 = obj.Q(1,1);
            v2 = obj.Q(2,2);
            v3 = obj.Q(3,3);

            Q = [v1*dt + v2/3*dt^3 + v3/20*dt^5 v2/2*dt^2 + v3/8*dt^4 v3/6*dt^3;
                 v2/2*dt^2 + v3/8*dt^4          v2*dt + v3/3*dt^3     v3/2*dt^2;
                 v3/6*dt^3                      v3/2*dt^2             v3*dt   ];
        end

        function assignclockdata(obj,name)
            %ASSIGNCLOCKDATA Sets a number of object properties based on
            %the name of a given oscillator, accessing data from JSON.
            %   Choose oscillator from list: "MicrochipCSAC", "MicrochipMAC",
            %   "SafranMAC", "SafranMiniRAFS", "RakonMiniUSO", "AccubeatUSO",
            %   "SafranRAFS", "ExcelitasRAFS"
            %
            %   Input:
            %    - name; string name of oscillator (must match file exactly)
            arguments
                obj     (1,1)   Clock
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
            obj.t_allan = data.stability.int';
            obj.s_allan = data.stability.dev';
            obj.f_noise = data.phase_noise.freq;
            obj.n_noise = data.phase_noise.noise;

            % organize variables for fitting process
            taus = obj.t_allan;
            stds = obj.s_allan;
            nrm = stds(end)^2;

            if strcmpi(obj.type, "Allan")
                % organize coefficients to solve Allan variance equation
                b = (stds.^2 - obj.a^2 / 2 * taus.^2) / nrm;
                
                % perform constrained linear least squares fit to Allan
                % variance equation
                opts = fitoptions('Method', 'LinearLeastSquares', ...
                    'Lower', [0 0], 'Robust', 'off');
                ftype = fittype({'1/x', 'x'}, 'options', opts);
                [curve, gof] = fit(taus, b, ftype);
                
                % convert and assign coefficients
                B = [coeffvalues(curve) * nrm obj.a^2/2]';
                s1 = sqrt(B(1));                % diffusion coefficient of white noise
                s2 = sqrt(3 * B(2));            % ^ of random walk frequency noise
                s3 = 0;                         % set to zero (aging rate is constant)

            elseif strcmpi(obj.type, "Hadamard")
                % organize coefficients to solve Hadamard variance equation
                b = stds.^2 / nrm;
                
                % perform constrained linear least squares fit to Hadamard
                % variance equation
                opts = fitoptions('Method', 'LinearLeastSquares', ...
                    'Lower', [0 0 0], 'Robust', 'off');
                ftype = fittype({'1/x', 'x', 'x^3'}, 'options', opts);
                [curve, gof] = fit(taus, b, ftype);

                % convert and assign coefficients
                B = coeffvalues(curve)' * nrm;
                s1 = sqrt(B(1));                % diffusion coefficient of white noise
                s2 = sqrt(B(2) * 6);            % ^ of random walk frequency noise
                s3 = sqrt(120 * B(3) / 11);     % set to zero (aging rate is constant)

            elseif strcmpi(obj.type, "constant")
                s1 = 0;
                s2 = 0;
                s3 = 0;
                gof.sse = 0;

            else
                error("stability2diffcoeff:invalidFit", ...
                    "%s not found, valid fit options found in documentation.", obj.type);
            end

            obj.R2 = gof.rsquare;           % summary statistic of fit quality
            obj.Q = diag([s1 s2 s3].^2);    % coefficients into process noise

            if obj.DEBUG
                % plot short-term stability fit to check quality
                ta = obj.t_allan(1):obj.t_allan(end);

                if strcmpi(name, "AccubeatUSO")
                    data2 = jsondecode(fileread("AltAccubeatUSO.json"));
                    t2 = data2.stability.int';
                    s2 = data2.stability.dev';
                    ta = t2(1):t2(end);
                end
                if strcmpi(obj.type, "Allan")
                    sa = sqrt(B(1) ./ ta + B(2) .* ta + obj.a^2/2 * ta .^ 2);
                else
                    sa = sqrt(B(1) ./ ta + B(2) * ta + B(3) * ta .^ 3);
                end

                entries = ["Model", "Data"];
            
                figure();
                loglog(ta, sa);
                hold on;
                scatter(obj.t_allan, obj.s_allan, 100, "rx");
                if strcmpi(name, "AccubeatUSO")
                    scatter(t2, s2, 50, "k*");
                    entries = ["Model", "Specified Data", "Sample Clock from Datasheet"];
                end
                hold off; grid on;
                xlabel("Interval (s)");
                ylabel("\sigma_Y(\tau)");
                title(sprintf("%s %s deviation fit", name, obj.type));
                legend(entries);
            end
        end

        function [err,var] = getjitter(obj,fc,Bn)
            %GETJITTER Returns the jitter noise of a clock at a specific noise
            %bandwidth, based on the phase noise statistics provided in the 
            %datasheets.
            %   Inputs:
            %    - fc; carrier frequency (to determine multiplication of
            %          clock frequency needed
            %    - Bn; carrier loop noise bandwidth
            %   Outputs:
            %    - err; sample error, in rad
            %    - var; variance of clock jitter, radians^2
            %
            %   Ref: Zucca, C. and Tavella, P.; doi.org/10.1109/TUFFC.2005.1406554
            arguments
                obj (1,1)   Clock
                fc  (1,1)   double {mustBePositive}
                Bn  (1,1)   double {mustBePositive}
            end

            Bn = Bn / 2;        % noise bandwidth presumed two-sided, so get one side
            N = fc / obj.f;
            noise = 10.^((obj.n_noise + 20*log10(N))/10);
            n_Bn = interp1(obj.f_noise, noise, Bn);
            ii = find(obj.f_noise < Bn);
            f_int = [obj.f_noise(ii) Bn];
            n_int = [noise(ii) n_Bn];
            A = trapz(f_int, n_int);
            A = 10*log10(A);
            
            var = 2*10^(A/10) * (Clock.c/(fc*2*pi))^2;
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
                obj (1,1)   Clock
            end

            C = obj.x0;
            fx = @(tau) [C(1) + (tau)*C(2) + (tau).^2/2*C(3); ...
                       C(2) + (tau)*C(3); ...
                       ones(size(tau))*C(3)];
        end

        function stability2diffcoeff(obj)
            %STABILITY2DIFFCOEFF Converts short-term stability data
            %(pre-loaded into instance properties) to diffusion
            %coefficients, stored in obj.Q.
            arguments
                obj (1,1)   Clock
            end

            taus = obj.t_allan;
            stds = obj.s_allan;
            nrm = stds(end)^2;

            if strcmpi(obj.type, "Allan")
                % organize coefficients to solve Allan variance equation
                b = (stds.^2 - obj.a^2 / 2 * taus.^2) / nrm;
                
                % perform constrained linear least squares fit to Allan
                % variance equation
                opts = fitoptions('Method', 'LinearLeastSquares', ...
                    'Lower', [0 0], 'Robust', 'off');
                ftype = fittype({'1/x', 'x'}, 'options', opts);
                [curve, gof] = fit(taus, b, ftype);
                
                % convert and assign coefficients
                B = [coeffvalues(curve) * nrm obj.a^2/2]';
                s1 = sqrt(B(1));                % diffusion coefficient of white noise
                s2 = sqrt(3 * B(2));            % ^ of random walk frequency noise
                s3 = 0;                         % set to zero (aging rate is constant)

            elseif strcmpi(obj.type, "Hadamard")
                % organize coefficients to solve Hadamard variance equation
                b = stds.^2 / nrm;
                
                % perform constrained linear least squares fit to Hadamard
                % variance equation
                opts = fitoptions('Method', 'LinearLeastSquares', ...
                    'Lower', [0 0 0], 'Robust', 'off');
                ftype = fittype({'1/x', 'x', 'x^3'}, 'options', opts);
                [curve, gof] = fit(taus, b, ftype);

                % convert and assign coefficients
                B = coeffvalues(curve)' * nrm;
                s1 = sqrt(B(1));                % diffusion coefficient of white noise
                s2 = sqrt(B(2) * 6);            % ^ of random walk frequency noise
                s3 = sqrt(120 * B(3) / 11);     % set to zero (aging rate is constant)

            elseif strcmpi(obj.type, "constant")
                s1 = 0;
                s2 = 0;
                s3 = 0;
                gof.sse = 0;

            else
                error("stability2diffcoeff:invalidFit", ...
                    "%s not found, valid fit options found in documentation.", obj.type);
            end

            obj.R2 = gof.rsquare;           % summary statistic of fit quality
            obj.Q = diag([s1 s2 s3].^2);    % coefficients into process noise

            if obj.DEBUG
                % plot short-term stability fit to check quality
                ta = obj.t_allan(1):obj.t_allan(end);
                if strcmpi(obj.type, "Allan")
                    sa = sqrt(B(1) ./ ta + B(2) .* ta + obj.a^2/2 * ta .^ 2);
                else
                    sa = sqrt(B(1) ./ ta + B(2) * ta + B(3) * ta .^ 3);
                end
            
                figure();
                loglog(ta, sa);
                hold on;
                scatter(obj.t_allan, obj.s_allan, 100, "rx");
                xlabel("Interval (s)");
                ylabel("\sigma_Y(\tau)");
                title(sprintf("%s %s deviation fit", name, obj.type));
                legend(["Model", "Data"]);
            end
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

