classdef NavSatellite < handle
    %NavSatellite Class for describing the properties and trajectory of a
    %satellite that is part of a radionavigation satellite system (whether
    %that's earth-based GNSS or lunar).
    
    properties
        % propagator instance (default one so MATLAB doesn't throw a fit)
        prop    (1,1)   OrbitPropagator = OrbitPropagator(0,1)
        % clock instance (default one again ^)
        clock   (1,1)   Clock = Clock(0,zeros(3,1),"none")
        % nav filter
        filter  (1,1)
        % antenna object
        ant     (1,1)   TransmitAntenna
        % reference state info (avoids recalling runto() on prop and clock
        % if data has already been requested before
        tr      (1,:)   double
        xr      (9,:)   double
        frame_r (1,:)   char = ''
        % navigation state info (avoids rerunning filter if data has already 
        % been requested before)
        tn      (1,:)   double
        x0      (:,:)   double      % doubles as default navigation info
        P0      (9,9)   double      % doubles as default navigation info
        xn      (9,:)   double
        Pn      (9,9,:) double
    end

    properties (Constant)
        % speed of light (m/s)
        c       (1,1)   double = 299792458
    end
    
    methods
        function obj = NavSatellite(prop,clock,filter,meas,ant,debug)
            %NAVSATELLITE Construct a NavSatellite instance.
            %   Inputs:
            %    - prop; propagator instance for generating future trajectories
            %            (OrbitPropagator, LunarPropagator, EarthPropagator)
            %    - clock; instance of s/c oscillator/clock (Clock)
            %    - filter; navigation filter type for finding s/c uncertainty,
            %              options are: "EKF"
            %    - meas; struct containing measurement info for the chosen
            %            filter -- names should correspond to filter constructor
            %            arguments (may also include varargin)
            %    - ant; TransmitAntenna object containing sat info
            %    - debug; should debug warnings be printed? true/false
            arguments
                prop    (1,1)   OrbitPropagator
                clock   (1,1)   Clock
                filter  (1,:)   {mustBeText}
                meas    (1,1)   struct
                ant     (1,1)   TransmitAntenna
                debug   (1,1)   = false
            end

            % assign instances
            obj.prop = prop;
            obj.clock = clock;
            obj.ant = ant;

            if strcmp(filter, "EKF")
                % CT process noise
                Q = [obj.prop.var zeros(6,3); zeros(3,6) obj.clock.var];

                % build filter, w/ or w/o optional args
                if isfield(meas, "varargin")
                    obj.filter = EKF("hybrid", @obj.dynamics, @obj.partials, ...
                        Q, @meas.h, meas.y, @meas.dhdx, meas.R, meas.t_meas, ...
                        meas.varargin);
                else
                    obj.filter = EKF("hybrid", @obj.dynamics, @obj.partials, ...
                        Q, @meas.h, meas.y, @meas.dhdx, meas.R, meas.t_meas);
                end
            else
                error("NavSatellite:invalidFilter", ...
                    'Filter must be "EKF". See documentation.');
            end

            if ~debug, warning('off', 'NavSatellite:debug'); end
        end

        function dxdt = dynamics(obj,t,x)
            %DYNAMICS Returns the full-state dynamics for a NavSatellite
            %instance. Incorporates both orbital and oscillator dynamics.
            %   Input:
            %    - t; seconds past J2000
            %    - x; s/c state
            arguments
                obj (1,1)   NavSatellite
                t   (1,1)   double
                x   (9,1)   double
            end

            dxdt = [obj.prop.dynamics(t,x(1:6)); obj.clock.dynamics(t,x(7:9))];
        end

        function A = partials(obj,t,x)
            %PARTIALS Invokes the partials of dynamics, with relevant settings, 
            %for this propagator instance. Incorporates both orbital and 
            %oscillator dynamics.
            %   Input:
            %    - t; simulation time in seconds past J2000
            %    - x; satellite state
        
            A1 = obj.prop.partials(t,x(1:6));
            A2 = obj.clock.partials(t,x(7:9));
            A = [A1 zeros(6,3); zeros(3,6) A2];
        end
        
        function xr = getrefstates(obj,ts,frame)
            %GETREFSTATES Returns reference trajectory information for the
            %satellite at the provided times and frame. Also referred to as
            %the true state.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - frame; reference frame to return data in
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                frame   (1,:)   char
            end

            if ts(1) < obj.prop.t0
                error("getrefstates:invalidInput", ...
                    "Simulation times must be after starting epoch of propagator.");
            end

            % avoid recalling if data has already been generated, as clock
            % will create a random new trajectory
            if all(size(ts) == size(obj.tr)) && ~any(ts - obj.tr) && strcmp(frame, obj.frame_r)
                % data has been generated previously
                warning("NavSatellite:debug", "Returning previously generated data.\n");
                xr = obj.xr;
                return;
            end

            % new data request
            ts = ts - obj.prop.t0;
            [~,xs] = obj.prop.runat(ts,frame);
            [~,xc] = obj.clock.runat(ts);
            xr = [xs; xc];

            % store for future calls
            obj.tr = ts + obj.prop.t0;
            obj.frame_r = frame;
            obj.xr = xr;
        end

        function [xn,Pn] = getnavstates(obj,ts,x0,P0)
            %GETNAVSTATES Returns estimated state and covariance from the 
            %navigation filter at the provided times in the J2000 frame.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - x0; state estimate at ts(1) in J2000 frame -- if empty,
            %          noise P0 will be applied to reference trajectory to get
            %          nav trajectory
            %    - P0; uncertainty of x0 at ts(1) in J2000 frame
            %    - frame; reference frame to return data in
            %   Output:
            %    - xn; state estimates at ts
            %    - Pn; state covariances at ts
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                x0      (:,:)   double
                P0      (9,9)   double
            end

            % avoid recalling if data has already been generated
            if all(size(ts) == size(obj.tn)) && ~any(ts - obj.tn) && ...
                    all(x0 == obj.x0) && all(P0 == obj.P0, 'all')
                % data has been generated previously
                warning("NavSatellite:debug", "Returning previously generated data.\n");
                xn = obj.xn;
                Pn = obj.Pn;
                return;
            end

            if isempty(x0)
                xn = obj.getrefstates(ts, 'J2000');
                % add noise to reference trajectory for nav states
                xn = mvnrnd(xn', P0)';
                Pn = repmat(P0, 1, 1, length(ts));
            elseif all(size(x0) == [9 1])
                error("getnavstates:noImplementedError", ...
                    "Navigation filtering feature is not yet implemented.");

                % add ts to obj.filter.t with union() (maybe strip everything
                % before ts(1)?)
                % run filter, get data at time steps, and return it
            else
                error("getnavstates:invalidInput", ...
                    "Input x0 must either be size [0 0] or [9 1].");
            end

            % store data in case called again
            obj.tn = ts;
            obj.x0 = x0;
            obj.P0 = P0;
            obj.xn = xn;
            obj.Pn = Pn;
        end

        function [xm,xp,Pm,Pp] = getmodelstates(obj,ts,x0,P0,cadence)
            %GETMODELSTATES Returns model-derived spacecraft states (PVT).
            %This is the information a user would be able to get from data
            %in the nav message. It encompasses ephemeris and clock model
            %errors when differenced from getrefstates().
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - x0; state estimate at ts(1) in J2000 frame -- if empty,
            %          noise P0 will be applied to reference trajectory to get
            %          nav trajectory
            %    - P0; uncertainty of x0 at ts(1) in J2000 frame
            %    - cadence; # of seconds after which models are updated, can be
            %               separate for ephemeris and clock -- 0 means don't
            %               update
            %   Output:
            %    - xm; model states, pos/vel in J2000. If Tm provided,
            %          becomes 2-page matrix where second page is xm(ts-Tm)
            %    - xp; propagated states that models are based on. If Tm given,
            %          becomes 2-page matrix where second page is xp(ts-Tm)
            %    - Pm; covariance of model fit from nav soln
            %    - Pp; covariance of propagation
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                x0      (:,:)   double
                P0      (9,9)   double
                cadence (1,2)   double {mustBeNonnegative} = zeros(1,2);
            end

            n = length(ts);
            [xnav, Pnav] = obj.getnavstates(ts, x0, P0);

            % find indices for update intervals of models
            propupdate = 1;
            clkupdate = 1;
            if any(cadence)
                t_prop = ts(1);
                t_clk = ts(1);
                for i=2:n-1
                    % as a quirk, require model updates to happen at odd
                    % indices to avoid poor data transitions for delta
                    % pseudoranges. i.e. the case where the first
                    % pseudorange is pre-update and the second is post
                    if ts(i) >= t_prop + cadence(1) && mod(i,2)
                        propupdate = [propupdate i];
                        t_prop = ts(i);
                    end
                    if ts(i) >= t_clk + cadence(2) && mod(i,2)
                        clkupdate = [clkupdate i];
                        t_clk = ts(i);
                    end
                end
            end
            propupdate = [propupdate n];
            clkupdate = [clkupdate n];

            % store propagator info for reassignment
            t0p = obj.prop.t0;
            x0p = obj.prop.x0;
            t0c = obj.clock.t0;
            x0c = obj.clock.x0;

            % generate models and get data
            xm = zeros(9, n);
            if nargout > 1, xp = zeros(9, n); end
            if nargout > 2, Pm = zeros(9, 9, n); end
            if nargout > 3, Pp = zeros(9, 9, n); end

            for i=1:length(propupdate)-1        % for orbit
                obj.prop.t0 = ts(propupdate(i));            % new start time
                obj.prop.x0 = xnav(1:6, propupdate(i));     % new start state
                tf = ts(propupdate(i+1)) - obj.prop.t0;     % end of validity interval
                fs = obj.prop.modelfit("Kepler", tf, 12);   % ephemeris model
                int = propupdate(i):propupdate(i+1);        % validity interval indices
                xm(1:6,int) = fs(ts(int));                  % assign model states

                if nargout > 1
                    % assign propagated states model is based on
                    [~,temp] = obj.prop.runat(ts(int) - ts(int(1)), 'J2000');
                    xp(1:6,int) = temp;
                end
                if nargout > 2
                    % Compute variance of states between model and propagation
                    diff = xm(1:6,int) - xp(1:6,int);
                    Pm(1:6,1:6,int) = repmat(diag(var(diff, 0, 2)), 1, 1, length(int));
                end
                if nargout > 3
                    % propagate nav uncertainty while model is valid
                    Pp(1:6,1:6,int) = ...
                        obj.prop.proplyapunov(ts(int), Pnav(1:6,1:6,int(1)));
                end
            end

            for i=1:length(clkupdate)-1         % for clock
                obj.clock.t0 = ts(clkupdate(i));            % new start time
                obj.clock.x0 = xnav(7:9, clkupdate(i));     % new start state
                fc = obj.clock.modelfit();                  % clock model
                int = clkupdate(i):clkupdate(i+1);          % validity interval indices
                xm(7:9,int) = fc(ts(int));                  % assign model states

                if nargout > 1
                    % assign propagated states model is based on (same for clock)
                    xp(7:9,int) = fc(ts(int));
                end
                if nargout > 3
                    % propagate nav uncertainty while model is valid
                    Pp(7:9,7:9,int) = ...
                        obj.clock.proplyapunov(ts(int), Pnav(7:9,7:9,int(1)));
                end
            end

            % reassign propagator info after it's been changed
            obj.prop.t0 = t0p;
            obj.prop.x0 = x0p;
            obj.clock.t0 = t0c;
            obj.clock.x0 = x0c;
        end

        % TODO: implement getsmartmodelstates() where updates are adaptive

        function [meas,var,true] = getmeasurements(obj,ts,user,opts)
            %GETMEASUREMENTS Computes pseudorange and pseudorange-rate
            %(Doppler) measurements at times ts between NavSatellite and
            %the user.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - user; User object instance
            %    - frame; reference frame user data is provided in
            %    - opts; s/c settings struct all fields must be set:
            %             - x0 (see getnavstates input)
            %             - P0 (see getnavstates input)
            %             - cadence (see getmodelstates input)
            %   Output:
            %    - meas; meas.psr is measured pseudorange measurements, 
            %            meas.psrr is doppler measurements
            %    - var; contains the overall variance error budget, incl 
            %           breakdowns by category, i.e. what the user could
            %           calculate
            %    - true; contains the true range and range-rate data, as
            %            well as errors broken down by source
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                user    (1,1)   User
                opts    (1,1)   struct
            end

            n = length(ts);
            % get true measurements
            [tt, r, dr] = obj.timeofflight(ts, user);

            % Whether it's a PLL or FLL, signal frequency is not observed
            % directly; rather it's measured as a difference of two phases.
            % Thus, clock states need to be obtained at measurement
            % spacings.
            % is the user capable of measuring range-rate? (has carrier tracking)
            rate = ~strcmpi(user.rec.carrier, "none");
            if rate
                % determine rate measurement spacing (carrier predetection
                % integration time, unless Tm is specified / larger)
                if user.rec.T_c > user.rec.Tm
                    Tm = user.rec.T_c;
                else
                    Tm = user.rec.Tm;
                end

                % catch if Tm is too large
                if Tm > min(ts(2:end)-ts(1:end-1))
                    error("getmeasurements:invalidSpacing", ...
                        "User Receiver.Tm must not be greater than minimum spacing between ts.");
                end

                % interleave times, i.e.
                % [ts(1)-Tm, ts(1), ts(2)-Tm, ts(2), ts(3)-Tm, ts(3), ... ]
                tt = [tt - Tm; tt];
                tt = tt(:)';
                ts = [ts - Tm; ts];
                ts = ts(:)';
            end

            % get user state info for measurements
            xref = obj.getrefstates(tt, 'J2000');
            xuser = user.getstate(ts, 'J2000');
            if nargout > 2
                [xmdl,xprop,Pmdl,Pprop] = obj.getmodelstates(tt, ...
                        opts.x0, opts.P0, opts.cadence);
            else
                [xmdl,xprop] = obj.getmodelstates(tt, opts.x0, ...
                        opts.P0, opts.cadence);
            end

            if rate
                % un-interleave times
                n2 = length(tt);
                tt = tt(2:2:n2);
                ts = ts(2:2:n2);

                % separate states
                % store concatenated data and create 2-page 3D matrices, then
                % un-interleave data, storing in 2 separate pages of matrix
                % first tt(i)-Tm comes, then tt(i), so store first value on 
                % second page and second on first page
                xref = cat(3, xref(:,2:2:n2), xref(:,1:2:n2));
                xuser = cat(3, xuser(:,2:2:n2), xuser(:,1:2:n2));
                xmdl = cat(3, xmdl(:,2:2:n2), xmdl(:,1:2:n2));
                xprop = cat(3, xprop(:,2:2:n2), xprop(:,1:2:n2));
                % remove tt-Tm data from covariances
                if nargout > 2
                    Pmdl = Pmdl(:,:,2:2:n2);
                    Pprop = Pprop(:,:,2:2:n2);
                end
            end

            psr_units = 1e3;    % convert distances from km to m
            psrr_units = 1e6;   % convert speeds from km/s to mm/s

            % create line-of-sight direction
            los = xref(1:3,:,:) - xuser(1:3,:,:);
            los = los ./ sqrt(sum(los.^2, 1));

            % computer state errors from propagation and model
            err_prop = xprop - xref;
            err_mdl = xmdl - xprop;

            % project error onto line-of-sight direction
            err.psr.prop = zeros(1,n);      % pseudorange error from nav uncertainty propagation
            err.psr.mdl = zeros(1,n);       % pseudorange error from model fitting
            err.psrr.prop = zeros(1,n);     % Doppler error from nav uncertainty propagation
            err.psrr.mdl = zeros(1,n);      % Doppler error from model fitting
            for i=1:n
                err.psr.prop(i) = err_prop(1:3,i,1)' * los(:,i,1) * psr_units;
                err.psr.mdl(i) = err_mdl(1:3,i,1)' * los(:,i,1) * psr_units;

                if rate         % range-rate measurements capable
                    err.psrr.prop(i) = (err_prop(1:3,i,1)'*los(:,i,1) - ...
                        err_prop(1:3,i,2)'*los(:,i,2)) / Tm * psrr_units;
                    err.psrr.mdl(i) = (err_mdl(1:3,i,1)'*los(:,i,1) - ...
                        err_mdl(1:3,i,2)'*los(:,i,2)) / Tm * psrr_units;
                end
            end
            % pseudorange error from clock model
            err.psr.clk = (err_mdl(7,:,1) + err_prop(7,:,1)) * obj.c;       % s to m
            % velocity error from clock model
            if rate     % receiver is capable of measuring rate
                err.psrr.clk = ((err_mdl(8,:,1) + err_prop(8,:,1)) + ...
                    (err_mdl(7,:,1)-err_mdl(7,:,2) + ...
                    err_prop(7,:,1)-err_prop(7,:,2)) / Tm) * obj.c * 1e3;   % s/s to mm/s
            end

            if nargout > 2
                % project uncertainty onto line-of-sight direction
                var.psr.prop = zeros(1,n);      % pseudorange uncertainty from nav
                var.psr.mdl = zeros(1,n);       % pseudorange uncertainty from model
                var.psrr.prop = zeros(1,n);     % Doppler uncertainty from nav
                var.psrr.mdl = zeros(1,n);      % Doppler uncertainty from model
                for i=1:n
                    var.psr.prop(i) = los(:,i,1)' * Pprop(1:3,1:3,i) * los(:,i,1) * psr_units^2;    % m^2
                    var.psr.mdl(i) = los(:,i,1)' * Pmdl(1:3,1:3,i) * los(:,i,1) * psr_units^2;      % m^2
                    var.psrr.prop(i) = los(:,i,1)' * Pprop(4:6,4:6,i) * los(:,i,1) * psrr_units^2;  % mm^2/s^2
                    var.psrr.mdl(i) = los(:,i,1)' * Pmdl(4:6,4:6,i) * los(:,i,1) * psrr_units^2;    % mm^2/s^2
                end
                % pseudorange uncertainty from clock model
                var.psr.clk = reshape(Pmdl(7,7,:) + Pprop(7,7,:), [1 n 1]) * obj.c^2;               % m^2
                % velocity uncertainty from clock model
                % Added frequency offset error as well!! We trust the LNSS
                % sats so when a receiver is comparing frequencies, it'll
                % assume its own is wrong. Thus the frequency offset is
                % mistakenly conferred onto the user or it throws off
                % velocity measurements. The allan variance is just noise
                % on top of frequency offset measurements, so since we are
                % directly attempting to measure the frequency that noise
                % is applied to our estimates + the offset!
                if rate
                    temp = obj.clock.dtcovariance(Tm);
                    var.psrr.clk = (reshape(Pmdl(8,8,:) + Pprop(8,8,:), [1 n 1]) + ...
                        ones(1,n) * temp(1,1) / Tm^2) * (obj.c * 1e3)^2;    % mm^2/s^2
                end
            end

            CN0 = obj.linkbudget(user, r);
            [erec, vrec] = user.rec.noise(CN0, ts, r, dr);

            % assign pseudorange outputs
            err.psr.rec = erec(1,:);
            var.psr.rec = vrec(1,:);
            err.psr.total = err.psr.prop + err.psr.mdl + err.psr.clk + err.psr.rec;
            true.psr = r;
            meas.psr = r + err.psr.total;

            if rate
                % assign pseudorange-rate outputs
                err.psrr.rec = erec(2,:);
                var.psrr.rec = vrec(2,:);
                err.psrr.total = err.psrr.prop + err.psrr.mdl + err.psrr.clk + err.psrr.rec;
                true.psrr = dr;
                meas.psrr = dr + err.psrr.total;
            end

            true.transmittime = tt;
            true.receivetime = ts;
            true.err = err;

            if nargout > 2
                var.psr.total = var.psr.prop + var.psr.mdl + var.psr.clk + var.psr.rec;
                if rate
                    var.psrr.total = var.psrr.prop + var.psrr.mdl + var.psrr.clk + var.psrr.rec;
                end
            end
        end

        function CN0 = linkbudget(obj,user,r)
            %LINKBUDGET Computes the receiver carrier to noise density ratio.
            %   Input:
            %    - user; User object instance
            %    - r; transmitter-receiver ranges (m)
            %   Output:
            %    - CN0; carrier to noise density ratio, dB-Hz
            arguments
                obj     (1,1)   NavSatellite
                user    (1,1)   User
                r       (1,:)   double {mustBePositive}
            end

            % link budget calculations to obtain C/N0
            freq = obj.ant.freq;
            Ad = 20 * log10((obj.c/freq)./(4*pi*r));    % dB, FSPL
            Ae = 0;                                     % dB, no atmospheric attenuation
            AP = obj.ant.P + obj.ant.gain + Ad + Ae;    % dBW, gain before receiver
            RP = AP + user.ant.gain + user.ant.As;      % dBW, gain before amps
            k = 1.3803e-23;                             % J/K, Boltzmann's constant
            N0 = 10*log10(k * user.ant.Ts);             % dBW/Hz, noise power spectral density
            CN0 = RP + user.ant.Nf + user.ant.L - N0;   % dB*Hz, carrier to noise density ratio
        end

        function [tt,r,dr] = timeofflight(obj,ts,user,tol)
            %TIMEOFFLIGHT Based on the provided receive times, find the
            %satellite transmit time and compute the range / range-rate.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - user; User module
            %    - tol; iteration tolerance for solving transmission time
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                user    (1,1)   User
                tol     (1,1)   double = 1e-9
            end

            n = length(ts);
            xuser = user.getstate(ts, 'J2000');
            xref0 = obj.getrefstates(ts, 'J2000');
            fref = obj.prop.statetotrajectory();
            fref = fref{1};
            % initial guess for transmit time is receive time
            tt = ts;
            % instantaneous range at measurement times, in m
            r = sqrt(sum((xref0(1:3,:) - xuser(1:3,:)).^2, 1)) * 1e3;
            % final s/c states at tt
            xref = zeros(6,n);
            % range-rate of measurements, in mm/s
            dr = zeros(1,n);

            for i=1:n
                rlast = r(i);

                for j=1:100
                    dt = rlast / obj.c;         % range to time-of-flight (s)
                    tj = ts(i) - dt;            % time offset guess
                    % updated range guess
                    xref(:,i) = fref(tj, 'J2000');
                    rj = norm(xref(1:3,i) - xuser(1:3,i)) * 1e3;

                    % if iteration is converging
                    if abs(rj - rlast) < tol
                        tt(i) = tj;
                        r(i) = rj;
                        break;
                    end

                    rlast = rj;                 % update iteration
                end
                if j == 100
                    error("timeofflight:notConverged", ...
                        "Time %d failed to converge in %d iterations.", i, j);
                end
                
                % compute range-rate by finding projection of relative
                % velocity along line-of-sight direction
                los = xref(1:3,i) - xuser(1:3,i);
                los = los / norm(los);
                vrel = xref(4:6,i) - xuser(4:6,i);
                dr(i) = vrel' * los * 1e6;      % convert km/s -> mm/s
            end
        end
    end

    methods (Static)
        function sandpile(t,data,labels,ytext,tlabel)
            %SANDPILE Plots a sandpile-looking breakdown of data on the
            %current figure.
            %   Input:
            %    - t; timestamps of data in seconds
            %    - data; rows are contributors, cols are timestamps t
            %    - labels; names of contributors (rows of data) for
            %              legend(), provide empty cell array {} to not display
            %    - ytext; text to display for ylabel
            %    - tlabel; (optional) include xlabel for time? Set to false
            %              if top plot of subplot, default true
            arguments
                t       (1,:)   double
                data    (:,:)   double {mustBeNonnegative}
                labels  (1,:)   {mustBeText}
                ytext   (1,:)   {mustBeText}
                tlabel  = true
            end

            colors = colororder;            % get colors for plotting
            c = size(colors,1);
            t = t - t(1);                   % normalize to 0
            units = "(s)";                  % default to seconds
            if t(end) >= 120 && t(end) < 120 * 60   % if 1min <= t < 120min, units are minutes
                t = t / 60;
                units = "(min)";
            elseif t(end) < 48 * 3600               % if 2hr <= t < 48hr, units are hours
                t = t / 3600;
                units = "(hrs)";
            elseif t(end) >= 2 * 86400              % if t >= 2d, units are days
                t = t / 86400;
                units = "(days)";
            end

            m = size(data,1);                       % number of contributors
            px = [t flip(t)];                       % x data
            data = [zeros(1,size(data,2)); data];   % add buffer row
            for i=1:m                               % plot patches
                data(i+1,:) = data(i+1,:) + data(i,:);  % add data
                py = [data(i,:) flip(data(i+1,:))];     % y data
                color = colors(mod(i-1,c)+1,:);         % get color
                patch(px, py, color, "EdgeColor", "none");
                if i == 1, hold on; end
            end

            % wrap up plot
            hold off; grid on;
            if tlabel, xlabel("Time " + units); end
            ylabel(ytext);
            if ~isempty(labels), legend(labels, "location", "best"); end              
        end
    end
end
