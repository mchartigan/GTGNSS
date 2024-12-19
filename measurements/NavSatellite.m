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
        filter
        % antenna object
        ant     (1,1)   TransmitAntenna
        % reference state info (avoids recalling runto() on prop and clock
        % if data has already been requested before
        tr      (1,:)   double = []
        xr      (9,:)   double
        frame_r (1,:)   char = ''
        fr      (1,1)   function_handle = @(~) 0
        % navigation state info (avoids rerunning filter if data has already 
        % been requested before)
        tn      (1,:)   double
        x0      (:,:)   double      % doubles as default navigation info
        P0      (9,9)   double      % doubles as default navigation info
        xn      (9,:)   double
        Pn      (9,9,:) double
        % should debug info be printed?
        DEBUG   (1,1)
    end

    properties (Access = private)
        % properties defining navigation message info (ephemeris and clock
        % models)
        tm_s    (1,:)   double      % model transition times for ephemeris
        % function handles for ephemeris models
        fm_s    (1,:)   cell    = {}
        tm_c    (1,:)   double      % model transition times for clock
        % function handles for clock models
        fm_c    (1,:)   cell    = {}
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
            %              options are: "EKF", "const"
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

            elseif strcmpi(filter, "const")
                % this option is for having a constant navigation uncertainty;
                % meas is only required to have a single property, P0, that
                % is a 9x9 covariance matrix
                if ~all(size(meas.P0) == [9 9])
                    error("NavSatellite:invalidMeas", ...
                        "For filter 'const', meas.P0 must be 9x9.");
                end
                obj.filter.P0 = meas.P0;

            else
                error("NavSatellite:invalidFilter", ...
                    'Filter must be "EKF". See documentation.');
            end

            if ~debug, warning('off', 'NavSatellite:debug'); end
            obj.DEBUG = debug;
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

            % avoid recalling if data has previously been generated
            if length(obj.tr) == length(ts) && ~any(obj.tr - ts) && strcmp(frame, obj.frame_r)
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
            obj.fr = @(t) ppval(spline(obj.tr, obj.xr), t);
        end

        function [xn,Pn] = getnavstates(obj,ts)
            %GETNAVSTATES Returns estimated state and covariance from the 
            %navigation filter at the provided times in the J2000 frame.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %   Output:
            %    - xn; state estimates at ts
            %    - Pn; state covariances at ts
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
            end

            if isa(obj.filter, 'struct')
                xn = obj.getrefstates(ts, 'J2000');
                % add noise to reference trajectory for nav states
                xn = mvnrnd(xn', obj.filter.P0)';
                Pn = repmat(obj.filter.P0, 1, 1, length(ts));
            else
                error("getnavstates:noImplementedError", ...
                    "Navigation filtering feature is not yet implemented.");

                % add ts to obj.filter.t with union() (maybe strip everything
                % before ts(1)?)
                % run filter, get data at time steps, and return it
            end

            % store data in case called again
            obj.tn = ts;
            obj.P0 = obj.filter.P0;
            obj.xn = xn;
            obj.Pn = Pn;
        end

        function [xm,xp,Pm,Pp,Po] = generatemodels(obj,ts,cadence)
            %GENERATEMODELS Returns model-derived spacecraft states (PVT).
            %This is the information a user would be able to get from data
            %in the nav message. It encompasses ephemeris and clock model
            %errors when differenced from getrefstates().
            %   Input:
            %    - ts; eval time steps, seconds past J2000
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
            %    - Po; covariance of nav soln at interval start
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                cadence (1,2)   double {mustBeNonnegative} = zeros(1,2);
            end

            n = length(ts);
            [xnav, Pnav] = obj.getnavstates(ts);

            I = NavSatellite.getupdateindices(ts, cadence);
            obj.tm_s = ts(I(1,:));
            obj.tm_c = ts(I(2,:));

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
            if nargout > 4, Po = zeros(9, 9, n); end

            for i=1:size(I, 2)-1        % for orbit
                int = I(1,i):I(1,i+1);                  % validity interval indices
                obj.prop.t0 = ts(int(1));               % new start time
                obj.prop.x0 = xnav(1:6, int(1));        % new start state
                tf = ts(int(end)) - obj.prop.t0;        % end of validity interval
                % ephemeris model, store in cell array
                obj.fm_s{i} = obj.prop.modelfit("Kepler", tf, 12);
                % assign model states relative to interval start
                xm(1:6,int) = obj.fm_s{i}(ts(int) - ts(int(1)));

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
                if nargout > 4
                    % determine error due to initial OD
                    Po(1:6,1:6,int) = repmat(Pnav(1:6,1:6,int(1)), 1, 1, length(int));
                end
            end

            for i=1:size(I, 2)-1        % for clock
                int = I(2,i):I(2,i+1);                  % validity interval indices
                obj.clock.t0 = ts(int(1));              % new start time
                obj.clock.x0 = xnav(7:9, int(1));       % new start state
                % clock model, store in cell array
                obj.fm_c{i} = obj.clock.modelfit();          
                % assign model states relative to interval start
                xm(7:9,int) = obj.fm_c{i}(ts(int) - ts(int(1)));     

                if nargout > 1
                    % assign propagated states model is based on (same for clock)
                    xp(7:9,int) = obj.fm_c{i}(ts(int) - ts(int(1)));
                end
                if nargout > 3
                    % propagate nav uncertainty while model is valid
                    Pp(7:9,7:9,int) = ...
                        obj.clock.proplyapunov(ts(int), Pnav(7:9,7:9,int(1)));
                end
                if nargout > 4
                    % determine error due to initial OD
                    Po(7:9,7:9,int) = repmat(Pnav(7:9,7:9,int(1)), 1, 1, length(int));
                end
            end

            % reassign propagator info after it's been changed
            obj.prop.t0 = t0p;
            obj.prop.x0 = x0p;
            obj.clock.t0 = t0c;
            obj.clock.x0 = x0c;
        end

        % TODO: implement getsmartmodelstates() where updates are adaptive

        function x = getmodelstates(obj,ts,frame)
            %GETMODELSTATES Given time(s) ts, return the model-computed
            %spacecraft state (ephemeris and time) as if from the
            %navigation message.
            %   Input:
            %    - ts; times past J2000 to get states at
            %    - frame; optional (default 'J2000'), frame to return data in
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                frame   (:,:)   {mustBeText} = 'J2000'
            end

            n = length(ts);
            x = zeros(9,n);

            % iterate over every time provided
            for i=1:n
                % iterate backwards through piecewise timesteps for ephemeris
                for j=length(obj.tm_s)-1:-1:1
                    % first start time that is <= given time, use that and break
                    if ts(i) >= obj.tm_s(j)
                        x(1:6,i) = cspice_sxform('J2000', frame, ts(i)) * ...
                                   obj.fm_s{j}(ts(i) - obj.tm_s(j));
                        break;
                    end
                end

                % iterate through piecewise timesteps for clock
                for j=length(obj.tm_c)-1:-1:1
                    % first start time that is <= given time, use that and break
                    if ts(i) >= obj.tm_c(j)
                        x(7:9,i) = obj.fm_c{j}(ts(i) - obj.tm_c(j));
                        break;
                    end

                end
            end
        end

        function [meas,var,true] = observemeas(obj,ts,user,opts)
            %OBSERVEMEAS Finds the real pseudorange and pseudorange-rate
            %(Doppler) measurements at times ts between NavSatellite and
            %the user.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - user; User object instance
            %    - frame; reference frame user data is provided in
            %    - opts; s/c settings struct all fields must be set:
            %             - cadence (see generatemodels input)
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
            % run propagator to get trajectory estimate
            obj.prop.runat(ts - obj.prop.t0, 'J2000');
            fref = obj.prop.statetotrajectory();
            [tt, r, dr] = obj.timeofflight(ts, fref{1}, user);

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
                if ts(1)-Tm < obj.prop.t0 || ts(1)-Tm < obj.clock.t0
                    error("observemeas:invalidStartMeas", ...
                        "First measurement must be >= User.rec.Tm past starting epoch.");
                end

                % interleave times
                tt_old = tt;
                tt = union(tt_old-Tm, tt_old);
                it = ismember(tt, tt_old);
                it_Tm = ismember(tt, tt_old-Tm);
                ts_old = ts;
                ts = union(ts_old-Tm, ts_old);
                is = ismember(ts, ts_old);
                is_Tm = ismember(ts, ts_old-Tm);
            end

            % get user state info for measurements
            xref = obj.getrefstates(tt, 'J2000');
            xuser = user.getstates(ts, 'J2000');
            if nargout > 2
                [xmdl,xprop,Pmdl,Pprop] = obj.generatemodels(tt, opts.cadence);
            else
                [xmdl,xprop] = obj.generatemodels(tt, opts.cadence);
            end

            if rate
                % get nav message update intervals as they were assigned in
                % .generatemodels()
                mask2 = ones(2,n);
                I = NavSatellite.getupdateindices(tt,opts.cadence);
                tI =  tt(union(I(1,:), I(2,:)));

                % un-interleave times
                tt = tt(it);
                ts = ts(is);

                % with measurement transmit times, compute which need to be
                % masked based on if they were within Tm of a nav message
                % update
                mat = tt - tI';
                mask2(2,:) = and(mask2(2,:), ~any(and(mat >= 0, mat <= Tm)));

                % separate states
                % store concatenated data and create 2-page 3D matrices, then
                % un-interleave data, storing in 2 separate pages of matrix
                % first tt(i)-Tm comes, then tt(i), so store first value on 
                % second page and second on first page
                xref = cat(3, xref(:,it), xref(:,it_Tm));
                xuser = cat(3, xuser(:,is), xuser(:,is_Tm));
                xmdl = cat(3, xmdl(:,it), xmdl(:,it_Tm));
                xprop = cat(3, xprop(:,it), xprop(:,it_Tm));
                % remove tt-Tm data from covariances
                if nargout > 2
                    Pmdl = Pmdl(:,:,it);
                    Pprop = Pprop(:,:,it);
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
                    err.psrr.prop(i) = err_prop(4:6,i,1)'*los(:,i,1) * psrr_units;
                    err.psrr.mdl(i) = err_mdl(4:6,i,1)'*los(:,i,1) * psrr_units;
                end
            end
            % pseudorange error from clock model
            err.psr.clk = (err_mdl(7,:,1) + err_prop(7,:,1)) * obj.c;       % s to m
            % velocity error from clock model
            if rate     % receiver is capable of measuring rate
                err.psrr.clk = (err_mdl(7,:,1)-err_mdl(7,:,2) + ...
                    err_prop(7,:,1)-err_prop(7,:,2)) / Tm * obj.c * 1e3;    % s/s to mm/s
                % err.psrr.clk = ((err_mdl(8,:,1) + err_prop(8,:,1))) * obj.c * 1e3;   % s/s to mm/s
            end

            if nargout > 2
                % project uncertainty onto line-of-sight direction
                var.psr.prop = zeros(1,n);      % pseudorange uncertainty from nav
                var.psr.mdl = zeros(1,n);       % pseudorange uncertainty from model
                var.psrr.prop = zeros(1,n);     % Doppler uncertainty from nav
                var.psrr.mdl = zeros(1,n);      % Doppler uncertainty from model
                for i=1:n
                    % pseudorange uncertainty from ephemeris, in m^2
                    var.psr.prop(i) = los(:,i,1)' * Pprop(1:3,1:3,i) * los(:,i,1) * psr_units^2;
                    var.psr.mdl(i) = los(:,i,1)' * Pmdl(1:3,1:3,i) * los(:,i,1) * psr_units^2;

                    if rate
                        % Doppler uncertainty from ephemeris, in mm^2/s^2
                        var.psrr.prop(i) = los(:,i,1)' * Pprop(4:6,4:6,i) * los(:,i,1) * psrr_units^2;
                        var.psrr.mdl(i) = los(:,i,1)' * Pmdl(4:6,4:6,i) * los(:,i,1) * psrr_units^2;
                    end
                end
                % pseudorange uncertainty from clock model, in m^2
                var.psr.clk_offset = reshape(Pmdl(7,7,:) + Pprop(7,7,:), [1 n 1]) * obj.c^2;
                % implement phase noise here if ya want
                var.psr.clk_stability = zeros(size(var.psr.clk_offset));
                var.psr.clk = var.psr.clk_offset + var.psr.clk_stability;

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
                    % Doppler uncertainty from s/c clock, in mm^2/s^2
                    var.psrr.clk_offset = reshape(Pmdl(8,8,:) + Pprop(8,8,:), [1 n 1]) * (obj.c * 1e3)^2;
                    var.psrr.clk_stability = ones(1,n) * obj.clock.stability(Tm) * (obj.c * 1e3)^2;
                    var.psrr.clk = var.psrr.clk_offset + var.psrr.clk_stability;
                end
            end

            % compute link budget and receiver noise
            CN0 = obj.linkbudget(user, r);
            % add mask to link budget %
            % get nadir direction at each time step
            nadir = xuser(1:3,:,1) ./ sqrt(sum(xuser(1:3,:,1).^2, 1));
            % compute angle between nadir and satellite
            tosat = acos(sum(nadir .* los(:,:,1), 1));
            % -300 dB/Hz if angle is over off-boresight mask angle
            CN0 = CN0 - 300 * (tosat > pi/2 - user.ant.mask);
            [erec, vrec, mask1] = user.rec.noise(CN0, ts, r, dr);

            % build and assign mask
            if rate
                meas.mask = and(mask1, mask2);
            else
                meas.mask = mask1;
            end

            % plot CN0 for simulation
            if obj.DEBUG
                figure();
                plotformat("APA", 0.25, "scaling", 1, "coloring", "science");
                elev = (pi/2 - tosat) * 180/pi;
                plot(elev(meas.mask(1,:)), CN0(meas.mask(1,:)), "LineWidth", 2);
                grid on;
                xlabel("Elevation angle (\circ)");
                ylabel("C/N0 (dB-Hz)");
                title("C/N0 vs. User Antenna Elevation Angle");
            end

            % assign pseudorange outputs
            err.psr.rec = erec(1,:);
            var.psr.rec_thermal = vrec.thermal(1,:);
            var.psr.rec_clk = vrec.clk(1,:);
            var.psr.rec_dyn = vrec.dyn(1,:);
            var.psr.rec = vrec.total(1,:);
            err.psr.total = err.psr.prop + err.psr.mdl + err.psr.clk + err.psr.rec;
            true.psr = r + obj.c * xuser(7,:,1);
            meas.psr = true.psr + err.psr.total;

            if rate
                % assign pseudorange-rate outputs
                err.psrr.rec = erec(2,:);
                var.psrr.rec_thermal = vrec.thermal(2,:);
                var.psrr.rec_clk = vrec.clk(2,:);
                var.psrr.rec_dyn = vrec.dyn(2,:);
                var.psrr.rec = vrec.total(2,:);
                err.psrr.total = err.psrr.prop + err.psrr.mdl + err.psrr.clk + err.psrr.rec;
                true.psrr = dr + obj.c * xuser(8,:,1);
                meas.psrr = true.psrr + err.psrr.total;
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

        function y = computemeas(obj,t,x,user)
            %COMPUTEMEAS Computes the ideal pseudorange and Doppler
            %measurements at time t between the user and NavSatellite
            %   Input:
            %    - t; time, seconds past J2000
            %    - x; state of USER at time t
            %    - user; User object
            arguments
                obj     (1,1)   NavSatellite
                t       (1,:)   double
                x       (9,:)   double
                user    (1,1)   User
            end

            r_u = x(1:3);           % user position
            v_u = x(4:6);           % user velocity

            % find transmission time w.r.t meas, based on nav msg knowledge
            % of satellite states
            state2eph = @(x) x(1:6,:);
            modeleph = @(t,~) state2eph(obj.getmodelstates(t));
            tt = obj.timeofflight(t, modeleph, user);
            x_s = obj.getmodelstates(tt, user.frame);
            r_s = x_s(1:3);
            v_s = x_s(4:6);

            dr = r_s - r_u;         % relative user->sat position
            dv = v_s - v_u;         % relative user->sat velocity
            rho = norm(dr);         % scalar range
            dvdr = dv'*dr;          % dot product of dv and dr
            y = [rho      + obj.c * (x(7) - x_s(7))
                 dvdr/rho + obj.c * (x(8) - x_s(8))];
        end

        function H = measpartials(obj,t,x,user)
            %MEASPARTIALS Computes the partial derivative of the measurement
            %model w.r.t. x at the current state.
            %   Inputs:
            %    - t; time of measurement
            %    - x; state of user
            %    - user; User object
            arguments
                obj     (1,1)   NavSatellite
                t       (1,1)   double
                x       (9,1)   double
                user    (1,1)   User
            end

            r_u = x(1:3);           % user position
            v_u = x(4:6);           % user velocity

            % find transmission time w.r.t meas, based on nav msg knowledge
            % of satellite states
            state2eph = @(x) x(1:6,:);
            modeleph = @(t,~) state2eph(obj.getmodelstates(t));
            tt = obj.timeofflight(t, modeleph, user);
            x_s = obj.getmodelstates(tt, user.frame);
            r_s = x_s(1:3);
            v_s = x_s(4:6);

            dr = r_s - r_u;         % relative user->sat position
            dv = v_s - v_u;         % relative user->sat velocity
            rho = norm(dr);         % scalar range
            dvdr = dv'*dr;          % dot product of dv and dr
            H = [-dr'/rho                     0 0 0    obj.c 0     0
                  (dr'*dvdr/rho^3 - dv'/rho) -dr'/rho  0     obj.c 0];
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

        function [tt,r,dr] = timeofflight(obj,ts,fref,user,tol)
            %TIMEOFFLIGHT Based on the provided receive times, find the
            %satellite transmit time and compute the range / range-rate.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - fref; handle to get NavSatellite reference state (valid
            %            over ts, +/- a couple seconds)
            %    - user; User module
            %    - tol; iteration tolerance for solving transmission time
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                fref    (1,1)   function_handle
                user    (1,1)   User
                tol     (1,1)   double = 1e-9
            end

            n = length(ts);
            xuser = user.getstates(ts, 'J2000');
            xref0 = fref(ts(1), 'J2000');

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
                    warning("timeofflight:notConverged", ...
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

    methods (Static, Access = private)
        function I = getupdateindices(ts,cadence)
            %GETUPDATEINDICES Returns indices of ts where ephemeris and
            %clock updates should happen, based on cadence.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - cadence; # of seconds after which models are updated, can be
            %               separate for ephemeris and clock -- 0 means don't
            %               update
            arguments
                ts      (1,:)   double
                cadence (1,2)   double {mustBeNonnegative} = zeros(1,2);
            end

            dt = ts(end) - ts(1);
            n = length(ts);
            I = ones(2, ceil(dt / cadence(1)) + 1);
            if any(cadence)
                t_prop = ts(1);
                ii = 2;
                t_clk = ts(1);
                jj = 2;
                for i=2:n
                    if ts(i) >= t_prop + cadence(1)
                        I(1,ii) = i;
                        t_prop = ts(i);
                        ii = ii + 1;
                    end
                    if ts(i) >= t_clk + cadence(end)
                        I(2,jj) = i;
                        t_clk = ts(i);
                        jj = jj + 1;
                    end
                end
            end

            I(:,end) = n;
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
