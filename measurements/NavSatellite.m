classdef NavSatellite < handle
    %NavSatellite Class for describing the properties and trajectory of a
    %satellite that is part of a radionavigation satellite system (whether
    %that's earth-based GNSS or lunar).
    
    properties
        % [satellite reference trajectory; clock reference trajectory]
        traj    (2,1)   Trajectory
        % propagator instance (default one so MATLAB doesn't throw a fit)
        prop    (1,1)   SatellitePropagator = SatellitePropagator(OrbitPropagator(1),Clock("none", zeros(4,1)) )
        % nav filter
        filter
        % antenna object
        ant     (1,1)   TransmitAntenna
        % nav message update cadence info
        cadence (1,1)   double {mustBeNonnegative}
        % nav message coefficient count (per axis, so will be 3x)
        ncoef   (1,1)   {mustBeNonnegative,mustBeInteger} = 6
        % reference state info (avoids recalling runto() on prop and clock
        % if data has already been requested before
        tr      (1,:)   double = []
        xr      (9,:)   double
        frame_r (1,:)   char = ''
        % navigation state info (avoids rerunning filter if data has already 
        % been requested before)
        tn      (1,:)   double
        x0      (:,:)   double      % doubles as default navigation info
        P0      (9,9)   double      % doubles as default navigation info
        xn      (9,:)   double
        Pn      (9,9,:) double
        % satellite ID number
        ID      (1,1)   {mustBeInteger,mustBeNonnegative}
        % should debug info be printed?
        DEBUG   (1,1)
        % various constant error variances
        group   (1,1)   double {mustBeNonnegative} = 0.1^2  % m^2
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
        % speed of light (m/s, m^2/s^2)
        c       (1,1)   double = 299792458
        c2      (1,1)   double = 299792458^2
        % speed of light (km/s, km^2/s^2)
        c_km    (1,1)   double = 299792.458
        c_km2   (1,1)   double = 299792.458^2
    end
    
    methods
        function obj = NavSatellite(traj,prop,filter,meas,ant,options)
            %NAVSATELLITE Construct a NavSatellite instance.
            %   Inputs:
            %    - traj (2,1) Trajectory; traj(1) is orbit trajectory,
            %       traj(2) is clock
            %    - prop; propagator instance for generating future trajectories
            %       (SatellitePropagator)
            %    - filter; navigation filter type for finding s/c uncertainty,
            %       options are: "EKF", "const"
            %    - meas; struct containing measurement info for the chosen
            %       filter -- names should correspond to filter constructor
            %       arguments (may also include varargin)
            %    - ant; TransmitAntenna object containing sat info
            %    - ID; satellite ID (should be unique)
            %    - cadence; navigation message update rate (in s)
            %    - debug; should debug warnings be printed? true/false
            arguments
                traj    (2,1)   Trajectory = Trajectory()
                prop    (1,1)   SatellitePropagator = SatellitePropagator()
                filter  (1,:)   {mustBeText} = "none"
                meas    (1,1)   struct = struct()
                ant     (1,1)   TransmitAntenna = TransmitAntenna()
                options.ID      (1,1)   {mustBeNonnegative,mustBeInteger} = 0
                options.cadence (1,1)   {mustBePositive} = 7200
                options.debug   (1,1)   double = false
            end
            % permit empty instantiation

            % assign instances
            obj.traj = traj;
            obj.prop = prop;
            obj.ant = ant;
            % initialize filter for default case of "none"
            obj.filter = [];

            if strcmpi(filter, "const")
                % this option is for having a constant navigation uncertainty;
                % meas is only required to have a single property, P0, that
                % is a 9x9 covariance matrix
                if ~all(size(meas.P0) == [9 9])
                    error("NavSatellite:invalidMeas", ...
                        "For filter 'const', meas.P0 must be 9x9.");
                end
                obj.filter.P0 = meas.P0;

            elseif ~strcmpi(filter, "none")
                error("NavSatellite:invalidFilter", ...
                    'Filter must be "EKF". See documentation.');
            end

            % assign options attributes
            obj.cadence = options.cadence;
            obj.ID = options.ID;
            if ~options.debug, warning('off', 'NavSatellite:debug'); end
            obj.DEBUG = options.debug;
        end

        function [T,dT,CN0,msg,err,var] = transmitsignal(obj,ts,user)
            %TRANSMITSIGNAL Computes the true transmit time (s) and Doppler
            %shift (s/s) between the satellite and the user. Navigation message
            %data necessary to reconstruct measurements is generated. Errors are
            %decomposed by source and provided as additional output.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - user; User object instance
            %    - frame; reference frame user data is provided in
            %    - opts; settings struct, fields include:
            %       - SISE; "NASA" or "custom"
            %       - cadence; see generatemodels input, only used if SISE
            %          is "custom"
            %   Output:
            %    - T; transmitter-receiver delay (s)
            %    - dT; transmitter-receiver Doppler (s/s)
            %    - CN0; carrier to noise density ratio, dB-Hz
            %    - msg; struct containing navigation message data
            %    - err; error applied to T and dT
            %    - var; variance of error err
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                user    (1,1)   User
            end

            % create copy of user (shallow, so referenced objects are same unfort)
            olduser = user;
            user = copy(user);
            % adjust user to be relative to obj.prop.body
            xtemp = user.motion.xs;
            xtemp(1:6,:) = xtemp(1:6,:) + cspice_spkezr(user.body, ...
                user.motion.ts, user.motion.frame, 'NONE', obj.prop.body) * 1e3;
            user.motion = Trajectory(user.motion.ts, xtemp, user.motion.frame);

            % get true measurements
            % run propagator to get trajectory estimate
            [tt, r, dr] = obj.timeofflight(ts, user);

            % SIGNAL IN SPACE ERROR CALCULATION %
            [err,var,msg,los] = obj.getSISE(tt,ts,user);

            % OUTPUT FORMATTING %
            % add applicable error to delay and Doppler
            xc = obj.traj(2).get(tt);
            T  = r  - xc(1,:) - err.clk_prop(1,:);
            dT = dr - xc(2,:) - err.clk_prop(2,:);

            % ANTENNA GAIN %
            % compute transmitter angle %
            % state of user w.r.t. obj.prop.body
            xuser = user.motion.getpos(ts, 'J2000');
            % state of sat w.r.t. obj.prop.body
            xsat  = obj.traj(1).get(tt, 'J2000');
            % get User->obj.prop.body direction at each time step
            u_u1 = -xuser ./ sqrt(sum(xuser.^2, 1));
            % get nadir direction at sat at each time step
            u_s1 = -xsat(1:3,:) ./ sqrt(sum(xsat(1:3,:).^2, 1));
            % compute angle between nadir and user
            touser = acos(sum(u_s1 .* -los, 1));
            % compute angle between User->obj.prop.body dir and satellite
            tosat = acos(sum(u_u1 .* los, 1));
            % determine received power at user antenna
            AP = obj.txlinkbudget(r,touser);
            CN0 = user.ant.getCN0(AP,tosat);

            % PLANET INTERSECTION CALCULATION %
            % primary body
            R1 = cspice_bodvrd(obj.prop.body, 'RADII', 3) * 1e3;
            R1 = max(R1);       % radius of obj.prop.body, m
            
            for i=1:length(ts)
                x_s = xsat(1:3,i);                  % Moon -> sat
                r_s = norm(x_s);                    % || Moon -> sat ||
                r_su = norm(xuser(1:3,i) - x_s);    % || sat -> User ||
                u_us = los(1:3,i);                  % User -> sat
                a_1st = asin(R1/r_s);               % < body-sat-bodyTangent angle
                a_1su = acos(u_us' * x_s / r_s);    % < body-sat-User angle
                r_t = sqrt(r_s^2 - R1^2);           % || sat -> bodyTangent ||
                % if the body center/body tangent angle from the sat
                % POV is bigger than the body center/user angle and the
                % range is > body tangent range, sat is out of view.
                if a_1su < a_1st && r_su > r_t
                    T(i)  = NaN;
                    dT(i) = NaN;
                    CN0(i) = CN0(i) - 300;
                end
            end

            % evaluate second body intersection if user and satellite
            % aren't around the same central body
            if ~strcmpi(obj.prop.body, user.body)
                % user trajectory relative to its own central body
                xuser = olduser.getstates(ts, 'J2000');
                % secondary body
                R2 = cspice_bodvrd(user.body, 'RADII', 3) * 1e3;
                R2 = max(R2);       % radius of user.body, m

                for i=1:length(ts)
                    r_u2 = norm(xuser(1:3,i));      % || User -> body ||
                    a_2ut = asin(R2/r_u2);          % < body-User-bodyTangent angle
                    % < body-User-sat angle
                    a_2us = acos(los(1:3,i)' * -xuser(1:3,i) / r_u2);
                    % if the body center/body tangent angle from the user POV
                    % is bigger than the body center/sat angle, sat is out of 
                    % view. Assumes sat is further away than body.
                    if a_2us < a_2ut
                        T(i)  = NaN;
                        dT(i) = NaN;
                        CN0(i) = CN0(i) - 300;
                    end
                end
            end
        end

        function [T,dT,CN0,msg,err,var] = transmitearthsignal(obj,ts,user)
            %TRANSMITEARTHSIGNAL Computes the true transmit time (s) and Doppler
            %shift (s/s) between the satellite and the user. Navigation message
            %data necessary to reconstruct measurements is generated. Errors are
            %decomposed by source and provided as additional output.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - user; User object instance
            %    - frame; reference frame user data is provided in
            %    - opts; settings struct, fields include:
            %       - SISE; "NASA" or "custom"
            %       - cadence; see generatemodels input, only used if SISE
            %          is "custom"
            %   Output:
            %    - T; transmitter-receiver delay (s)
            %    - dT; transmitter-receiver Doppler (s/s)
            %    - CN0; carrier to noise density ratio, dB-Hz
            %    - msg; struct containing navigation message data
            %    - err; error applied to T and dT
            %    - var; variance of error err
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                user    (1,1)   User
            end

            % get true measurements
            % run propagator to get trajectory estimate
            [tt, r, dr] = obj.timeofflight(ts, user);

            % SIGNAL IN SPACE ERROR CALCULATION %
            [err,var,msg,los] = obj.getSISE(tt,ts,user);

            % OUTPUT FORMATTING %
            % add applicable error to delay and Doppler
            xc = obj.traj(2).get(tt);
            T  = r  - xc(1,:) - err.clk_prop(1,:);
            dT = dr - xc(2,:) - err.clk_prop(2,:);

            

            % PLANET INTERSECTION CALCULATION %
            R = cspice_bodvrd('MOON', 'RADII', 3) * 1e3;
            R = R(1);           % radius of moon
            Re = cspice_bodvrd('EARTH', 'RADII', 3) * 1e3;
            Re = max(Re) + 1e6;
            for i=1:length(ts)
                x_s = xsat(1:3,i);                  % Moon -> GNSS
                r_s = norm(x_s);                    % || Moon -> GNSS ||
                r_su = norm(xuser(1:3,i) - x_s);    % || GNSS -> User ||
                u_us = los(1:3,i);                  % User -> GNSS
                a_mst = asin(R/r_s);                % < Moon-GNSS-MoonTangent angle
                a_msu = acos(u_us' * x_s / r_s);    % < Moon-GNSS-User angle
                r_t = sqrt(r_s^2 - R^2);            % || GNSS -> MoonTangent ||
                % if the moon center/moon tangent angle from the GNSS
                % POV is bigger than the moon center/user angle and the
                % range is > moon tangent range, GNSS is out of view.
                if a_msu < a_mst && r_su > r_t
                    T(i)  = NaN;
                    dT(i) = NaN;
                    CN0(i) = CN0(i) - 300;
                end

                r_ue = norm(x_ue(:,i));         % || User -> Earth ||
                a_eut = asin(Re/r_ue);          % < Earth-User-EarthTangent angle
                % < Earth-User-GNSS angle
                a_eus = acos(u_us' * x_ue(:,i) / r_ue);
                r_t2 = sqrt(r_ue^2 - Re^2);     % || User -> EarthTangent ||
                % if the earth center/earth tangent angle from the user POV
                % is bigger than the earth center/GNSS angle and the range
                % is > earth tangent range, GNSS is out of view
                if a_eus < a_eut && r_su > r_t2
                    T(i)  = NaN;
                    dT(i) = NaN;
                    CN0(i) = CN0(i) - 300;
                end
            end

            figure();
            plotformat("APA", 0.6);
            yyaxis left;
            tplot = (ts - ts(1)) / 60;
            plot(tplot, touser * 180/pi);
            hold on;
            plot(tplot, tosat * 180/pi);
            hold off;
            xlabel("Time (min)");
            ylabel("Angle (deg)");
            legend(["GNSS Boresight", "User Boresight"], location="best");

            yyaxis right;
            plot(tplot, CN0);
            ylabel("C/N0 (dB-Hz)");
        end

        function [err,var,msg,los,bias] = getSISE(obj,tt,ts,user)
            %GETSISE Return the signal in space error for the satellite to
            %the given user.
            %   SISE Position: 13.43 m 3-sigma
            %   SISE Velocity: 1.2 mm/s 3-sigma @ 10s
            %
            %   References:
            %    - Speciale, N., Lunar Relay Services Requirements Document 
            %       (SRD), ESC-LCRNS-REQ-0090, NASA.
            %
            %   Input:
            %    - tt; signal transmission times, seconds past J2000
            %    - ts; signal reception times, seconds past J2000
            %    - user; User object instance
            arguments
                obj     (1,1)   NavSatellite
                tt      (1,:)   double
                ts      (1,:)   double
                user    (1,1)   User
            end

            n = length(ts);     % no. of measurements
            frame = 'J2000';
            % get reference trajectory of satellite and clock
            xref = zeros(9,n);
            xref(1:6,:) = obj.traj(1).get(tt, frame);
            xref(7:9,:) = obj.traj(2).get(tt);
            % get reference trajectory of user
            xuser = user.getstates(ts, frame);
            
            % create line-of-sight direction
            los = xref(1:3,:,1) - xuser(1:3,:,1);
            los = los ./ sqrt(sum(los.^2, 1));

            % generate nav update times
            % start whenever the satellite starts :)
            tmsg = obj.traj(1).ts(1):obj.cadence:tt(end);
            if tmsg(end) ~= tt(end), tmsg = [tmsg tt(end)]; end
            [xmsg, Pmsg] = obj.getnavstates(tmsg, frame);
            msg = zeros(length(tmsg)-1, 23+3*obj.ncoef);
            Pnav = zeros(9,9,n);        % store starting uncertainty
            xprop = zeros(9,n);         % store the propagated states
            Pprop = zeros(9,9,n);
            xmdl = zeros(9,n);

            % compute errors and variances
            err.eph_prop = zeros(2,n);      % error due to initial OD and propagation
            err.eph_mdl  = zeros(2,n);      % error due to ephemeris parameterization
            err.clk_prop = zeros(2,n);      % error due to clock est. and propagation
            err.clk_mdl  = zeros(2,n);      % error due to clock parameterization
            err.group    = zeros(2,n);      % residual of group delay calibration error
            var.eph_prop = zeros(2,n);      % variance of initial OD and propagation
            var.eph_mdl  = zeros(2,n);      % variance of ephemeris parameterization
            var.clk_prop = zeros(2,n);      % variance of clock est. and propagation
            var.clk_mdl  = zeros(2,n);      % variance of clock parameterization
            var.group    = zeros(2,n);      % variance of group delay calibration error
            var.phase    = zeros(2,n);      % variance due to clock phase noise
            var.ODTS     = zeros(2,n);      % variance of initial ODTS

            % iterate over update times
            for i=1:length(tmsg)-1
                % nav states are applicable starting at tt(1)-1, so this logic
                % should cover all tt
                jj = and(tt > tmsg(i), tt <= tmsg(i+1));
                % % store nav uncertainty
                Pnav(:,:,jj) = repmat(Pmsg(:,:,i), 1, 1, sum(jj));
                % propagate states over given times. tmsg(i) provided so
                % trajectory starts at appropriate time
                tsub = [tmsg(i) tmsg(i+1)];
                [ts,xsub] = obj.prop.run(tsub, xmsg(:,i), 1000, frame, false);
                % provide these propagated states (plus initial one) as a
                % Trajectory and create a navigation message about it.
                subeph = Trajectory(ts, xsub(1:6,:), frame);
                subclk = Trajectory(ts, xsub(7:9,:));
                msg(i,:) = obj.generatenavmsg([subeph; subclk]);

                if sum(jj)
                    % cut out xmsg(:,i) since it may not align with tt
                    xprop(:,jj) = [subeph.get(tt(jj), frame); subclk.get(tt(jj))];
                    Ptemp = obj.prop.proplyapunov([ts(1) tt(jj)], xsub(:,1), Pmsg(:,:,i));
                    Pprop(:,:,jj) = Ptemp(:,:,2:end);
    
                    % compute model states and all errors/variances
                    for k=find(jj)
                        [xmdl(:,k),T] = RadiometricObsSim.geteph(tt(k), ...
                            obj.ID, msg(i,:), obj.prop.orbit.pri.GM);
    
                        % rotate to inertial since that's where we're handling
                        xmdl(1:6,k) = T \ xmdl(1:6,k);
                        % compute time step errors
                        err_prop = xprop(:,k) - xref(:,k);
                        err_mdl = xmdl(:,k) - xprop(:,k);
    
                        % range and range-rate error (in m and m/s) due to propagation
                        err.eph_prop(1,k) = err_prop(1:3)' * los(:,k); % / obj.c_km;
                        err.eph_prop(2,k) = err_prop(4:6)' * los(:,k); % / obj.c_km;
                        % " due to ephemeris model
                        err.eph_mdl(1,k) = err_mdl(1:3)' * los(:,k); % / obj.c_km;
                        err.eph_mdl(2,k) = err_mdl(4:6)' * los(:,k); % / obj.c_km;
                        % " due to onboard clock offset from proper time
                        % (negative to account for how it impacts the measurement)
                        err.clk_prop(:,k) = err_prop(7:8);
                        % " due to clock model
                        err.clk_mdl(:,k) = err_mdl(7:8);
                        % variance from OD and state propagation (in m^2 and m^2/s^2)
                        var.eph_prop(1,k) = los(:,k)' * Pprop(1:3,1:3,k) * los(:,k); % / obj.c_km2;
                        var.eph_prop(2,k) = los(:,k)' * Pprop(4:6,4:6,k) * los(:,k); % / obj.c_km2;
                        var.clk_prop(:,k) = diag(Pprop(7:8,7:8,k));
                        % " from ODTS
                        var.ODTS(1,k) = los(:,k)' * Pnav(1:3,1:3,k) * los(:,k) + Pnav(7,7,k);
                        var.ODTS(2,k) = los(:,k)' * Pnav(4:6,4:6,k) * los(:,k) + Pnav(8,8,k);
                    end
    
                    % " from parameterization (computed in batch)
                    var.eph_mdl(:,jj) = repmat(sum(err.eph_mdl(:,jj).^2, 2)/(sum(jj)-1), 1, sum(jj));
                    var.clk_mdl(:,jj) = repmat(sum(err.clk_mdl(:,jj).^2, 2)/(sum(jj)-1), 1, sum(jj));
                end
            end

            % work group delays (calibration error, not evolving over time)
            var.group(1,:) = obj.group;
            err.group(1,:) = mvnrnd(0, obj.group);
            % phase noise and frequency stability (already in error but not
            % variance budget)
            % get loop bandwidth
            if user.rx.PLL, Bn = user.rx.Bn_PLL;
            elseif user.rx.FLL, Bn = user.rx.Bn_FLL;
            else, Bn = user.rx.Bn;
            end
            % jitter, in rad^2
            [~,var.phase(1,:)] = obj.prop.clock.getjitter(user.rx.freq,Bn);
            % convert to m^2
            var.phase(1,:) = var.phase(1,:) * (2*pi*user.rx.freq)^(-2) * obj.c^2;
            % (s/s)^2 to (m/s)^2
            var.phase(2,:) = obj.prop.clock.stability(user.rx.T_FLL) * obj.c^2;

            % Apply only the errors that will occur due to signal
            % transmission. We're making this realistic here!
            % (i.e. let user make the modeling errors themselves)
            err.total = -err.clk_prop + err.eph_mdl + err.eph_prop + ...
                        err.clk_mdl + err.group;
            % Keep all the variance terms just for budgeting
            var.total = var.eph_prop + var.eph_mdl + var.clk_prop + ...
                        var.clk_mdl + var.group + var.phase;
        end

        function [xn,Pn] = getnavstates(obj,ts,frame)
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
                frame   (1,:)   {mustBeText} = 'J2000'
            end

            xn = zeros(9, length(ts));
            if isa(obj.filter, 'struct')
                xn(1:6,:) = obj.traj(1).get(ts, frame);
                xn(7:9,:) = obj.traj(2).get(ts);
                % add noise to reference trajectory for nav states
                xn = mvnrnd(xn', obj.filter.P0)';
                Pn = repmat(obj.filter.P0, 1, 1, length(ts));

                % store data in case called again
                obj.tn = ts;
                obj.P0 = obj.filter.P0;
                obj.xn = xn;
                obj.Pn = Pn;

            elseif isempty(obj.filter)
                % no filter at all, return truth states
                xn(1:6,:) = obj.traj(1).get(ts, frame);
                xn(7:9,:) = obj.traj(2).get(ts);
                Pn = repmat(zeros(9,9), 1, 1, length(ts));

                % store data in case called again
                obj.tn = ts;
                obj.P0 = Pn(:,:,1);
                obj.xn = xn;
                obj.Pn = Pn;
            else
                error("getnavstates:noImplementedError", ...
                    "Navigation filtering feature is not yet implemented.");

                % add ts to obj.filter.t with union() (maybe strip everything
                % before ts(1)?)
                % run filter, get data at time steps, and return it
            end
        end

        function msg = generatenavmsg(obj,traj)
            %GENERATENAVMSG Returns all coefficients and info needed to
            %compute satellite ephemeris and clock offsets from Trajectory
            %instances provided.
            %   MESSAGE COLUMN INDICES
            %    1   2   3     4     5     6    7    8    9     10    11 12 13 ...
            %    ts  ID  IODC  IODE  t_oc  af0  af1  af2  T_GD  t_oe  A  e  i0 ...
            %    14    15  16  17    18       19  
            %    RAAN0 w0  M0  idot  RAANdot  wdot
            %   Input:
            %    - traj (2,1) Trajectory; 2-element Trajectory vector, first
            %       being ephemeris and second is clock
            %   Output:
            %    - msg; Formatted array containing info to compute
            %       ephemeris and clock at any time
            arguments
                obj     (1,1)   NavSatellite
                traj    (2,1)   Trajectory
            end

            % get navigation states at update times
            traj_eph = traj(1);
            traj_clk = traj(2);
            t0 = traj_eph.t0;

            % known fixed values
            T_GD = 0;       % group delay offset b/n broadcast freqs (normal group delay can be added to af0)

            % update ephemeris
            eph = obj.prop.orbit.AFSfit(traj_eph, obj.ncoef);
            IODE = mod(floor(eph.t_oe), 1024);      % simple hash

            % update clock offsets
            t_oc = t0;
            IODC = mod(floor(t_oc), 1024);          % simple hash
            [~,D] = obj.prop.clock.modelfit(traj_clk, t0);
            af0 = D(1);
            af1 = D(2);
            af2 = D(3)/2;

            % store in array for quick access
            msg = [t0 obj.ID IODC IODE t_oc af0 af1 af2 T_GD ...
                   eph.t_oe eph.a eph.e eph.i eph.RAAN eph.w eph.M0 ...
                   eph.A eph.VP eph.Cx eph.Cy eph.Cz];
        end

        function AP = txlinkbudget(obj,r,beta)
            %TXLINKBUDGET Computes the received power at the user antenna.
            %   Input:
            %    - user; User object instance
            %    - r; transmitter-receiver ranges (m)
            %    - beta; transmitter-user angle (rad)
            %   Output:
            %    - AP; power at the user antenna, dBW
            arguments
                obj     (1,1)   NavSatellite
                r       (1,:)   double {mustBePositive}
                beta    (1,:)   double
            end

            % link budget calculations to obtain C/N0
            freq = obj.ant.freq;
            Ad = 20 * log10((obj.c/freq)./(4*pi*r));    % dB, FSPL
            Ae = 0;                                     % dB, no atmospheric attenuation
            AP = obj.ant.getEIRP(beta) + Ad + Ae;       % dBW, gain before receiver
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
                tol     (1,1)   double = 1e-3
            end

            n = length(ts);
            xuser = user.getstates(ts, 'J2000');
            xref = obj.traj(1).get(ts, 'J2000');

            % initial guess for transmit time is receive time
            tt = ts;
            % instantaneous range at measurement times, in m
            r = sqrt(sum((xref(1:3,:) - xuser(1:3,:)).^2, 1));
            % range-rate of measurements, in m/s
            dr = zeros(1,n);

            for i=1:n
                rlast = r(i);

                for j=1:20
                    dt = rlast / obj.c;         % range to time-of-flight (s)
                    tj = ts(i) - dt;            % time offset guess
                    % updated range guess
                    xref(:,i) = obj.traj(1).get(tj, 'J2000');
                    rj = norm(xref(1:3,i) - xuser(1:3,i));

                    % if iteration is converging
                    if abs(rj - rlast) < tol
                        tt(i) = tj;
                        r(i) = rj;
                        break;
                    elseif j == 10
                        error("timeofflight:notConverged", ...
                            "Time %d failed to converge in %d iterations.", i, j);
                    end

                    rlast = rj;                 % update iteration
                end
                
                % compute range-rate by finding projection of relative
                % velocity along line-of-sight direction
                los = xref(1:3,i) - xuser(1:3,i);
                los = los / norm(los);
                vrel = xref(4:6,i) - xuser(4:6,i);
                dr(i) = vrel' * los;
            end
        end
    end


    methods (Static)
        function tlabel = sandpile3D(t,data,item)
            %SANDPILE3D Plots a sandpile-looking breakdown of data on the
            %current 3D figure. Part of a larger figure, so does not format
            %the current figure at all.
            %   Input:
            %    - t; timestamps of data in seconds
            %    - data; rows are contributors, cols are timestamps t
            %    - item; number of item to plot on y axis
            arguments
                t       (1,:)   double
                data    (:,:)   double {mustBeNonnegative}
                item    (1,1)   {mustBeNonnegative,mustBeInteger}
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
            py = item * ones(1, length(px));        % position
            data = [zeros(1,size(data,2)); data];   % add buffer row
            for i=1:m                               % plot patches
                data(i+1,:) = data(i+1,:) + data(i,:);  % add data
                pz = [data(i,:) flip(data(i+1,:))];     % z data
                color = colors(mod(i-1,c)+1,:);         % get color
                patch(px, py, pz, color, "EdgeColor", "k");
                if i == 1, hold on; end
            end

            % return what the time label should be
            tlabel = "Time " + units;           
        end

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
