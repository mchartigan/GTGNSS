classdef NavSatellite < handle
    %NavSatellite Class for describing the properties and trajectory of a
    %satellite that is part of a radionavigation satellite system (whether
    %that's earth-based GNSS or lunar).
    
    properties
        % [satellite reference trajectory; clock reference trajectory]
        traj    (2,1)   Trajectory
        % propagator instance (default one so MATLAB doesn't throw a fit)
        prop    (1,1)   SatellitePropagator
        % nav filter
        filter
        % antenna object
        ant     (1,1)   TransmitAntenna
        % nav message update cadence info
        cadence (1,1)   double {mustBeNonnegative}
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
        function obj = NavSatellite(traj,prop,filter,meas,ant,ID,cadence,debug)
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

            % permit empty instantiation
            if nargin ~=0
                % assign instances
                obj.traj = traj;
                obj.prop = prop;
                obj.ant = ant;
                obj.cadence = cadence;
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
    
                obj.ID = ID;
                if ~debug, warning('off', 'NavSatellite:debug'); end
                obj.DEBUG = debug;
            end
        end

        function [T,dT,AP,msg,err,var] = transmitsignal(obj,ts,user)
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
            %    - AP; power at the user antenna, dBW
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
            T  = r / obj.c  + err.total(1,:);
            dT = dr / obj.c + err.total(2,:);

            % ANTENNA MASK %
            AP = obj.txlinkbudget(r);
            % add transmitter mask to link budget %
            xuser = user.getstates(ts, 'J2000');
            xsat  = obj.traj(1).get(tt, 'J2000');
            % get nadir direction at user at each time step
            nadir = xuser(1:3,:) ./ sqrt(sum(xuser(1:3,:).^2, 1));
            % get zenith direction at sat at each time step
            zenith = -xsat(1:3,:) ./ sqrt(sum(xsat(1:3,:).^2, 1));
            % compute angle between nadir and satellite
            tosat = acos(sum(nadir .* los, 1));
            % compute angle between zenith and user
            touser = acos(sum(zenith .* -los, 1));
            % -300 dB if angles are over off-boresight mask angle
            AP = AP - 300 * (tosat > pi/2 - user.ant.mask);
            AP = AP - 300 * (touser > pi/2 - obj.ant.mask);

            % PLANET INTERSECTION CALCULATION %
            R = cspice_bodvrd('MOON', 'RADII', 3);
            R = R(1);           % radius of moon
            for i=1:length(ts)
                x_s = xsat(1:3,i);
                r_s = norm(x_s);
                r_su = norm(xuser(1:3,i) - x_s);
                u_us = los(1:3,i);
                a_sm = asin(R/r_s);
                a_su = acos(u_us' * x_s / r_s);
                r_t = sqrt(r_s^2 - R^2);
                % if the moon center/moon tangent angle from the satellite
                % POV is bigger than the moon center/sat-user angle and the
                % range is > moon tangent range, satellite is out of view.
                if a_su < a_sm && r_su > r_t
                    T(i)  = NaN;
                    dT(i) = NaN;
                end
            end
        end

        function [err,var,msg,los] = getSISE(obj,tt,ts,user)
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
            % get reference trajectory of satellite and clock
            xref = zeros(9,n);
            xref(1:6,:) = obj.traj(1).get(tt, 'J2000');
            xref(7:9,:) = obj.traj(2).get(tt);
            % get reference trajectory of user
            xuser = user.getstates(ts, 'J2000');
            
            % create line-of-sight direction
            los = xref(1:3,:,1) - xuser(1:3,:,1);
            los = los ./ sqrt(sum(los.^2, 1));

            % generate nav update times
            tmsg = [tt(1)-1:obj.cadence:tt(end) tt(end)];
            [xmsg, Pmsg] = obj.getnavstates(tmsg);
            msg = zeros(length(tmsg)-1, 19);
            Pnav = zeros(9,9,n);        % store starting uncertainty
            xprop = zeros(9,n);         % store the propagated states
            Pprop = zeros(9,9,n);

            % iterate over update times
            for i=1:length(tmsg)-1
                % nav states are applicable starting at tt(1)-1, so this logic
                % should cover all tt
                jj = and(tt > tmsg(i), tt <= tmsg(i+1));
                % store nav uncertainty
                Pnav(:,:,jj) = repmat(Pmsg(:,:,i), 1, 1, length(jj));
                % propagate states over given times. tmsg(i) provided so
                % trajectory starts at appropriate time
                tsub = [tmsg(i) tt(jj)];
                xsub = obj.prop.runat(tsub, xmsg(:,i), 'J2000');
                % cut out xmsg(:,i) since it may not align with tt
                xprop(:,jj) = xsub(:,2:end);
                Pprop(:,:,jj) = ...
                    obj.prop.proplyapunov(tt(jj), xmsg(:,i), Pmsg(:,:,i));
                % provide these propagated states (plus initial one) as a
                % Trajectory and create a navigation message about it.
                subeph = Trajectory(tsub, xsub(1:6,:), 'J2000');
                subclk = Trajectory(tsub, xsub(7:9,:));
                msg(i,:) = obj.generatenavmsg([subeph; subclk]);
            end
            
            % compute errors and variances
            err_prop = xprop - xref;
            Pprop = Pprop - Pnav;           % separate errors from initial est. and propagation
            err.eph_prop = zeros(2,n);      % error due to initial OD and propagation
            err.clk_prop = zeros(2,n);      % error due to clock est. and propagation
            var.eph_est  = zeros(2,n);      % variance of initial OD
            var.eph_prop = zeros(2,n);      % variance of state propagation (-OD)
            var.clk_est  = zeros(2,n);      % variance of initial clock est.
            var.clk_prop = zeros(2,n);      % variance of clock propagation (-est.)
            % project onto line-of-sight direction
            for i=1:n
                % range and range-rate error (in s and s/s) due to propagation
                err.eph_prop(1,i) = err_prop(1:3,i)' * los(:,i) / obj.c_km;
                err.eph_prop(2,i) = err_prop(4:6,i)' * los(:,i) / obj.c_km;
                % " due to onboard clock offset from proper time
                err.clk_prop(:,i) = err_prop(7:8,i) / obj.c;
                % variance from OD (in s^2 and s^2/s^2)
                var.eph_est(1,i) = los(:,i)' * Pnav(1:3,1:3,i) * los(:,i) / obj.c_km2;
                var.eph_est(2,i) = los(:,i)' * Pnav(4:6,4:6,i) * los(:,i) / obj.c_km2;
                var.clk_est(:,i) = diag(Pnav(7:8,7:8,i)) / obj.c2; 
                % " from state propagation
                var.eph_prop(1,i) = los(:,i)' * Pprop(1:3,1:3,i) * los(:,i) / obj.c_km2;
                var.eph_prop(2,i) = los(:,i)' * Pprop(4:6,4:6,i) * los(:,i) / obj.c_km2;
                var.clk_prop(:,i) = diag(Pprop(7:8,7:8,i)) / obj.c2;
            end

            % total everything up
            err.total = err.eph_prop + err.clk_prop;
            var.total = var.eph_est + var.eph_prop + var.clk_est + var.clk_prop;
        end

        function [erec,vrec,mask] = getUEE(obj,ts,r,dr,user,los,mask)
            %GETUEE Returns the user equipment error for a given set of
            %parameters.
            %   Input:
            %    - ts; measurement times in seconds past J2000
            %    - r; transmitter-receiver ranges (m) to compute error for
            %    - dr; transmitter-receiver range-rates (mm/s)
            %    - user; User object instance
            %    - los; line-of-sight direction from user to satellite
            %    - mask; pre-existing mask for range and -rate measurements
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                r       (1,:)   double {mustBePositive}
                dr      (1,:)   double
                user    (1,1)   User
                los     (3,:)   double
                mask    (2,:)   double = []
            end

            % compute link budget and receiver noise
            CN0 = obj.txlinkbudget(user, r);
            % add mask to link budget %
            % get nadir direction at user at each time step
            nadir = user.xs(1:3,:) ./ sqrt(sum(user.xs(1:3,:).^2, 1));
            % get zenith direction at sat at each time step
            zenith = -obj.xr(1:3,:) ./ sqrt(sum(obj.xr(1:3,:).^2, 1));
            % compute angle between nadir and satellite
            tosat = acos(sum(nadir .* los, 1));
            touser = acos(sum(zenith .* -los, 1));
            % -300 dB/Hz if angle is over off-boresight mask angle
            CN0 = CN0 - 300 * (tosat > pi/2 - user.ant.mask);
            CN0 = CN0 - 300 * (touser > pi/2 - obj.ant.mask);
            % if 1
            %     RP = CN0(CN0 > -200);
            %     figure();
            %     plot(RP);
            %     fprintf("max: %f\nmin: %f\n", max(RP), min(RP));
            % end
            [erec, vrec, mask1] = user.rec.noise(CN0, ts, r, dr, obj.clock);

            % build and assign mask
            if ~isempty(mask)
                mask = and(mask1, mask);
            else
                mask = mask1;
            end

            % plot CN0 for simulation
            if obj.DEBUG
                figure();
                plotformat("APA", 0.4, "scaling", 1, "coloring", "science");
                elev = (pi/2 - tosat) * 180/pi;
                plot(elev(mask(1,:)), CN0(mask(1,:)), "LineWidth", 2);
                grid on;
                xlabel("Elevation angle (\circ)");
                ylabel("C/N0 (dB-Hz)");
                title("C/N0 vs. User Antenna Elevation Angle");

                figure();
                plotformat("APA", 0.4, "scaling", 1, "coloring", "science");
                plot(tosat(mask(1,:)) * 180/pi, touser(mask(1,:)) * 180/pi, "LineWidth", 2);
                grid on;
                xlabel("User off-boresight angle (\circ)");
                ylabel("Satellite off-boresight angle (\circ)");
            end
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

            xn = zeros(9, length(ts));
            if isa(obj.filter, 'struct')
                xn(1:6,:) = obj.traj(1).get(ts, 'J2000');
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
                xn(1:6,:) = obj.traj(1).get(ts, 'J2000');
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
            eph = obj.prop.orbit.AFSfit(traj_eph, 10);
            IODE = mod(floor(eph.t_oe), 1024);      % simple hash

            % update clock offsets
            t_oc = t0;
            IODC = mod(floor(t_oc), 1024);          % simple hash
            [~,D] = obj.prop.clock.modelfit(traj_clk, t0);
            af0 = D(1);
            af1 = D(2);
            af2 = D(3);

            % store in array for quick access
            msg = [t0 obj.ID IODC IODE t_oc af0 af1 af2 T_GD ...
                   eph.t_oe eph.a eph.e eph.i0 eph.RAAN0 eph.w0 eph.M0 ...
                   eph.idot eph.RAANdot eph.wdot];
        end

        function AP = txlinkbudget(obj,r)
            %TXLINKBUDGET Computes the received power at the user antenna.
            %   Input:
            %    - user; User object instance
            %    - r; transmitter-receiver ranges (m)
            %   Output:
            %    - AP; power at the user antenna, dBW
            arguments
                obj     (1,1)   NavSatellite
                r       (1,:)   double {mustBePositive}
            end

            % link budget calculations to obtain C/N0
            freq = obj.ant.freq;
            Ad = 20 * log10((obj.c/freq)./(4*pi*r));    % dB, FSPL
            Ae = 0;                                     % dB, no atmospheric attenuation
            AP = obj.ant.P + obj.ant.gain + Ad + Ae;    % dBW, gain before receiver
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
            xuser = user.getstates(ts, 'J2000');
            xref0 = obj.traj(1).get(ts(1), 'J2000');

            % initial guess for transmit time is receive time
            tt = ts;
            % instantaneous range at measurement times, in m
            r = sqrt(sum((xref0(1:3,:) - xuser(1:3,:)).^2, 1)) * 1e3;
            % final s/c states at tt
            xref = zeros(6,n);
            % range-rate of measurements, in m/s
            dr = zeros(1,n);

            for i=1:n
                rlast = r(i);

                for j=1:100
                    dt = rlast / obj.c;         % range to time-of-flight (s)
                    tj = ts(i) - dt;            % time offset guess
                    % updated range guess
                    xref(:,i) = obj.traj(1).get(tj, 'J2000');
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
                dr(i) = vrel' * los * 1e3;      % convert km/s -> m/s
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
