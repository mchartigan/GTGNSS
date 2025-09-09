 classdef RadiometricObsSim < Measurement
    %RADIOMETRICOBSSIM Main class for simulating GNSS-like radiometric
    %measurements (GPS, Galileo, LCRNS, etc.).
    %   This class interfaces with NavSatellite and Receiver to create
    %   pseudorange / pseudorange-rate / Doppler measurements, provide
    %   measurement models, handle navigation message parsing, and provide
    %   measurement error estimates. Takes as input a constellation of
    %   NavSatellites and a Receiver design.

    properties
        sats    (1,:)   NavSatellite
        % No. of satellites in constellation
        nsats   (1,1)   {mustBeInteger,mustBePositive} = 1
        user    (1,1)   User
        % No. of measurements user can make (DLL,PLL,FLL)
        m       (1,1)   {mustBeInteger,mustBePositive} = 1
        dim             = 1
        % Navigation message frames
        msgs    (:,:)   double
        % reference frame of measurements (for geteph())
        frame   (1,:)   {mustBeText} = 'MOON_ME'
    end

    properties (Constant)
        GM = 4.902800e3;                % km^3/s^2, DE440 gravitation parameter of Moon
        O_dot = 2.6618857610e-6;        % rad/s, rotation rate of Moon
        R = 1738100;                    % meters, equatorial radius
        ee = 0.049142364109;            % eccentricity of lunar ellipsoid
        c = 299792458                   % m/s, speed of light
    end

    methods
        function obj = RadiometricObsSim(sats,user)
            %RADIOMETRICOBSSIM Construct a RadiometricObsSim instance.
            %   Input:
            %    - sats; row vector of NavSatellite instances comprising
            %       constellation
            %    - user; User instance detailing their receiver,
            %       trajectory, etc.

            obj.sats = sats;
            obj.nsats = length(sats);
            obj.user = user;
            obj.m = 1 + user.rx.PLL + user.rx.FLL;
            obj.dim = obj.nsats * obj.m;
        end

        function [y,R,var] = getmeas(obj,ts)
            %GETMEAS Returns measurements, that the receiver can make, for
            %the user at times ts.
            %   Input:
            %    - ts; time steps ts of measurements (s past J2000)
            %    - user; User instance describing their state
            arguments
                obj     (1,1)   RadiometricObsSim
                ts      (1,:)   double
            end

            n = length(ts);
            y = nan(3*obj.nsats,n);
            x_user = obj.user.getstates(ts, obj.frame);
            obj.msgs = [];
            var = cell(obj.nsats,1);
            R = zeros(3*obj.nsats,3*obj.nsats,n);

            for i=1:obj.nsats
                % FIELD INCOMING MEASUREMENTS -- y_raw(ts) %
                [T,dT,AP,msg,~,var1] = obj.sats(i).transmitsignal(ts, obj.user);
                [y_raw,~,var2] = obj.user.rx.tracksat(ts, T, dT, AP);
                obj.msgs = [obj.msgs; msg];

                % MERGE VARIANCES %
                % this function hardcodes delta-pseudorange measurements as
                % a difference of the previous and current ones.
                % TODO: variable measurement spacing?

                % update total based on new measurement model
                % delta-pseudorange measurements
                var2.total(2,:) = [nan var1.total(2,2:end).*(ts(2:end) - ts(1:end-1)).^2 + ...
                    var2.total(2,1:end-1) + var2.total(2,2:end)];
                % pseudorange measurements
                var2.total(1,:) = var2.total(1,:) + var1.total(1,:);
                % Doppler measurements
                var2.total(3,:) = var2.total(3,:) + var1.total(2,:);
                % expand SISE variances to DLL, PLL, FLL measurements
                fields1 = fieldnames(var1);
                for j=1:numel(fields1)
                    var1.(fields1{j}) = [var1.(fields1{j})(1,:); ...
                        var1.(fields1{j})(1,:); var1.(fields1{j})(2,:)];
                end
                % merge fields over from var2
                fields2 = fieldnames(var2);
                for j=1:numel(fields2)
                    var1.(fields2{j}) = var2.(fields2{j});
                end
                % assign to cell array
                var{i} = var1;


                % BUILD PSEUDORANGE AND DOPPLER %
                % All the while adding the user clock bias and drift
                % m, pseudorange
                y(i,:) = y_raw(1,:) + x_user(7,:);
                % Measurement underweighting?
                R(i,i,:) = var1.total(1,:);
                % m, delta-pseudorange
                if obj.user.rx.PLL
                    y(i+obj.nsats,:) = y_raw(2,:) + x_user(7,:);
                    % change measurements to delta-pseudorange
                    y(i+obj.nsats,:) = [nan ...
                        y(i+obj.nsats,2:end) - y(i+obj.nsats,1:end-1)];
                    R(i+obj.nsats,i+obj.nsats,:) = var1.total(2,:);
                end
                % m/s, Doppler
                if obj.user.rx.FLL
                    y(i+2*obj.nsats,:) = y_raw(3,:) + x_user(8,:);
                    R(i+2*obj.nsats,i+2*obj.nsats,:) = var1.total(3,:);
                end
            end
        end

        function [y,xs] = computemeas(obj,tr,x,tprev,xprev)
            %COMPUTEMEAS Computes the ideal pseudorange and Doppler
            %measurements at time t between the user and NavSatellites
            %   Input:
            %    - tr; measurement time (or user's best estimate thereof),
            %       seconds past J2000
            %    - x; (best est. of) state of USER at time t in MOON_ME
            %       [pos (km); vel (km/s); t bias (s); drift (s/s); rate (s/s^2)]
            arguments
                obj     (1,1)   RadiometricObsSim
                tr      (1,1)   double
                x       (9,:)   double
                tprev   (1,1)   double = NaN
                xprev   (9,:)   double = zeros(9,1)
            end

            r_u = x(1:3);           % user position
            v_u = x(4:6);           % user velocity

            y = nan(3*obj.nsats,1);
            xs = zeros(9,obj.nsats);

            for i=1:obj.nsats
                % find transmission time w.r.t meas, based on nav msg knowledge
                % of satellite states
                tt = obj.timeofflight(tr,x,obj.sats(i).ID,obj.msgs);
                [x_s,T] = obj.geteph(tt,obj.sats(i).ID,obj.msgs);
                % rotate to inertial if that's how we're managing things
                if strcmp(obj.frame, 'J2000')
                    x_s(1:6) = cspice_invstm(T) * x_s(1:6);
                elseif ~strcmp(obj.frame, 'MOON_ME')
                    error("computemeas:invalidFrame", ...
                        "Cannot recognize frame %s. Supported frames are 'J2000', 'MOON_ME'.", obj.frame);
                end
                xs(:,i) = x_s;
                r_s = x_s(1:3); 
                v_s = x_s(4:6);
    
                dr = (r_s - r_u);       % relative user->sat position (m)
                dv = (v_s - v_u);       % relative user->sat velocity (m/s)
                rho = norm(dr);         % scalar range
                dvdr = dv'*dr;          % dot product of dv and dr

                % m, pseudorange (DLL)
                y(i) = rho + x(7) - x_s(7);
                % m/s, Doppler
                if obj.user.rx.FLL
                    y(i + 2*obj.nsats) = dvdr/rho + x(8) - x_s(8);
                end
                % m, pseudorange (PLL)
                if obj.user.rx.PLL && ~isnan(tprev)
                    y(i + obj.nsats) = rho + x(7) - x_s(7);

                    % compute previous pseudorange and difference them
                    % find transmission time w.r.t meas, based on nav msg knowledge
                    % of satellite states
                    tt = obj.timeofflight(tprev,xprev,obj.sats(i).ID,obj.msgs);
                    [x_s,T] = obj.geteph(tt,obj.sats(i).ID,obj.msgs);
                    % rotate to inertial if that's how we're managing things
                    if strcmp(obj.frame, 'J2000')
                        x_s(1:6) = cspice_invstm(T) * x_s(1:6);
                    elseif ~strcmp(obj.frame, 'MOON_ME')
                        error("computemeas:invalidFrame", ...
                            "Cannot recognize frame %s. Supported frames are 'J2000', 'MOON_ME'.", obj.frame);
                    end
                    dr = x_s(1:3) - xprev(1:3);     % relative user->sat position (m)
                    rho = norm(dr);                 % scalar range
                    % change it to delta-pseudorange
                    y(i+obj.nsats) = y(i+obj.nsats) - (rho + xprev(7) - x_s(7));
                end
            end
        end

        function [H,J] = measpartials(obj,tr,x,tprev,xprev)
            %MEASPARTIALS Computes the partial derivative of the measurement
            %model w.r.t. x at the current state.
            %   Input:
            %    - tr; measurement time (or user's best estimate thereof),
            %       seconds past J2000
            %    - x; (best est. of) state of USER at time t in MOON_ME
            %       [pos (km); vel (km/s); t bias (s); drift (s/s); rate (s/s^2)]
            arguments
                obj     (1,1)   RadiometricObsSim
                tr      (1,:)   double
                x       (9,:)   double
                tprev   (1,1)   double = NaN
                xprev   (9,:)   double = zeros(9,1)
            end

            r_u = x(1:3);           % user position
            v_u = x(4:6);           % user velocity

            H = zeros(3*obj.nsats,9);
            J = zeros(3*obj.nsats,9);

            for i=1:obj.nsats
                % find transmission time w.r.t meas, based on nav msg knowledge
                % of satellite states
                tt = obj.timeofflight(tr,x,obj.sats(i).ID,obj.msgs);
                [x_s,T] = obj.geteph(tt,obj.sats(i).ID,obj.msgs);
                % rotate to inertial if that's how we're managing things
                if strcmp(obj.frame, 'J2000')
                    x_s(1:6) = cspice_invstm(T) * x_s(1:6);
                elseif ~strcmp(obj.frame, 'MOON_ME')
                    error("measpartials:invalidFrame", ...
                        "Cannot recognize frame %s. Supported frames are 'J2000', 'MOON_ME'.", obj.frame);
                end
                r_s = x_s(1:3); 
                v_s = x_s(4:6);
    
                dr = (r_s - r_u);       % relative user->sat position (m)
                dv = (v_s - v_u);       % relative user->sat velocity (m/s)
                rho = norm(dr);         % scalar range
                dvdr = dv'*dr;          % dot product of dv and dr

                % m, pseudorange (DLL)
                H(i,:) = [-dr'/rho 0 0 0 1 0 0];
                % m/s, Doppler
                if obj.user.rx.FLL
                    H(i + 2*obj.nsats,:) = ...
                        [(dr'*dvdr/rho^3 - dv'/rho) -dr'/rho 0 1 0];
                end
                % m, pseudorange (PLL)
                if obj.user.rx.PLL && ~isnan(tprev)
                    H(i + obj.nsats,:) = [-dr'/rho 0 0 0 1 0 0];

                    % find transmission time w.r.t meas, based on nav msg knowledge
                    % of satellite states
                    tt = obj.timeofflight(tprev,xprev,obj.sats(i).ID,obj.msgs);
                    [x_s,T] = obj.geteph(tt,obj.sats(i).ID,obj.msgs);
                    % rotate to inertial if that's how we're managing things
                    if strcmp(obj.frame, 'J2000')
                        x_s(1:6) = cspice_invstm(T) * x_s(1:6);
                    elseif ~strcmp(obj.frame, 'MOON_ME')
                        error("computemeas:invalidFrame", ...
                            "Cannot recognize frame %s. Supported frames are 'J2000', 'MOON_ME'.", obj.frame);
                    end
                    dr = x_s(1:3) - xprev(1:3);     % relative user->sat position (m)
                    rho = norm(dr);                 % scalar range

                    % find second measurement matrix
                    J(i + obj.nsats,:) = -[-dr'/rho 0 0 0 1 0 0];
                end
            end
        end

        function tt = timeofflight(obj,tr,user,ID,msg,tol)
            %TIMEOFFLIGHT Based on the provided receive times, find the
            %satellite transmit time and compute the range / range-rate.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - user; (best est. of) user state at time  in MOON_ME
            %    - ID; ID # of satellite to find
            %    - msg; matrix of navigation message data
            %    - tol; iteration tolerance for solving transmission time
            arguments
                obj     (1,1)   RadiometricObsSim
                tr      (1,1)   double
                user    (9,1)   double
                ID      (1,1)   {mustBeInteger,mustBeNonnegative}
                msg     (:,:)   double
                tol     (1,1)   double = 1e-10
            end

            % initial guess at signal time-of-flight
            [xSV,T] = obj.geteph(tr, ID, msg);
            % rotate to inertial if that's how we're managing things
            if strcmp(obj.frame, 'J2000')
                xSV(1:6) = cspice_invstm(T) * xSV(1:6);
            end
            tt = tr;
            dt_old = norm(xSV(1:3) - user(1:3)) / obj.c;

            for i=1:10
                tt = tr - dt_old;           % time offset guess
                % updated time-of-flight guess
                [xSV,T] = obj.geteph(tt, ID, msg);
                % rotate to inertial if that's how we're managing things
                if strcmp(obj.frame, 'J2000')
                    xSV(1:6) = cspice_invstm(T) * xSV(1:6);
                end
                dt = norm(xSV(1:3) - user(1:3)) / obj.c;

                % if iteration is converging
                if abs(dt - dt_old) < tol
                    tt = tr - dt;
                    break;
                end
                dt_old = dt;                % update iteration
            end
            if i == 10
                warning("timeofflight:notConverged", ...
                    "Failed to converge in %d iterations.", i);
            end
        end

        function ploterror(obj,t,yobs,ycomp,var)
            %PLOTERROR Plots psuedorange, delta-pseudorange, and/or Doppler
            %error over time -- along with the corresponding variance.
            %   Measurement model is
            %      yobs = ycomp + error(ts)
            %
            %   Input:
            %    - t; measurement times (s past J2000)
            %    - yobs; observed measurements corresponding to ts
            %    - ycomp; computed measurements corresponding to ts
            %    - var; variance struct output by RadiometricObsSim.getmeas()
            arguments
                obj     (1,1)   RadiometricObsSim
                t       (1,:)   double
                yobs    (:,:)   double
                ycomp   (:,:)   double
                var     (1,1)   struct
            end

            % get plotting timescale and units
            t = t - t(1);
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

            figure();
            plotformat("APA", 0.25*obj.m + 0.25, "coloring", "greyscale");
            colors = colororder;
            tiledlayout(obj.m,1);

            % plot pseudorange error
            nexttile();
            mask = ~isnan(yobs(1,:));
            dt = t(mask);
            p_err = yobs(1,mask) - ycomp(1,mask);
            p_std = 3 * sqrt(var.total(1,mask));
            plot(dt, p_err);
            hold on;
            patch([dt flip(dt)], [-p_std flip(p_std)], colors(2,:), ...
                "FaceAlpha", 0.6, "EdgeColor", "none");
            hold off; grid on;
            axis([t(1) t(end) -inf inf]);
            ylabel("Error (m)");
            title("Pseudorange measurement error");

            % compute delta-pseudorange at measurement spacing and plot
            if obj.user.rx.PLL
                nexttile();
                mask = ~isnan(yobs(2,:));
                dt = t(mask);
                dp_err = yobs(2,mask) - ycomp(2,mask);
                dp_std = 3 * sqrt(var.total(2,mask));
                msgbound = abs(dp_err) > 0.1;
                dp_err(msgbound) = [];
                dp_std(msgbound) = [];
                dt(msgbound) = [];
                plot(dt, dp_err);
                hold on;
                patch([dt flip(dt)], [-dp_std flip(dp_std)], colors(2,:), ...
                    "FaceAlpha", 0.6, "EdgeColor", "none");
                hold off; grid on;
                axis([t(1) t(end) -inf inf]);
                ylabel("Error (m)");
                title("Delta-pseudorange measurement error");
            end

            % plot Doppler error
            if obj.user.rx.FLL
                nexttile();
                mask = ~isnan(yobs(3,:));
                dt = t(mask);
                f_err = yobs(3,mask) - ycomp(3,mask);
                f_std = 3 * sqrt(var.total(3,mask));
                plot(dt, f_err);
                hold on;
                patch([dt flip(dt)], [-f_std flip(f_std)], colors(2,:), ...
                    "FaceAlpha", 0.6, "EdgeColor", "none");
                hold off; grid on;
                axis([t(1) t(end) -inf inf]);
                ylabel("Error (m/s)");
                title("Doppler measurement error");
            end

            % ending data
            xlabel(sprintf("Time %s", units));
        end
    end

    methods (Static)
        function [x,T] = geteph(t,ID,msg)
            %GETEPH Read the navigation message and return the satellite
            %ephemeris at time t.
            %   Input:
            %    - t   (1,1) double; time, seconds past J2000
            %    - ID  (1,1) double; satellite ID
            %    - msg (:,:) double; navigation message data, in array. See
            %       NavSatellite.generatenavmsg() for format.
            %   Output:
            %    - x; satellite state at time t
            %       [pos (km); vel (km/s); t bias (s); drift (s/s); rate (s/s^2)]

            % find line
            for line=size(msg,1):-1:1
                if msg(line,2) == ID && msg(line,1) <= t
                    msg = msg(line,:);
                    break;
                end
                if line == 1
                    error("geteph:invalidMSG", ...
                        "Applicable message could not be found for t=%.0f, ID=%d", ...
                        t, ID);
                end
            end

            tau = t - msg(10);                  % time since ephemeris epoch
            RAAN = msg(14);                     % rad, right ascension
            i  = msg(13);                       % rad, inclination
            w = msg(15);                        % rad, arg. of perilune
            n = sqrt(RadiometricObsSim.GM/msg(11)^3);   % rad/s, mean motion
            M = msg(16) + n*tau;                % rad, mean anomaly
            [f, ~] = mean2true(M, msg(12));     % rad, true anomaly
            
            x = oe2rv(msg(11), msg(12), i, RAAN, w, f, RadiometricObsSim.GM);
            % get rotation matrix from inertial to body-fixed
            eul = msg(17:22)';
            eul(1:3) = eul(1:3) + eul(4:6) * tau;
            T = cspice_eul2xf(eul,3,1,3);
        
            % compute clock offsets
            tau = t - msg(5);
            xc = msg(6) + msg(7)*tau + msg(8)*tau^2;
            yc = msg(7) + 2*msg(8)*tau;
            zc = 2*msg(8);
            x = [x; xc; yc; zc];

            % compute differential orbital corrections
            NC = (length(msg) - 23)/3;  % No. coefficients per axis

            if NC > 0   % if differential corrections provided
                VP = msg(23);                % s, validity period
                Cx = msg(24:23+NC);
                Cy = msg(24+NC:23+2*NC);
                Cz = msg(24+2*NC:end);
                G = [Cx; Cy; Cz];
                tau = 2*(tau)/VP - 1;
                x(1:3) = x(1:3) + G * (tau.^(0:NC-1))';
                x(4:6) = x(4:6) + G * ((0:NC-1).*(tau.^([0 0:NC-2])) * 2/VP)';
            end

            % apply transformation
            x(1:6) = T * x(1:6);
        end
    end
end