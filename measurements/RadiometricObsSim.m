classdef RadiometricObsSim < handle
    %RADIOMETRICOBSSIM Main class for simulating GNSS-like radiometric
    %measurements (GPS, Galileo, LCRNS, etc.).
    %   This class interfaces with NavSatellite and Receiver to create
    %   pseudorange / pseudorange-rate / Doppler measurements, provide
    %   measurement models, handle navigation message parsing, and provide
    %   measurement error estimates. Takes as input a constellation of
    %   NavSatellites and a Receiver design.

    properties
        sats    (1,:)   NavSatellite
    end

    properties (Constant)
        GM = 4.902800e3;                % km^3/s^2, DE440 gravitation parameter of Moon
        O_dot = 2.6618857610e-6;        % rad/s, rotation rate of Moon
        R = 1738100;                    % meters, equatorial radius
        ee = 0.049142364109;            % eccentricity of lunar ellipsoid
    end

    methods
        function obj = RadiometricObsSim(sats)
            %RADIOMETRICOBSSIM Construct a RadiometricObsSim instance.
            %   Detailed explanation goes here
            obj.sats = sats;
        end

        function [y,x_s] = computemeas(obj,t,x,user)
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
            state2eph = @(x_) x_(1:6,:);
            modeleph = @(t_,~) state2eph(obj.getmodelstates(t_));
            tt = obj.timeofflight(t, modeleph, user);
            x_s = obj.getmodelstates(tt, user.frame);
            r_s = x_s(1:3);
            v_s = x_s(4:6);

            dr = (r_s - r_u)*1e3;   % relative user->sat position (m)
            dv = (v_s - v_u)*1e6;   % relative user->sat velocity (mm/s)
            rho = norm(dr);         % scalar range
            dvdr = dv'*dr;          % dot product of dv and dr
            y = [rho      + (x(7) - x_s(7))
                 dvdr/rho + (x(8) - x_s(8))*1e3];
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

            dr = (r_s - r_u)*1e3;   % relative user->sat position (m)
            dv = (v_s - v_u)*1e6;   % relative user->sat velocity (mm/s)
            rho = norm(dr);         % scalar range
            dvdr = dv'*dr;          % dot product of dv and dr
            H = [-dr'/rho*1e3                     0 0 0        1 0   0
                  (dr'*dvdr/rho^3 - dv'/rho)*1e3 -dr'/rho*1e6  0 1e3 0];
        end
    end

    methods (Static)
        function x = geteph(t,ID,msg)
            %GETEPH Read the navigation message and return the satellite
            %ephemeris at time t.
            %   Input:
            %    - t   (1,1) double; time, seconds past J2000
            %    - ID  (1,1) double; satellite ID
            %    - msg (:,:) double; navigation message data, in array. See
            %       NavSatellite.generatenavmsg() for format.

            tau = t - msg(10);                  % time since ephemeris epoch
            RAAN = msg(14) + msg(18)*tau;       % rad, right ascension
            i  = msg(13) + msg(17)*tau;         % rad, inclination
            w = msg(15) + msg(19)*tau;          % rad, arg. of perilune
            n = sqrt(RadiometricObsSim.GM/msg(11)^3);   % rad/s, mean motion
            M = msg(16) + n*tau;                % rad, mean anomaly
            [f, E] = mean2true(M, msg(12));     % rad, true anomaly
            
            u = f + w;                          % rad, argument of latitude
            r = msg(11) * (1 - msg(12)*cos(E)); % km, radius
            % get position in perifocal frame
            xp = r*cos(u);
            yp = r*sin(u);
            % get position in body-fixed frame
            xk = xp*cos(RAAN) - yp*cos(i)*sin(RAAN);
            yk = xp*sin(RAAN) + yp*cos(i)*cos(RAAN);
            zk = yp*sin(i);
        
            % get SV velocity parameters
            Edot = n / (1 - msg(12)*cos(E));
            fdot = Edot*sqrt(1-msg(12)^2)/(1-msg(12)*cos(E));
            idot = msg(17);
            udot = fdot + msg(19);
            RAANdot = msg(18);
            rdot = msg(12)*msg(11)*Edot*sin(E);
            xpdot = rdot*cos(u) - r*udot*sin(u);
            ypdot = rdot*sin(u) + r*udot*cos(u);
            xkdot = -xp*RAANdot*sin(RAAN) + xpdot*cos(RAAN) - ypdot*sin(RAAN)*cos(i) - yp*(RAANdot*cos(RAAN)*cos(i) - idot*sin(RAAN)*sin(i));
            ykdot =  xp*RAANdot*cos(RAAN) + xpdot*sin(RAAN) + ypdot*cos(RAAN)*cos(i) - yp*(RAANdot*sin(RAAN)*cos(i) + idot*cos(RAAN)*sin(i));
            zkdot = ypdot*sin(i) + yp*idot*cos(i);
        
            x = [xk yk zk xkdot ykdot zkdot]';
        end
    end
end