classdef IMU < handle
    %IMU Model of an inertial measurement unit. Used to provide body
    %acceleration info to propagators and/or filters.

    properties
        % Trajectory instance of s/c acceleration over time
        a_m     (1,1)   Trajectory = Trajectory([0 inf],zeros(3,2),'J2000')
        % Velocity random walk specral density, m^2/s^3
        Sn      (1,1)   double {mustBeNonnegative} = 0
        % sample step size, s
        step    (1,1)   double = 1
    end

    methods
        function obj = IMU(accel,N)
            %IMU Construct an IMU instance.
            %   Input:
            %    - accel; Trajectory instance that describes the true body
            %       acceleration of the IMU
            %    - N; velocity random walk, in units of m/s^(3/2)

            if nargin ~= 0
                obj.a_m = accel;
                obj.Sn = N^2;
            end
        end

        function a = read(obj,t)
            %READIMU Returns sensed acceleration (in ICRF frame) at time t.
            %   Input:
            %    - t; time in s past J2000

            a = obj.a_m.get(t,'J2000') + mvnrnd([0 0 0], eye(3)*obj.Sn*obj.step)';
        end

        function Q = noise(obj,dt)
            %NOISE Returns the discrete-time process noise (pos,vel) over the
            %interval.
            Qvel = eye(3)*obj.Sn*obj.step;
            dt = dt/obj.step;
            Q = [Qvel * dt^3/3 Qvel * dt^2/2; Qvel * dt^2/2 Qvel * dt];
        end
    end
end