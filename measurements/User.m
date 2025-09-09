classdef User < handle
    %USER Stores trajectory, receiver, and antenna information for user
    %navigation.
    
    properties
        % user trajectory 
        motion  (1,1)   Trajectory
        clock   (1,1)   Trajectory
        % receiver
        rx      (1,1)   Receiver
        % receiver antenna
        ant     (1,1)   ReceiveAntenna
    end
    
    methods
        function obj = User(motion,clock,rx,ant)
            %USER Construct a User instance.
            %   Input:
            %    - motion; user motion Trajectory instance
            %    - clock; clock offset Trajectory instance (same times as
            %       motion)
            %    - rec; Receiver object instance
            %    - ant; ReceiveAntenna object instance

            if nargin ~= 0
                obj.motion = motion;
                obj.clock = clock;
                obj.rx = rx;
                obj.ant = ant;
            end
        end
        
        function x = getstates(obj,ts,frame)
            %GETSTATE Returns state data in the requested frame.
            %   Input:
            %    - ts; time steps, seconds past J2000
            %    - frame; reference frame to provide data in
            
            x = zeros(obj.motion.dim + obj.clock.dim, length(ts));
            % Populate state data based on motion and clock trajectories
            x(1:obj.motion.dim,:) = obj.motion.get(ts, frame);
            x(obj.motion.dim+1:end,:) = obj.clock.get(ts);
        end
    end
end

