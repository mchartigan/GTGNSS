classdef User < handle
    %USER Stores trajectory, receiver, and antenna information for user
    %navigation.
    
    properties
        % user trajectory 
        traj    (1,1)   function_handle = @(varargin) disp([])
        % reference frame for user trajectory
        frame   (1,:)   {mustBeText} = 'J2000'
        % receiver
        rec     (1,1)   Receiver = Receiver("none")
        % receiver antenna
        ant     (1,1)   ReceiveAntenna
    end
    
    methods
        function obj = User(fx,frame,rec,ant)
            %USER Construct a User instance.
            %   Input:
            %    - fx; function handle for user trajectory @(t) in seconds
            %          past J2000, returns [pos (km); vel (km/s)] in frame
            %    - frame; reference frame of fx
            %    - rec; Receiver object instance
            %    - ant; ReceiveAntenna object instance
            
            if (nargin(fx) ~= 1 && nargin(fx) ~= -1) || (nargout(fx) ~= 1 && nargout (fx) ~= -1)
                error("User:invalidInput", "fx must take (t) as an input and output (x)");
            end

            obj.traj = fx;
            obj.frame = frame;
            obj.rec = rec;
            obj.ant = ant;
        end
        
        function x = getstate(obj,ts,frame)
            %GETSTATE Returns state data in the requested frame.
            %   Input:
            %    - ts; time steps, seconds past J2000
            %    - frame; reference frame to provide data in
            
            x = obj.traj(ts);
            if ~strcmp(frame,obj.frame)
                for i=1:length(ts)
                    x(:,i) = cspice_sxform(obj.frame,frame,ts(i)) * x(:,i);
                end
            end
        end
    end
end

