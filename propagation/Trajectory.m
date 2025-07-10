classdef Trajectory < handle
    %TRAJECTORY Generic trajectory class, which is used to describe, store, and
    %manipulate time-evolving data. General operations supported are
    %storing, accessing, and spline interpolation into function handles.

    properties
        t0      (1,1)   double  % start time
        x0      (:,1)   double  % start state
        ts      (1,:)   double  % time series
        xs      (:,:)   double  % states over time
        frame   (1,:)   char    % reference frame of states (if applicable)
        dim     (1,1)   {mustBeInteger,mustBePositive} = 1  % state dimension
    end
    properties (Access = private)
        pp      (1,1)   struct  % piecewise polynomial struct
    end

    methods
        function obj = Trajectory(ts,xs,frame)
            %TRAJECTORY Construct a Trajectory instance.
            %   Input:
            %    - ts; time series of data
            %    - xs; array, rows are states and cols correspond to ts
            %    - frame; optional, SPICE reference frame of states
            arguments
                ts = []
                xs = []
                frame = ''
            end
            
            if nargin ~= 0
                obj.t0 = ts(1);
                obj.x0 = xs(:,1);
                obj.ts = ts;
                obj.xs = xs;
                obj.frame = frame;
                obj.dim = size(xs, 1);
            end
        end

        function x = get(obj,t,outframe)
            %GET Returns the state at t. To get interpolated values,
            %obj.interp() must be called first.
            %   Input:
            %    - t; time(s) to access
            %    - outframe; optional, frame to provide data in (if prev.
            %       specified)
            arguments
                obj         (1,1) Trajectory
                t           (1,:) double
                outframe    (1,:) char = ''
            end

            % check if values are explicitly defined
            [~,ii] = ismember(t, obj.ts);
            if any(ii == 0)     % if not, interpolate
                if isempty(fieldnames(obj.pp))  % error if obj.interp not called yet
                    % error("Trajectory:invalidGet", ...
                    %     "obj.interp must be called before requesting interpolated data.");
                    obj.interp();
                end

                x = ppval(obj.pp, t);   % get values from piecewise polynomial
            else                % else, get values directly
                x = obj.xs(:, ii);
            end

            % check if frame transformations requested and necessary
            if ~isempty(outframe) && ~strcmpi(outframe, obj.frame)
                if isempty(obj.frame)   % error if obj.frame not defined
                    error("Trajectory:invalidFrame", ...
                        "Output frame specified but not set in the object.");
                end

                % check if state dimension is valid (3 or 6)
                n = size(obj.xs, 1);
                if n == 6       % do pos,vel transformation
                    for i=1:length(t)
                        x(:,i) = cspice_sxform(obj.frame, outframe, t(i)) * x(:,i);
                    end
                elseif n == 3   % do pos transformation
                    for i=1:length(t)
                        x(:,i) = cspice_pxform(obj.frame, outframe, t(i)) * x(:,i);
                    end
                else            % can't transform data (wrong size)
                    error("Trajectory:invalidSize", ...
                        "State must be 3D or 6D, rather than %dD.", n);
                end
            end
        end

        function interp(obj)
            %INTERP Generates the spline interpolation of the trajectory at t.
            obj.pp = spline(obj.ts, obj.xs);
        end

        function traj = sub(obj,t0,tf)
            %SUB Returns a Trajectory object that's a subset of the current.
            arguments
                obj (1,1)   Trajectory
                t0  (1,1)   double
                tf  (1,1)   double
            end

            ii = and(obj.ts >= t0, obj.ts <= tf);
            traj = Trajectory(obj.ts(ii), obj.xs(:,ii), obj.frame);
        end
    end
end