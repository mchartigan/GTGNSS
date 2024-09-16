classdef GNSSmeasurements < handle
    %GNSSMEASUREMENTS Handler for determining satellites visible to user,
    %computing measurements, determining measurement partials, etc.
    
    properties
        tmeas   (1,:) double        % measurement times
        fsats   (:,1) cell          % cell array of handles to satellites
        R       (:,:) double        % measurement covariance
        earth   (1,1) struct        % structure containing info about earth (R)
        moon    (1,1) struct        % structure containing info about moon (R)
        visible (:,:) double        % visibility of n sats at m time steps
        ymeas   (:,:) double        % measurement array
        user    (9,:) double        % true state of user
        max_a   (1,1) double        % maximum satellite beam width
        range   (1,1) logical = false                           % boolean, is range desired?
        rate    (1,1) logical = false                           % boolean, is range-rate desired?
        nsats   (1,1) {mustBePositive, mustBeInteger} = 1       % number of GNSS satellites
        m       (1,1) {mustBePositive, mustBeInteger} = 1       % # of steps in tmeas

        c = 299792.458              % km/s, speed of light
    end
    
    methods
        function obj = GNSSmeasurements(obs,tmeas,fsats,R,earth,moon,user,max_a)
            %GNSSMEASUREMENTS Construct an GNSS measurement object
            %instance.
            %   Inputs:
            %    - obs; what measurements to use ("RANGE", "RATE", or "BOTH")
            %    - tmeas; times at which to provide measurements
            %    - fsats; cell array of satellite state function handles;
            %             {@(t) x(t) [pos (km); vel (km/s)]; ... }
            %    - R; measurement covariance matrix, (nsats,nsats) or (2*nsats,2*nsats) depending on obs
            %    - earth; struct containing info about earth (radius at minimum)
            %    - moon; struct containing info about moon (radius at minimum)
            %    - user; true state of user at measurement times [pos (km); vel (km/s); clk]
            %    - max_a; optional, maximum antenna beam-width of GNSS (in
            %             radians)
            
            if nargin < 8
                max_a = 30 * pi / 180;
            elseif nargin < 7
                error("GNSSmeasurements:nargin", ...
                    "Incorrect number of input arguments; accepts 7.")
            end

            % assign measurements desired to booleans
            mustBeOption(obs);
            if strcmpi(obs,"RANGE") || strcmpi(obs,"BOTH")
                obj.range = true;
            end
            if strcmpi(obs,"RATE") || strcmpi(obs,"BOTH")
                obj.rate = true;
            end

            obj.tmeas = tmeas;
            obj.m = length(tmeas);
            obj.fsats = fsats;
            obj.R = R;
            obj.nsats = length(fsats);
            obj.earth = earth;
            obj.moon = moon;
            obj.user = user;
            obj.max_a = max_a;

            % initialize visibility and measurement arrays
            obj.visible = zeros(obj.nsats, obj.m);
            obj.ymeas = zeros((obj.range+obj.rate) * obj.nsats, obj.m);

            % compute visibility and measurements
            obj.computevisible();
        end

        function y = compute(obj,t,x)
            %COMPUTE Computes measurements for user at time t.
            %   Input:
            %    - t; time of measurement
            %    - x; state of user
            arguments
                obj (1,1) GNSSmeasurements
                t   (1,1) double
                x   (9,1) double
            end

            if ~ismember(t,obj.tmeas)
                error("GNSSmeasurements:timestep", ...
                    "t is not in measurement times provided upon initialization.")
            end

            % initialize measurement vector
            y = zeros((obj.range+obj.rate) * obj.nsats,1);
            % get visibility
            visi = obj.visible(:,find(obj.tmeas == t,1));
            % user properties
            r_u = x(1:3);       % user position vector
            v_u = x(4:6);       % user velocity vector

            for i=1:obj.nsats   % consider every satellite
                if visi(i)          % if satellite is visible
                    x_s = obj.fsats{i}(t);          % satellite state
                    r_s = x_s(1:3);                 % satellite position vector
                    v_s = x_s(4:6);                 % satellite velocity vector
                    r_us = r_s - r_u;               % user->sat position vector

                    dr = (v_s - v_u)' * r_us / norm(r_us);  % range-rate

                    if obj.range        % obs == "RANGE"
                        y(i) = norm(r_us) + obj.c*x(7);
                        if obj.rate     % obs == "BOTH"
                            y(obj.nsats+i) = dr + obj.c*x(8);
                        end
                    elseif obj.rate     % obs == "RATE"
                        y(i) = dr + obj.c*x(8);
                    end
                end
            end
        end

        function H = partials(obj,t,x)
            %PARTIALS Computes the partial derivative of the measurement
            %model w.r.t. x at the current state.
            %   Inputs:
            %    - t; time of measurement
            %    - x; state of user
            arguments
                obj (1,1) GNSSmeasurements
                t   (1,1) double
                x   (9,1) double
            end

            if ~ismember(t,obj.tmeas)
                error("GNSSmeasurements:timestep", ...
                    "t is not in measurement times provided upon initialization.")
            end

            % get visibility
            visi = obj.visible(:,find(obj.tmeas == t,1));
            % initialize partials matrix
            H = zeros((obj.range+obj.rate) * obj.nsats, length(x));
            % user properties
            r_u = x(1:3);       % user position vector
            v_u = x(4:6);       % user velocity vector

            for i=1:obj.nsats   % consider every satellite
                if visi(i)          % if satellite is visible
                    x_s = obj.fsats{i}(t);          % satellite state
                    r_s = x_s(1:3);                 % satellite position vector
                    v_s = x_s(4:6);                 % satellite velocity vector

                    dr = r_s - r_u;     % relative user->sat position
                    dv = v_s - v_u;     % relative user->sat velocity
                    rho = norm(dr);     % scalar range
                    dvdr = dv'*dr;      % dot product of dv and dr

                    if obj.range        % obs == "RANGE"
                        H(i,:) = [-dr'/rho 0 0 0 obj.c 0 0];
                        if obj.rate     % obs == "BOTH"
                            H(obj.nsats+i,:) = [(dr'*dvdr/rho^3 - dv'/rho) -dr'/rho 0 obj.c 0];
                        end
                    elseif obj.rate     % obs == "RATE"
                        H(i,:) = [(dr'*dvdr/rho^3 - dv'/rho) -dr'/rho 0 obj.c 0];
                    end
                end
            end
        end
    end

    methods (Access = private)
        function computevisible(obj)
            %COMPUTEVISIBLE Called at instantiation; creates array of
            %satellite visibilities at each time step and generates
            %pseudorange and -range-rate measurements while it's at it.

            for i=1:obj.m       % iterate over time steps
                ti = obj.tmeas(i);          % current time
                x_u = obj.user(:,i);        % user state
                r_u = x_u(1:3);             % user position vector
                v_u = x_u(4:6);             % user velocity vector
                x_e = obj.earth.x(ti);      % earth state
                r_e = x_e(1:3) - r_u;       % earth position w.r.t. user

                % figure();
                % scatter3(r_u(1), r_u(2), r_u(3), 'b', 'x');
                % hold on;
                % scatter3(0, 0, 0, 'k', 'filled', 'o');
                % scatter3(x_e(1), x_e(2), x_e(3), 'b', 'filled', 'o');

                for j=1:obj.nsats       % iterate over each satellite
                    try
                        x_s = obj.fsats{j}(ti);         % satellite state
                    catch               % catch SPICE issues (no data at time)
                        % sat at center of earth (will fail checks)
                        x_s = [x_e; 0; 0; 0];
                    end

                    r_s = x_s(1:3);                 % satellite position vector
                    v_s = x_s(4:6);                 % satellite velocity vector
                    r_us = r_s - r_u;               % user->sat position vector

                    % scatter3(x_s(1), x_s(2), x_s(3), 'r', '+');

                    % compute earth-center / earth-tangent angle, earth-center /
                    % GNSS sat angle, and GNSS-to-earth / GNSS-to-user angle
                    ang_e = asin(obj.earth.R / norm(r_e));
                    ang_s = acos(dot(r_us,r_e) / (norm(r_us)*norm(r_e)));
                    ang_a = acos(dot(r_e-r_us,-r_us) / (norm(r_e-r_us)*norm(-r_us)));

                    % compute moon-center / moon-tangent angle, GNSS-to-moon / 
                    % GNSS-to-user angle
                    ang_m = asin(obj.moon.R / norm(-r_s));
                    ang_n = acos(dot(-r_s,-r_us) / (norm(-r_s)*norm(-r_us)));

                    dr = (v_s - v_u)' * r_us / norm(r_us);  % range-rate

                    % satellite is NOT (further than earth AND within its view angle)
                    % AND user is within the GNSS satellite beam width AND user is NOT
                    % (past the moon and within its view angle from GNSS sat)
                    if ~(norm(r_us) > norm(r_e) && ang_s < ang_m) && ang_a < obj.max_a && ...
                       ~(norm(-r_us) > norm(-r_s) && ang_n < ang_m)
                        obj.visible(j,i) = 1;       % spacecraft is visible

                        % go thru skill tree to assign measurements
                        if obj.range
                            obj.ymeas(j,i) = norm(r_us) + obj.c*x_u(7) + mvnrnd(0,obj.R(j,j));
                        
                            if obj.rate
                                obj.ymeas(obj.nsats+j,i) = dr + obj.c*x_u(8) + ...
                                    mvnrnd(0,obj.R(obj.nsats+j,obj.nsats+j));
                            end
                        elseif obj.rate
                            obj.ymeas(j,i) = dr + obj.c*x_u(8) + mvnrnd(0,obj.R(j,j));
                        end

                        % scatter3(x_s(1), x_s(2), x_s(3), 'g', '+');
                    end
                end

                % hold off; axis equal; grid on;
            end
        end
    end
end

% custom validation
function mustBeOption(txt)
%MUSTBEOPTION Tests that string option is allowable

if ~strcmpi(txt,"RANGE") && ~strcmpi(txt,"RATE") && ~strcmpi(txt,"BOTH")
    error("GNSSmeasurements:obsType", ...
                    "Observation type obs must be 'RANGE', 'RATE', or 'BOTH'.")
end
end