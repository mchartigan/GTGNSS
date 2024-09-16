classdef UKF < handle
    %UKF Standard representation of an unscented Kalman filter with
    %[continuous/discrete] dynamics and discrete measurements.
    
    properties
        t      (1,:)   double               % simulation times
        t_meas (1,:)   double               % measurement times
        x      (:,:)   double               % state history
        y      (:,:)   double               % measurements
        P      (:,:,:) double               % state covariance history
        Q      (:,:)   double               % process noise
        R      (:,:)   double               % measurement noise
        opts   (1,1)   struct               % ODE45 propagation options
        n      (1,1)   {mustBePositive, mustBeInteger} = 1      % # of states
        m      (1,1)   {mustBePositive, mustBeInteger} = 1      % # of steps in t_meas
        s      (1,1)   {mustBePositive, mustBeInteger} = 1      % # of steps in t
        type   (1,1)   {mustBeOption} = "hybrid"                % UKF dynamics type
        f      (1,1)   function_handle = @(varargin) disp([])   % system dynamics
        h      (1,1)   function_handle = @(varargin) disp([])   % measurement model
    end
    
    methods
        function obj = UKF(type, f, Q, h, y, R, t_meas, varargin)
            %UKF Construct a UKF instance.
            %   Inputs:
            %    - type; either "discrete" or "hybrid"
            %    - f; system dynamics (continuous, integrable w/ ODE45)
            %    - Q; process noise matrix (continuous)
            %    - h; measurement model
            %    - y; measurements
            %    - R; measurement noise matrix (discrete)
            %    - t_meas; time stamps where measurements are taken
            %    - t_sim; optional (if different from t_meas) name-value
            %             pair, time stamps to get state between measurements
            %    - opts; optional name-value pair, ODE45 propagation options

            % default values for things
            obj.t = t_meas;
            obj.t_meas = t_meas;
            obj.opts = odeset();

            if nargin > 7
                for i=1:2:nargin-7
                    if strcmp(varargin{i}, "t_sim")
                        obj.t = union(t_meas, varargin{i+1});
                    elseif strcmp(varargin{i}, "opts")
                        obj.opts = varargin{i+1};
                    end
                end
            elseif nargin < 7
                error("UKF:nargin", "Too few arguments.")
            end

            % check if I/O is correct (nargout == -1 acceptable, means
            % variable which is usually from anonymous funcs)
            if nargin(f) ~= 2 || (nargout(f) ~= 1 && nargout(f) ~= -1)
                error("UKF:functionIO", ...
                    "f must accept (t,x) and return (dyn).")
            end
            if nargin(h) ~= 2 || (nargout(h) ~= 1 && nargout(h) ~= -1)
                error("UKF:functionIO", ...
                    "h must accept (t,x) and return (y).")
            end
            if size(y,2) ~= length(t_meas)
                error("UKF:measurementNum", ...
                    "# of columns of y must equal length of t_meas.")
            end

            obj.type = lower(type);
            obj.f = f;
            obj.h = h;
            obj.y = y;
            obj.Q = Q;
            obj.R = R;

            obj.n = size(Q,1);
            obj.m = length(t_meas);
            obj.s = length(obj.t);

            obj.x = zeros(obj.n, obj.s);
            obj.P = zeros(obj.n, obj.n, obj.s);
        end
        
        function run(obj, x0, P0)
            %RUN Execute simulation of UKF, provided initial conditions.
            %   Inputs:
            %    - x0; initial state
            %    - P0; initial state covariance
            arguments
                obj
                x0  (:,1) double {mustBeNx1(obj,x0)}
                P0  (:,:) double {mustBeNxN(obj,P0)}
            end

            % initialize state info
            obj.x(:,1)   = x0;
            obj.P(:,:,1) = P0;
            ni = obj.n;
            
            for k=2:obj.s
                % Time Update
                % (1) choose sigma points
                xis = obj.getSigmaPoints(obj.x(:,k-1), obj.P(:,:,k-1));
                % (2) Propagate sigma points forward using model
                for i=1:2*ni
                    xis(:,i) = obj.propDynamics(k-1, xis(:,i));
                end
                % (3) Generate a priori estimate of x
                x_ = 1/(2*ni) * sum(xis, 2);
                % (4) Generate a priori estimate of P
                xis = xis - x_;
                if strcmp(obj.type, "hybrid")       % continuous dynamics
                    Qi = obj.Q*(obj.t(k)-obj.t(k-1));
                else                                % discrete dynamics
                    Qi = obj.Q;
                end
                P_ = 1/(2*ni) * (xis*xis') + Qi;

                if ismember(obj.t(k), obj.t_meas)   % step with measurement
                    j = find(obj.t_meas == obj.t(k));
                    % Measurement Update
                    % (1) Choose sigma points again
                    xis = obj.getSigmaPoints(x_, P_);
                    x_ = 1/(2*ni) * sum(xis, 2);
                    % (2) Generate measurements from sigma points
                    yis = zeros(size(obj.y,1), 2*ni);
                    for i=1:2*ni
                        yis(:,i) = obj.h(obj.t(k),xis(:,i));
                    end
                    % (3) combine to create predicted measurement
                    y_ = 1/(2*ni) * sum(yis, 2);
                    % (4) Generate measurement variance
                    yis = yis - y_;
                    Py = 1/(2*ni) * (yis*yis') + obj.R;
                    % (5) Estimate cross covariance
                    xis = xis - x_;
                    Pxy = 1/(2*ni) * xis*yis';
                    % (6) Compute gain
                    K = Pxy / Py;
                    % (7) Compute measurement update
                    yj = obj.y(:,j);
                    obj.x(:,k) = x_ + K*(yj - y_);
                    obj.P(:,:,k) = P_ - K*Py*K';
                else                                % step without measurement
                    obj.x(:,k) = x_;
                    obj.P(:,:,k) = P_;
                end
            end
        end

        function xi = getSigmaPoints(obj, x, P)
            %GETSIGMAPOINTS Generate 2N sigma points around the given state
            %and covariance.
            %    - obj; self, UKF object
            %    - k; index to generate about
            arguments
                obj (1,1) UKF
                x   (:,1) double {mustBeNx1(obj,x)}
                P   (:,:) double {mustBeNxN(obj,P)}
            end

            ni = obj.n;
            xi = repmat(x, 1, 2*ni);    % vector of 2N sigma points
            try
                rootNP = chol(ni*P)';
            catch
                0;
            end

            for i=1:ni
                xi(:,i) = xi(:,i) + rootNP(:,i);
                xi(:,ni+i) = xi(:,ni+i) - rootNP(:,i);
            end
        end

        function x = propDynamics(obj, k0, xi)
            %PROPDYNAMICS Propagates dynamics (obj.f) from given to next 
            %time and provides final state.
            %   Input:
            %    - obj; self, UKF object
            %    - k0; starting index
            %    - xi; optional different starting state
            arguments
                obj (1,1) UKF
                k0  (1,1) {mustBePositive,mustBeInteger}
                xi  (:,1) double {mustBeNx1(obj,xi)} = obj.x(:,k0)
            end

            if strcmp(obj.type, "hybrid")   % continuous dynamics
                [~,X] = ode45(obj.f, obj.t([k0 k0+1]), xi, obj.opts);
                x = X(end,:)';
            else                            % discrete dynamics
                x = obj.f(obj.t([k0 k0+1]), xi);
            end
        end
    end
end

% custom validation
function mustBeOption(txt)
%MUSTBEOPTION Tests that string option is allowable

if ~sum(contains(["discrete","hybrid"], txt))
    error("UKFtype:notOption", "UKF type must be 'discrete' or 'hybrid'.")
end
end

function mustBeNx1(obj,x)
%MUSTBENX1 Tests that state is nx1

if size(x,1) ~= obj.n || size(x,2) ~= 1
    error("vector:wrongSize", "Vector must be N x 1.")
end
end

function mustBeNxN(obj,x)
%MUSTBENXN Tests that state is nxn

if size(x,1) ~= obj.n || size(x,2) ~= obj.n
    error("matrix:wrongSize", "Matrix must be N x N.")
end
end
