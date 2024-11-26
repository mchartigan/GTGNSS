classdef LUMVE < handle
    %LUMVE Standard representation of a linear unbiased minimum variance
    %estimator (LUMVE) with [discrete/continuous] dynamics and discrete
    %measurements.
    
    properties
        t      (1,:)   double               % simulation times
        t_meas (1,:)   double               % measurement times
        x      (:,:)   double               % state history
        y      (:,:)   double               % measurements
        P0     (:,:)   double               % initial state covariance
        Q      (:,:)   double               % process noise
        R      (:,:,:) double               % measurement noise
        opts   (1,1)   struct               % ODE45 propagation options
        n      (1,1)   {mustBePositive, mustBeInteger} = 1      % # of states
        m      (1,1)   {mustBeNonnegative, mustBeInteger} = 1      % # of steps in t_meas
        s      (1,1)   {mustBePositive, mustBeInteger} = 1      % # of steps in t
        iter   (1,1)   {mustBePositive, mustBeInteger} = 1      % # of LUMVE iterations
        type   (1,1)   {mustBeOption} = "hybrid"                % dynamics type
        f      (1,1)   function_handle = @(varargin) disp([])   % system dynamics
        h      (1,1)   function_handle = @(varargin) disp([])   % measurement model
        dfdx   (1,1)   function_handle = @(varargin) disp([])   % Jacobian of system dynamics
        dhdx   (1,1)   function_handle = @(varargin) disp([])   % Jacobian of measurement model
    end
    
    methods
        function obj = LUMVE(type, n, f, dfdx, h, y, dhdx, R, iter, t_meas, varargin)
            %LUMVE Construct a LUMVE instance (either discrete or hybrid).
            %   Inputs:
            %    - type; either "discrete" or "hybrid"
            %    - n; dimension of state
            %    - f; system dynamics (continuous, integrable w/ ODE45 or discrete)
            %    - dfdx; Jacobian of dynamics w.r.t. state
            %    - h; measurement model
            %    - y; measurements
            %    - dhdx; Jacobian of measurement w.r.t. state
            %    - R; measurement noise matrix (discrete)
            %    - iter; number of iterations for LUMVE to take
            %    - t_meas; time stamps where measurements are taken
            %    - t_sim; optional (if different from t_meas) name-value
            %             pair, time stamps to get state between measurements
            %    - opts; optional name-value pair, ODE45 propagation options

            % default values for things
            obj.t = t_meas;
            obj.t_meas = t_meas;
            obj.opts = odeset();

            if nargin > 10
                for i=1:2:nargin-10
                    if strcmp(varargin{i}, "t_sim")
                        obj.t = union(t_meas, varargin{i+1});
                    elseif strcmp(varargin{i}, "opts")
                        obj.opts = varargin{i+1};
                    end
                end
            elseif nargin < 10
                error("LUMVE:nargin", "Too few arguments.")
            end

            % check if I/O is correct (nargout == -1 acceptable, means
            % variable which is usually from anonymous funcs)
            if (nargin(f) ~= 2 && nargin(f) ~= -1) || (nargout(f) ~= 1 && nargout(f) ~= -1)
                error("LUMVE:functionIO", ...
                    "f must accept (t,x) and return (dyn).")
            end
            if (nargin(dfdx) ~= 2 && nargin(dfdx) ~= -1) || (nargout(dfdx) ~= 1 && nargout(dfdx) ~= -1)
                error("LUMVE:functionIO", ...
                    "dfdx must accept (t,x) and return (dfdx).")
            end
            if nargin(h) ~= 2 || (nargout(h) ~= 1 && nargout(h) ~= -1)
                error("LUMVE:functionIO", ...
                    "h must accept (t,x) and return (y).")
            end
            if nargin(dhdx) ~= 2 || (nargout(dhdx) ~= 1 && nargout(dhdx) ~= -1)
                error("LUMVE:functionIO", ...
                    "dhdx must accept (t,x) and return (dhdx).")
            end
            % check number of supplied measurements is correct
            if size(y,2) ~= length(t_meas)
                error("LUMVE:measurementNum", ...
                    "# of columns of y must equal length of t.")
            end

            obj.type = lower(type);
            obj.f = f;
            obj.dfdx = dfdx;
            obj.h = h;
            obj.dhdx = dhdx;
            obj.y = y;
            obj.R = R;
            obj.iter = iter;

            % dimensions of various items
            obj.n = n;
            obj.m = length(t_meas);
            obj.s = length(obj.t);

            % initialize state dimensions
            obj.x = zeros(obj.n, obj.s);
            obj.P0 = zeros(obj.n, obj.n);
        end
        
        function run(obj, x0, P0)
            %RUN Execute simulation of LUMVE, provided initial conditions
            %(or not).
            %   Inputs:
            %    - x0; optional, initial state
            %    - P0; optional, initial state covariance
            arguments
                obj
                x0 (:,1) double {mustBeNx1(obj,x0)} = zeros(obj.n,1);
                P0 (:,:) double {mustBeNxN(obj,P0)} = zeros(obj.n,obj.n);
            end
            
            % check appropriate number of inputs supplied
            if nargin ~= 1 && nargin ~= 3
                error("run:nargin", "Either provide both initial conditions or neither.")
            end

            
            % initialize state info
            obj.x   = repmat(x0, 1, obj.s);
            xbar    = zeros(size(x0));   % state deviation
            obj.P0  = P0;
            
            for i=1:obj.iter    % iterate specified # of times
                STM = eye(obj.n);

                % use a-priori info if available
                if i ~= 1 || nargin == 3
                    HTWH = inv(obj.P0);
                    HTWy = obj.P0 \ xbar;
                else
                    HTWH = zeros(obj.n, obj.n);
                    HTWy = zeros(obj.n, 1);
                end

                % initialize Aprev
                Aprev = obj.dfdx(obj.t(1), obj.x(:,1));

                for k=2:obj.s       % iterate through simulation times
                    % store next nominal state if we have an x0
                    if i ~= 1 || nargin == 3
                        obj.x(:,k) = obj.propDynamics(k-1);
                    end

                    tk = obj.t(k);
                    dt = tk - obj.t(k-1);       % time step size

                    x_ = obj.x(:,k);            % computed state

                    % from current state estimate, get STM using Lear's
                    A = obj.dfdx(tk, x_);       % state matrix
                    % STM = numericallyGetNextSTM(STM, A, dtk, obj.opts);
                    % propogate STM using Lear's method, as recommended by
                    % NASA Filter Best Practice eqn H
                    STM_tkm1_tk = eye(obj.n) + (A + Aprev)*dt/2 ...
                        + sign(dt)*(A*Aprev)*dt^2/2 ;
                    STM = STM_tkm1_tk*STM;
                    % update previous jacobian
                    Aprev = A;
                
                    if ismember(tk, obj.t_meas)     % step with measurement
                        j = find(obj.t_meas == tk);

                        yj = obj.y(:,j);            % observed measurement
                        mask = ~isnan(yj);          % generate mask of any NaN
                        Y = yj - obj.h(tk, x_);     % measurement residual (O - C)
                        Y = Y(mask);                % remove NaN
                        H_ = obj.dhdx(tk, x_);      % get full measurement matrix
                        H = H_(mask,:) * STM;       % measurement matrix (mask NaN)
                
                        % compute rank-1 update
                        R_ = obj.R(mask,mask,j);    % mask out NaN meas
                        HTWH = HTWH + H'*(R_\H);
                        HTWy = HTWy + H'*(R_\Y);
                    end
                end


                obj.P0 = inv(HTWH);
                dX = obj.P0 * HTWy;             % update to state deviation

                % Update nominal trajectory and deviation vector xbar --
                % note that we require nominal trajectory obj.x(:,1) + xbar
                % to be a constant for all iterations
                obj.x(:,1) = obj.x(:,1) + dX;   % get final estimate
                xbar = xbar - dX;
            end

            for k=2:obj.s       % get final state propagation
                obj.x(:,k) = obj.propDynamics(k-1);
            end
        end

        function x = propDynamics(obj, k0)
            %PROPDYNAMICS Propagates dynamics (obj.f) from given to next 
            %time and provides final state.
            %   Input:
            %    - obj; self, LUMVE object
            %    - k0; starting index
            arguments
                obj (1,1) LUMVE
                k0  (1,1) {mustBePositive,mustBeInteger}
            end
            
            if strcmp(obj.type, "hybrid")   % continuous dynamics
                [~,X] = ode45(obj.f, obj.t([k0 k0+1]), obj.x(:,k0), obj.opts);
                x = X(end,:)';
            else                            % discrete dynamics
                x = obj.f(obj.t([k0 k0+1]), obj.x(:,k0));
            end
        end
    end
end

% custom validation
function mustBeOption(txt)
%MUSTBEOPTION Tests that string option is allowable

if ~sum(contains(["discrete","hybrid"], txt))
    error("LUMVEtype:notOption", "LUMVE type must be 'discrete' or 'hybrid'.")
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


