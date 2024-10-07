classdef EKF < handle
    %EKF Standard representation of a [discrete/hybrid] extended Kalman 
    %filter with [discrete/continuous] dynamics and discrete measurements.
    
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
        m      (1,1)   {mustBeNonnegative, mustBeInteger} = 1   % # of steps in t_meas
        s      (1,1)   {mustBePositive, mustBeInteger} = 1      % # of steps in t
        type   (1,1)   {mustBeOption} = "hybrid"                % EKF dynamics type
        f      (1,1)   function_handle = @(varargin) disp([])   % system dynamics
        h      (1,1)   function_handle = @(varargin) disp([])   % measurement model
        dfdx   (1,1)   function_handle = @(varargin) disp([])   % Jacobian of system dynamics
        dhdx   (1,1)   function_handle = @(varargin) disp([])   % Jacobian of measurement model
    end
    
    methods
        function obj = EKF(type, f, dfdx, Q, h, y, dhdx, R, t_meas, varargin)
            %EKF Construct an EKF instance (either discrete or hybrid).
            %   Inputs:
            %    - type; either "discrete" or "hybrid"
            %    - f; system dynamics (continuous, integrable w/ ODE45 or discrete)
            %    - dfdx; Jacobian of dynamics w.r.t. state
            %    - Q; process noise matrix (continuous or discrete, must be same as dyn)
            %    - h; measurement model
            %    - y; measurements
            %    - dhdx; Jacobian of measurement w.r.t. state
            %    - R; measurement noise matrix (discrete)
            %    - t_meas; time stamps where measurements are taken
            %    - t_sim; optional (if different from t_meas) name-value
            %             pair, time stamps to get state between measurements
            %    - opts; optional name-value pair, ODE45 propagation options

            % default values for things
            obj.t = t_meas;
            obj.t_meas = t_meas;
            obj.opts = odeset();

            if nargin > 9
                for i=1:2:nargin-9
                    if strcmp(varargin{i}, "t_sim")
                        obj.t = union(t_meas, varargin{i+1});
                    elseif strcmp(varargin{i}, "opts")
                        obj.opts = varargin{i+1};
                    end
                end
            elseif nargin < 9
                error("EKF:nargin", "Too few arguments.")
            end

            % check if I/O is correct (nargout == -1 acceptable, means
            % variable which is usually from anonymous funcs)
            if (nargin(f) ~= 2 && nargin(f) ~= -1) || (nargout(f) ~= 1 && nargout(f) ~= -1)
                error("EKF:functionIO", ...
                    "f must accept (t,x) and return (dyn).")
            end
            if (nargin(dfdx) ~= 2 && nargin(dfdx) ~= -1) || (nargout(dfdx) ~= 1 && nargout(dfdx) ~= -1)
                error("EKF:functionIO", ...
                    "dfdx must accept (t,x) and return (dfdx).")
            end
            if (nargin(h) ~= 2 && nargin(h) ~= -1) || (nargout(h) ~= 1 && nargout(h) ~= -1)
                error("EKF:functionIO", ...
                    "h must accept (t,x) and return (y).")
            end
            if (nargin(dhdx) ~= 2 && nargin(dhdx) ~= -1) || (nargout(dhdx) ~= 1 && nargout(dhdx) ~= -1)
                error("EKF:functionIO", ...
                    "dhdx must accept (t,x) and return (dhdx).")
            end
            % check number of supplied measurements is correct
            if size(y,2) ~= length(t_meas)
                error("EKF:measurementNum", ...
                    "# of columns of y must equal length of t_meas.")
            end

            obj.type = lower(type);
            obj.f = f;
            obj.dfdx = dfdx;
            obj.h = h;
            obj.dhdx = dhdx;
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
            %RUN Execute simulation of EKF, provided initial conditions.
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
            
            for k=2:obj.s
                % store next state/cov matrix to be included in filter
                [x_, P_] = obj.propagate(k-1);
                tk = obj.t(k);
            
                if ismember(tk, obj.t_meas)   % step with measurement
                    j = find(obj.t_meas == obj.t(k));
                
                    H = obj.dhdx(tk,x_);                    % measurement partials matrix
                    K = P_*H' / (H*P_*H' + obj.R);          % Kalman gain
                    yj = obj.y(:,j);                        % get state measurement
                    obj.x(:,k) = x_ + K*(yj - obj.h(tk,x_));% post-fit state estimate
                    % post-fit est. error covariance
                    obj.P(:,:,k) = (eye(obj.n) - K*H)*P_;   %*(eye(obj.n)-K*H)' + K*obj.R*K';
                else                                % step without measurement
                    obj.x(:,k) = x_;
                    obj.P(:,:,k) = P_;
                end
            end
        end

        function [x,P] = propagate(obj, k0)
            %PROPDYNAMICS Propagates dynamics (obj.f) and lyapunov equation
            %(obj.Lyapunov) one time step and provides final state.
            %   Input:
            %    - obj; self, EKF object
            %    - k0; starting index
            arguments
                obj (1,1) EKF
                k0  (1,1) {mustBePositive,mustBeInteger}
            end

            if strcmp(obj.type, "hybrid")   % continuous dynamics
                % reshape starting P to correct format
                P0 = reshape(obj.P(:,:,k0), obj.n*obj.n, 1);
                v0 = [obj.x(:,k0); P0];
                [~,V] = ode45(@(t,v) obj.fulldyn(t, v, obj.dfdx(t,obj.x(:,k0)), ...
                    obj.Q), obj.t([k0 k0+1]), v0, obj.opts);
    
                % store covariance matrices in appropriate structure
                x = V(end,1:obj.n)';
                P = reshape(V(end,obj.n+1:end)', obj.n, obj.n);
            else                            % discrete dynamics
                x = obj.f(obj.t([k0 k0+1]), obj.x(:,k0));
                F = obj.dfdx(obj.t([k0 k0+1]), obj.x(:,k0));
                P = F*obj.P(:,:,k0)*F' + obj.Q;
            end
        end

        function x = propDynamics(obj, k0)
            %PROPDYNAMICS Propagates dynamics (obj.f) from given to next 
            %time and provides final state.
            %   Input:
            %    - obj; self, EKF object
            %    - k0; starting index
            arguments
                obj (1,1) EKF
                k0  (1,1) {mustBePositive,mustBeInteger}
            end
            
            if strcmp(obj.type, "hybrid")   % continuous dynamics
                [~,X] = ode45(obj.f, obj.t([k0 k0+1]), obj.x(:,k0), obj.opts);
                x = X(end,:)';
            else                            % discrete dynamics
                x = obj.f(obj.t([k0 k0+1]), obj.x(:,k0));
            end
        end

        function P = propLyapunov(obj, k0)
            %PROPLYAPUNOV Propagates Lyapunov equations (obj.lyapunov) from
            %given to next time and provides final covariance matrix.
            %   Input:
            %    - obj; self, EKF object
            %    - k0; starting index
            arguments
                obj (1,1) EKF
                k0  (1,1) {mustBePositive,mustBeInteger}
            end
            
            if strcmp(obj.type, "hybrid")   % continuous dynamics
                % reshape starting P to correct format
                P0 = reshape(obj.P(:,:,k0), obj.n*obj.n, 1);
                [~,Y] = ode45(@(t,p) obj.lyapunov(p, obj.dfdx(t,obj.x(:,k0)), ...
                    obj.Q), obj.t([k0 k0+1]), P0, obj.opts);
    
                % store covariance matrices in appropriate structure
                P = reshape(Y(end,:)', obj.n, obj.n);
            else                            % discrete dynamics
                F = obj.dfdx(obj.t([k0 k0+1]), obj.x(:,k0));
                P = F*obj.P(:,:,k0)*F' + obj.Q;
            end
        end

        function dv = fulldyn(obj,t,v,A,Q)
            %FULLDYN Describes the joint state dynamics and continuous-time
            %Lyapunov equations for a hybrid EKF. For use with ODE45 or
            %other propagator.
            %
            %   Inputs: (dims)
            %    - obj; self, EKF object
            %    - t; time step
            %    - v; vectorized full state
            %    - A; (nxn) linearized state dynamics matrix
            %    - Q; (nxn) process noise matrix (may also take form LQL')

            x_ = v(1:obj.n);
            P_ = v(obj.n+1:end);
            dv = [obj.f(t,x_); obj.lyapunov(P_,A,Q)];
        end

        function dP = lyapunov(obj,P,A,Q)
            %LYAPUNOV Describes the continuous-time Lyapunov equations for 
            %a hybrid EKF. For use with ODE45 or other propagator.
            %   Takes the form
            %       dP = AP + PA' + Q
            %   P must be passed in as a column vector rather than a matrix
            %   to play nice with ODE45.
            %
            %   Inputs: (dims)
            %    - obj; self, EKF object
            %    - P; ((nxn)x1) Vectorized covariance matrix
            %    - A; (nxn) linearized state dynamics matrix
            %    - Q; (nxn) process noise matrix (may also take form LQL')
            arguments
                obj (1,1) EKF
                P   (:,1) double
                A   (:,:) double {mustBeNxN(obj,A)}
                Q   (:,:) double {mustBeNxN(obj,Q)}
            end

            % column vector to square matrix
            P = reshape(P, obj.n, obj.n);   
            dP = A*P + P*A' + Q;            % Lyapunov equation
            % square matrix to column vector
            dP = reshape(dP, obj.n*obj.n, 1);   
        end
    end
end

% custom validation
function mustBeOption(txt)
%MUSTBEOPTION Tests that string option is allowable

if ~sum(contains(["discrete","hybrid"], txt))
    error("EKFtype:notOption", "EKF type must be 'discrete' or 'hybrid'.")
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
