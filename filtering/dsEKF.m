classdef dsEKF < handle
    %DSEKF Construction of a delayed-state extended Kalman filter with
    %specified dynamics and discrete measurements.
    
    properties
        t      (1,:)   double               % simulation times
        t_meas (1,:)   double               % measurement times
        x      (:,:)   double               % state history
        y      (:,:)   double               % measurements
        R      (:,:,:) double               % measurement noise
        P      (:,:,:) double               % state covariance history
        U      (:,:)   double               % measurement underweighting
        opts   (1,1)   struct               % ODE45 propagation options
        n      (1,1)   {mustBePositive, mustBeInteger} = 1      % # of states
        m      (1,1)   {mustBeNonnegative, mustBeInteger} = 1   % # of steps in t_meas
        s      (1,1)   {mustBePositive, mustBeInteger} = 1      % # of steps in t
        prop   (1,1)   Propagator = OrbitPropagator(1)          % describes system dynamics
        meas   (1,1)   Measurement = EmptyMeasurement()         % describes system<->measurement interface

        % DEBUG ONLY PROPERTY %
        % truth trajectory so we can compare during execution
        truth   (1,1)   User
    end
    
    methods
        function obj = dsEKF(prop,meas,y,R,t_meas,options)
            %DSEKF Construct a dsEKF instance (either discrete or hybrid).
            %   Inputs:
            %    - prop; dynamics propagator instance (inherits Propagator)
            %    - meas; measurement instance
            %    - y; measurements
            %    - R; measurement covariance matrix
            %    - t_meas; time stamps where measurements are taken
            %    - t_sim; optional (if different from t_meas) name-value
            %             pair, time stamps to get state between measurements
            %    - opts; optional name-value pair, ODE45 propagation options
            arguments
                prop    (1,1)   Propagator
                meas    (1,1)   Measurement
                y       (:,:)   double
                R       (:,:,:) double
                t_meas  (1,:)   double
                options.t_sim   (1,:)   double = []
                options.opts    (1,1)   struct = odeset()
            end

            % check number of supplied measurements is correct
            if size(y,2) ~= length(t_meas)
                error("dsEKF:measurementNum", ...
                    "# of columns of y must equal length of t_meas.")
            end

            % assign values
            obj.prop = prop;
            obj.meas = meas;
            obj.y = y;
            obj.R = R;
            obj.t = union(t_meas, options.t_sim);
            obj.t_meas = t_meas;

            obj.n = prop.dim;
            obj.m = length(t_meas);
            obj.s = length(obj.t);
            obj.x = zeros(obj.n, obj.s);
            obj.P = zeros(obj.n, obj.n, obj.s);
        end
        
        function run(obj,x0,P0)
            %RUN Execute simulation of EKF, provided initial conditions.
            %   Inputs:
            %    - x0; initial state
            %    - P0; initial state covariance
            arguments
                obj (1,1) dsEKF
                x0  (:,1) double {mustBeNx1(obj,x0)}
                P0  (:,:) double {mustBeNxN(obj,P0)}
            end

            % initialize state info
            obj.x(:,1)   = x0;
            obj.P(:,:,1) = P0;
            Ak_1 = obj.prop.partials(obj.t(1), x0);
            tprev = obj.t(1);
            xprev = obj.x(:,1);
            Pprev = obj.P(:,:,1);

            for k=2:obj.s
                tk = obj.t(k);
                dt = tk - obj.t(k-1);
                % store next state/cov matrix to be included in filter
                x_ = obj.prop.runat([obj.t(k-1) tk], obj.x(:,k-1));
                x_ = x_(:,2);
                Ak = obj.prop.partials(tk, x_);
                % compute STM using 2nd-order Runge-Kutta (one of Lear's methods)
                Phi = eye(obj.n) + (Ak_1 + Ak)/2*dt + Ak_1*Ak*dt^2/2;
                % get DT process noise
                S = obj.prop.noise(dt, obj.x(:,k-1));
                % a-priori state covariance
                P_ = Phi*obj.P(:,:,k-1)*Phi' + S;
            
                if ismember(tk, obj.t_meas)     % step with measurement
                    j = find(obj.t_meas == obj.t(k));
                
                    yj = obj.y(:,j);                % get state measurement
                    mask = ~isnan(yj);              % generate mask of any NaN

                    if sum(mask) > 0
                        0;
                    end
                    % get computed measurement from Measurement model
                    ycomp = obj.meas.computemeas(tk,x_,tprev,xprev);
                    Y = yj - ycomp;                 % measurement residual (O - C)
                    Y = Y(mask);                    % mask out invalid meas
                    % measurement partials matrix
                    [H,J] = obj.meas.measpartials(tk,x_,tprev,xprev);
                    H = H(mask,:);                  % mask out invalid meas
                    J = J(mask,:);
                    % get appropriate measurement noise
                    if size(obj.R, 3) > 1       % time-varying
                        Rk = obj.R(:,:,j);
                    else                        % time-invariant
                        Rk = obj.R(:,:);
                    end

                    % underweight the pseudorange measurements
                    % ns = length(mask);
                    % Rk(mask(1:ns/3),mask(1:ns/3)) = Rk(mask(1:ns/3),mask(1:ns/3));
                    % Rk(mask(ns/3+1:end),mask(ns/+1:end)) = Rk(mask(ns/3+1:end),mask(ns/+1:end));
                    Rk = Rk(mask,mask);
                    
                    % post-fit est. error covariance
                    L = H*P_*H' + Rk + J*Pprev*Phi'*H' + H*Phi*Pprev*J' + J*Pprev*J';
                    K = (P_*H' + Phi*Pprev*J') / L; % Kalman gain
                    obj.x(:,k) = x_ + K*Y;          % post-fit state estimate
                    obj.P(:,:,k) = P_ - K*L*K';

                    % store states for next time
                    tprev = tk;
                    xprev = obj.x(:,k);
                    Pprev = obj.P(:,:,k);

                else                            % step without measurement
                    obj.x(:,k) = x_;
                    obj.P(:,:,k) = P_;
                end
            end
        end
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
