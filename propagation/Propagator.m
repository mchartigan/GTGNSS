classdef (Abstract) Propagator < handle
    %PROPAGATOR Abstract class for defining various types of propagators.
    %This applies most to orbital propagators and clock propagators.
    
    properties (Abstract)
        dim (1,1)   {mustBeInteger,mustBePositive}  % dimension of state
    end

    methods
        function dP = lyapunov(obj,P,A,Q)
            %LYAPUNOV Describes the continuous-time Lyapunov equations. 
            %For use with ODE45 or other propagator.
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
                obj (1,1) Propagator
                P   (:,1) double
                A   (:,:) double
                Q   (:,:) double
            end

            % column vector to square matrix
            P = reshape(P, obj.dim, obj.dim);   
            dP = A*P + P*A' + Q;            % Lyapunov equation
            % square matrix to column vector
            dP = reshape(dP, obj.dim*obj.dim, 1);   
        end
    end
    
    methods (Abstract)
        out = run(obj,ts,x0,n)      % run to final time w/ equispaced steps
        out = runat(obj,ts,x0)      % run to specific times
        out = modelfit(obj)         % fit an approximation model to the result
        dxdt = dynamics(obj,t,x)    % partial differential equations governing dynamics
        A = partials(obj,t,x)       % jacobian of dynamics
        P = proplyapunov(obj,ts,P0) % function to propagate lyapunov equations
    end
end

