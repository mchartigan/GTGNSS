classdef NavSatellite < handle
    %NavSatellite Class for describing the properties and trajectory of a
    %satellite that is part of a radionavigation satellite system (whether
    %that's earth-based GNSS or lunar).
    
    properties
        % propagator instance (default one so MATLAB doesn't throw a fit)
        prop    (1,1)   OrbitPropagator = OrbitPropagator(0,1)
        % nav filter (default one so MATLAB doesn't throw a fit)
        filter  (1,1)   EKF = EKF("discrete", @(~) 0, @(~) 0, [], @(~) 0, @(~) 0, @(~) 0, [], [])
    end
    
    methods
        function obj = NavSatellite(prop,filter)
            %NAVSATELLITE Construct a NavSatellite instance.
            %   Inputs:
            %    - prop; propagator instance for generating future trajectories
            %            (OrbitPropagator, LunarPropagator, EarthPropagator)
            %    - filter; navigation filter for finding s/c uncertainty (EKF)

            obj.prop = prop;
            obj.filter = filter;
        end
        
        function xs = getrefstates(obj,ts,frame)
            %GETREFSTATES Returns reference trajectory information for the
            %satellite at the provided times and frame. Also referred to as
            %the true state.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - frame; reference frame to return data in
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                frame   (1,:)   char
            end

            if ts(1) > obj.prop.t0
                error("getrefstates:invalidInput", ...
                    "Simulation times must be after starting epoch of propagator.");
            end

            ts = ts - obj.prop.t0;
            [~,xs] = obj.prop.runat(ts,frame);
        end

        function [xs,Ps] = getnavstates(obj,ts,x0,P0,frame)
            %GETNAVSTATES Returns estimated state and covariance from the 
            %navigation filter at the provided times and frame.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - x0; state estimate at ts(1) in J2000 frame
            %    - P0; uncertainty of x0 at ts(1) in J2000 frame
            %    - frame; reference frame to return data in
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                x0      (9,1)   double
                P0      (9,9)   double
                frame   (1,:)   char
            end

            % add ts to obj.filter.t with union() (maybe strip everything
            % before ts(1)?)
            % run filter, get data at time steps, and return it
        end

        function [] = getmodelstates(obj,ts,frame)
            %GETMODELSTATES description here
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                frame   (1,:)   char
            end
    
            % generally should be given starting epoch, validity window. take
            % that and update the propagator and use ephemerisfit. Then take
            % the returned function and get state estimate at provided ts.
            % Transform frame if necessary.

            % alternatively, be given a whole window and call getnavstates,
            % compute the model states with sliding windows scaled to meet
            % some tolerance (maybe provided?). Return that, along w/
            % getnavstates output if desired to break down error budget.
        end
    end
end

