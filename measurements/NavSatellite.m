classdef NavSatellite < handle
    %NavSatellite Class for describing the properties and trajectory of a
    %satellite that is part of a radionavigation satellite system (whether
    %that's earth-based GNSS or lunar).
    
    properties
        % propagator instance (default one so MATLAB doesn't throw a fit)
        prop    (1,1)   OrbitPropagator = OrbitPropagator(0,1)
        % clock instance (default one again ^)
        clock   (1,1)   Clock = Clock(0,zeros(3,1),"none")
        % nav filter
        filter  (1,1)
        % s/c nav default uncertainty
        P0      (9,9)   double
    end
    
    methods
        function obj = NavSatellite(prop,clock,filter,meas)
            %NAVSATELLITE Construct a NavSatellite instance.
            %   Inputs:
            %    - prop; propagator instance for generating future trajectories
            %            (OrbitPropagator, LunarPropagator, EarthPropagator)
            %    - clock; instance of s/c oscillator/clock (Clock)
            %    - filter; navigation filter type for finding s/c uncertainty,
            %              options are: "EKF"
            %    - meas; struct containing measurement info for the chosen
            %            filter -- names should correspond to filter constructor
            %            arguments (may also include varargin)
            arguments
                prop    (1,1)   OrbitPropagator
                clock   (1,1)   Clock
                filter  (1,:)   {mustBeText}
                meas    (1,1)   struct
            end

            % assign instances
            obj.prop = prop;
            obj.clock = clock;

            if strcmp(filter, "EKF")
                % CT process noise
                Q = [obj.prop.var zeros(6,3); zeros(3,6) obj.clock.var];

                % build filter, w/ or w/o optional args
                if isfield(meas, "varargin")
                    obj.filter = EKF("hybrid", @obj.dynamics, @obj.partials, ...
                        Q, @meas.h, meas.y, @meas.dhdx, meas.R, meas.t_meas, ...
                        meas.varargin);
                else
                    obj.filter = EKF("hybrid", @obj.dynamics, @obj.partials, ...
                        Q, @meas.h, meas.y, @meas.dhdx, meas.R, meas.t_meas);
                end
            else
                error("NavSatellite:invalidFilter", ...
                    'Filter must be "EKF". See documentation.');
            end
        end

        function dxdt = dynamics(obj,t,x)
            %DYNAMICS Returns the full-state dynamics for a NavSatellite
            %instance. Incorporates both orbital and oscillator dynamics.
            %   Input:
            %    - t; seconds past J2000
            %    - x; s/c state
            arguments
                obj (1,1)   NavSatellite
                t   (1,1)   double
                x   (9,1)   double
            end

            dxdt = [obj.prop.dynamics(t,x(1:6)); obj.clock.dynamics(t,x(7:9))];
        end

        function A = partials(obj,t,x)
            %PARTIALS Invokes the partials of dynamics, with relevant settings, 
            %for this propagator instance. Incorporates both orbital and 
            %oscillator dynamics.
            %   Input:
            %    - t; simulation time in seconds past J2000
            %    - x; satellite state
        
            A1 = obj.prop.partials(t,x(1:6));
            A2 = obj.clock.partials(t,x(7:9));
            A = [A1 zeros(6,3); zeros(3,6) A2];
        end
        
        function xn = getrefstates(obj,ts,frame)
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
            [~,xc] = obj.clock.runat(ts);

            xn = [xs; xc];
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

        function [xm,xn] = getmodelstates(obj,ts,fromnav,cadence)
            %GETMODELSTATES Returns model-derived spacecraft states (PVT).
            %This is the information a user would be able to get from data
            %in the nav message. It encompasses ephemeris and clock model
            %errors when differenced from getrefstates().
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - fromnav; boolean, should data model is interpolating be
            %               propagated from a live nav solution (true) or 
            %               reference info (false)? If false, obj.P0 must be 
            %               populated with default nav uncertainty.
            %    - cadence; # of seconds after which models are updated, can be
            %               separate for ephemeris and clock -- 0 means don't
            %               update
            %   Output:
            %    - xm; model states, pos/vel in J2000
            arguments
                obj     (1,1)   NavSatellite
                ts      (1,:)   double
                fromnav (1,1)   = 0
                cadence (1,2)   double {mustBeNonnegative} = zeros(1,2);
            end

            t0 = ts(1);
            n = length(ts);

            if fromnav
                error("getmodelstates:notImplementedError", ...
                    "Models from nav data has yet to be implemented.");
            else
                xn = obj.getrefstates(ts, 'J2000');
                % add noise to reference trajectory for nav states
                xn = mvnrnd(xn', obj.P0)';
            end

            % find indices for update intervals of models
            propupdate = 1;
            clkupdate = 1;
            if any(cadence)
                t_prop = ts(1);
                t_clk = ts(1);
                for i=2:n
                    if ts(i) > t_prop + cadence(1)
                        propupdate = [propupdate i];
                        t_prop = ts(i);
                    end
                    if ts(i) > t_clk + cadence(2)
                        clkupdate = [clkupdate i];
                        t_clk = ts(i);
                    end
                end
            end
            propupdate = [propupdate n];
            clkupdate = [clkupdate n];

            % generate models and get data
            xm = zeros(9, n);
            for i=1:length(propupdate)-1        % for orbit
                obj.prop.t0 = ts(propupdate(i));            % new start time
                obj.prop.x0 = xn(1:6, propupdate(i));       % new start state
                tf = ts(propupdate(i+1)) - obj.prop.t0;     % end of validity interval
                fs = obj.prop.modelfit("Kepler", tf, 12);   % ephemeris model
                int = propupdate(i):propupdate(i+1);        % validity interval indices
                xm(1:6,int) = fs(ts(int));                  % assign model states
            end
            for i=1:length(clkupdate)-1         % for clock
                obj.clock.t0 = ts(clkupdate(i));            % new start time
                obj.clock.x0 = xn(7:9, clkupdate(i));       % new start state
                fc = obj.clock.modelfit();                  % clock model
                int = clkupdate(i):clkupdate(i+1);          % validity interval indices
                xm(7:9,int) = fc(ts(int));                  % assign model states
            end
        end

        % TODO: implement getsmartmodelstates() where updates are adaptive
    end
end

