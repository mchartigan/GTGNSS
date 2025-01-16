classdef LunarPropagator < OrbitPropagator
    %LUNARPROPAGATOR Generic propagation class for lunar satellites.
    %   Can propagate lunar orbits for various lengths of time and starting
    %   conditions. Mainly created to reduce repetition / verbosity of
    %   scripts.
    
    properties
        % Add any additional properties here
        var (6,6)   double      % continuous-time noise for filters
    end
    
    methods
        function obj = LunarPropagator(t0,x0,ord,nbods,varargin)
            %LUNARPROPAGATOR Construct a LunarPropagator instance.
            %   Inputs:
            %    - t0; character string, 'DD-MMM-YYYY XX:XX:XX'
            %    - x0; starting states -- either array of OE structs, (6,n)
            %          array of starting states (MOON_OP frame), or xopt output
            %          from Conopt(2)
            %    - ord; maximum degree and order of gravity model to use
            %    - nbods; what secondary bodies to include (1: Earth,
            %            2:+sun, 3:+jupiter)
            %    - opts; optional name-value arg, ODE45 integration tolerances
            arguments
                t0      (1,:)
                x0      (:,:)   
                ord     (1,1)   {mustBeInteger,mustBePositive}
                nbods   (1,1)   {mustBeInteger,mustBeNonnegative}
            end
            arguments (Repeating)
                varargin
            end

            % call superclass constructor
            obj = obj@OrbitPropagator(t0, ord, varargin);
            
            cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
            [R,C,S,norms] = cofloader("LP165P.cof", false);
            
            % planetary info
            bods = getplanets('MOON', "MOON", "EARTH", "SUN", "JUPITER");
            bods(1).R = R * 1e-3;           % convert from m to km
            bods(1).C = C;                  % store in moon struct for orbitaldynamics
            bods(1).S = S;                  % store in moon struct for orbitaldynamics
            bods(1).norms = norms;          % store in moon struct for orbitaldynamics
            bods(1).frame = 'MOON_PA';      % body-fixed frame of coefficients
            obj.pri = bods(1);              % primary body
            obj.sec = bods(2:nbods+1);      % secondary bodies

            % estimate noise, Q
            obj.var = diag([[0 0 0] [1 1 1]*3e-5/ord].^2);

            % Parse x0
            errmsg = "x0 must be either an array of OE structs, (6,n) " + ...
                    "array of starting states, or xopt output from Conopt(2)";
            if isa(x0, 'struct')        % oes provided
                oes = x0;
                obj.nsats = length(x0);
            elseif isa(x0, 'double')    % not oes
                % xopt output or states
                if all(size(x0) == [1 14]) || all(size(x0) == [1 24])
                    oes = xopt2oes(x0);
                    obj.nsats = 6;
                elseif size(x0,1) == 6
                    obj.x0 = x0;
                    obj.nsats = size(x0,2);

                    return;     % return early to avoid oes2x0
                else                            % invalid input
                    error("LunarPropagator:invalidInput", errmsg);
                end
            else                                % invalid input
                error("LunarPropagator:invalidInput", errmsg);
            end

            % Convert oes to states
            obj.x0 = zeros(6,obj.nsats);
            for i=1:obj.nsats
                [r,v] = oe2rv(oes(i).a,oes(i).e,oes(i).i,oes(i).RAAN,oes(i).w,oes(i).f,obj.pri.GM);
                obj.x0(:,i) = cspice_sxform('MOON_OP', 'J2000', obj.t0) * [r; v];
            end

            % preallocate information if it was requested
            if ~isempty(obj.t_pre)
                % generate rotation matrices
                T_J2PA = cspice_pxform('J2000', obj.pri.frame, obj.t_pre);

                % convert rotation matrices to quaternions
                l = length(obj.t_pre);
                temp = zeros(l,4);
                for j=1:l
                    temp(j,:) = rotm2quat(T_J2PA(:,:,j));
                end
                quats = quaternion(temp(:,1)', temp(:,2)', temp(:,3)', temp(:,4)');
            end
        end

        function plotlastorbits(obj,frame)
            %PLOTLASTORBITS Generates a plot of the most recently created
            %satellite trajectories in the provided frame.
            %   Input:
            %    - frame; reference frame to plot trajectories in
            arguments
                obj     (1,1)   LunarPropagator
                frame   (1,:)   char
            end

            data = obj.xs;

            % convert data to new frame if required
            if isempty(obj.frame)
                error("plotlastorbits:noData", ...
                    "No data has been generated yet!");
            elseif ~strcmpi(frame, obj.frame)
                for i=1:size(data,2)
                    T = cspice_sxform(obj.frame, frame, obj.ts(i));
                    for j=1:obj.nsats
                        data(:,i,j) = T * data(:,i,j);
                    end
                end
            end

            % % not supported utility
            % plotformat("APA", 1, "coloring", "science");
            plotLunarOrbit(obj.ts, permute(data, [2,1,3]), frame, "Satellite trajectories");
        end

        function [comp,exp] = computedriftrates(obj)
            %COMPUTEDRIFTRATES Computes the drift of right ascension for
            %each orbit over the propagation period.
            %   Output:
            %    - comp; (1,nsats) drift rates computed from frozen orbit eqs
            %    - exp; (1,nsats) drift rates calculated from propagation
            arguments
                obj (1,1)   LunarPropagator
            end
            
            % throw error if there hasn't been a propagation yet
            if isempty(obj.frame)
                error("plotlastorbits:noData", ...
                    "No data has been generated yet!");
            end

            comp = zeros(1,obj.nsats);
            exp = zeros(1,obj.nsats);
            for j=1:obj.nsats
                xo = cspice_sxform(obj.frame, 'MOON_OP', obj.ts(1)) * obj.xs(:,1,j);
                xf = cspice_sxform(obj.frame, 'MOON_OP', obj.ts(end)) * obj.xs(:,end,j);
                [a,e,i,r0,~,~] = rv2oe(xo(1:3), xo(4:6), obj.pri.GM);
                [~,~,~,rf,~,~] = rv2oe(xf(1:3), xf(4:6), obj.pri.GM);

                comp(j) = ascendingnodedrift(a,e,i);
                exp(j) = (rf - r0) / (obj.ts(end) - obj.ts(1));
            end
        end

        function [dRAAN0,dRAANf] = changeinrightascension(obj)
            %CHANGEINRIGHTASCENSION Find the change in relative right ascension
            %spacing over the propagation interval.

            % throw error if there hasn't been a propagation yet
            if isempty(obj.frame)
                error("plotlastorbits:noData", ...
                    "No data has been generated yet!");
            end

            dRAAN0 = zeros(1,obj.nsats);
            dRAANf = zeros(1,obj.nsats);
            RAAN0 = zeros(1,obj.nsats);
            RAANf = zeros(1,obj.nsats);

            % assign starting and ending right ascensions
            for j=1:obj.nsats
                xo = cspice_sxform(obj.frame, 'MOON_OP', obj.ts(1)) * obj.xs(:,1,j);
                xf = cspice_sxform(obj.frame, 'MOON_OP', obj.ts(end)) * obj.xs(:,end,j);
                [~,~,~,RAAN0(j),~,~] = rv2oe(xo(1:3), xo(4:6), obj.pri.GM);
                [~,~,~,RAANf(j),~,~] = rv2oe(xf(1:3), xf(4:6), obj.pri.GM);
            end

            % expand arrays to make finding adjacent spacing easier
            RAAN0 = [RAAN0 RAAN0(1)];
            RAANf = [RAANf RAANf(1)];

            for j=1:obj.nsats
                dRAAN0(j) = angdiff(RAAN0(j+1), RAAN0(j));
                dRAANf(j) = angdiff(RAANf(j+1), RAANf(j));
            end
        end
    end
end

