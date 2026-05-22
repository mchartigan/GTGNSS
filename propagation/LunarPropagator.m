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
        function obj = LunarPropagator(ord,nbods,options)
            %LUNARPROPAGATOR Construct a LunarPropagator instance.
            %   Inputs:
            %    - t0; character string, 'DD-MMM-YYYY XX:XX:XX', or double
            %       of seconds past J2000
            %    - x0; starting states -- either array of OE structs, (6,n)
            %       array of starting states (MOON_OP frame), or xopt output
            %       from Conopt(2)
            %    - ord; maximum degree and order of gravity model to use
            %    - nbods; what secondary bodies to include (1: Earth,
            %       2:+sun, 3:+jupiter)
            %    - opts; optional name-value arg, ODE45 integration tolerances
            arguments
                ord     (1,1)   {mustBeInteger,mustBePositive}
                nbods   (1,1)   {mustBeInteger,mustBeNonnegative}
                options.opts    (1,1)   struct = odeset("RelTol", 1e-9, "AbsTol", 1e-11)
                options.Cr      (1,1)   double {mustBeNonnegative} = 0
                options.Am      (1,1)   double {mustBeNonnegative} = 0
                options.pre     (1,1)   double = 0
                options.units   (1,1)   {mustBeText} = "km"
            end

            % call superclass constructor
            passargs = namedargs2cell(options);
            obj = obj@OrbitPropagator(ord,passargs{:});
            
            % cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
            % [R,C,S,norms] = cofloader("LP165P.cof", false);
            [R,C,S,norms] = sha_loader("gggrx_0900c_sha.tab", 200);
            
            % planetary info
            bods = getplanets('MOON', obj.unit, "MOON", "EARTH", "SUN", "JUPITER");
            bods(1).R = R * 1e-3 * obj.unit;    % convert from m to km if necessary
            bods(1).C = C;                  % store in moon struct for orbitaldynamics
            bods(1).S = S;                  % store in moon struct for orbitaldynamics
            bods(1).norms = norms;          % store in moon struct for orbitaldynamics
            bods(1).frame = 'MOON_PA';      % body-fixed frame of coefficients

            obj.pri = bods(1);              % primary body
            obj.sec = bods(2:nbods+1);      % secondary bodies

            % estimate noise, Q
            obj.var = diag([[0 0 0] [1 1 1]*3e-5/ord].^2);

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
                obj.q_I2F = quaternion(temp(:,1)', temp(:,2)', temp(:,3)', temp(:,4)');

                % change planet positions to piecewise polynomials
                x_pri = obj.pri.x(obj.t_pre);
                pp_pri = spline(obj.t_pre, x_pri);
                obj.pri.x = @(tau) ppval(pp_pri, tau);
                
                for j=1:length(obj.sec)
                    x_secj = obj.sec(j).x(obj.t_pre);
                    pp_secj = spline(obj.t_pre, x_secj);
                    obj.sec(j).x = @(tau) ppval(pp_secj, tau);
                end
            end
        end

        function [comp,exp] = computedriftrates(obj,traj)
            %COMPUTEDRIFTRATES Computes the drift of right ascension for
            %each orbit over the propagation period.
            %   Input:
            %    - traj; Trajectory instance(s)
            %   Output:
            %    - comp; (1,nsats) drift rates computed from frozen orbit eqs
            %    - exp; (1,nsats) drift rates calculated from propagation
            arguments
                obj     (1,1)   LunarPropagator
                traj    (1,:)   Trajectory
            end
            
            nsats = length(traj);
            comp = zeros(1,nsats);
            exp = zeros(1,nsats);

            for j=1:nsats
                xo = traj(j).get(traj.ts(1), outframe='MOON_OP');
                xf = traj(j).get(traj.ts(end), outframe='MOON_OP');
                [a,e,i,r0,~,~] = rv2oe(xo(1:3), xo(4:6), obj.pri.GM);
                [~,~,~,rf,~,~] = rv2oe(xf(1:3), xf(4:6), obj.pri.GM);

                comp(j) = ascendingnodedrift(a,e,i);
                exp(j) = (rf - r0) / (traj(j).ts(end) - traj(j).ts(1));
            end
        end

        function [dRAAN0,dRAANf] = changeinrightascension(obj,traj)
            %CHANGEINRIGHTASCENSION Find the change in relative right ascension
            %spacing over the propagation interval.
            %   Input:
            %    - traj; Trajectory instance(s)
            arguments
                obj     (1,1)   LunarPropagator
                traj    (1,:)   Trajectory
            end

            nsats = length(traj);
            dRAAN0 = zeros(1,nsats);
            dRAANf = zeros(1,nsats);
            RAAN0  = zeros(1,nsats);
            RAANf  = zeros(1,nsats);

            % assign starting and ending right ascensions
            for j=1:nsats
                xo = traj(j).get(traj.ts(1), outframe='MOON_OP');
                xf = traj(j).get(traj.ts(end), outframe='MOON_OP');
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

        function fit = AFSfit(obj,traj,n)
            %MODELFIT Fit the AFS navigation message ephemeris format to the
            %provided trajectory (not well defined atm lol).
            arguments
                obj     (1,1)   LunarPropagator
                traj    (1,1)   Trajectory
                n       (1,1)   {mustBeInteger,mustBeNonnegative}
            end

            diff = n > 0;
            if ~diff, n = 10; end

            % use chebichev nodes for fitting
            t0 = traj.t0;
            span = chebichev(n-1);
            fit.VP = traj.ts(end) - t0;
            teval = (span + 1) * fit.VP / 2 + t0;
            xeval = traj.get(teval, outframe='J2000');
            
            % change to orbital elements
            as = zeros(1,n); es = zeros(1,n); is = zeros(1,n);
            RAANs = zeros(1,n); ws = zeros(1,n); fs = zeros(1,n);
            A_I2ME = zeros(6,n);
            
            for k=1:n
                [as(k),es(k),ik,Ok,wk,fs(k)] = rv2oe(xeval(:,k),obj.pri.GM);
                % % get perifocal to ICRF rotation matrix
                % T_P2I = rotz(-Ok) * rotx(-ik) * rotz(-wk);
                % % get perifocal to MOON ME rotation matrix
                % T_P2ME = cspice_pxform('J2000', 'MOON_ME', teval(k)) * T_P2I;
                % [Ok, ik, wk] = cspice_m2eul(T_P2ME, 3, 1, 3);
                RAANs(k) = Ok;
                is(k) = ik;
                ws(k) = wk;
                T = cspice_sxform('J2000', 'MOON_ME', teval(k));
                A_I2ME(:,k) = cspice_xf2eul(T,3,1,3);
            end
            
            fit.t_oe = traj.t0;
            fit.a = mean(as);
            fit.e = mean(es);
            fit.i = mean(is);
            fit.RAAN = mean(RAANs);
            fit.w = mean(ws);
            fit.M0 = true2mean(fs(1), fit.e);
            fit.A = [A_I2ME(1:3,1)' mean(A_I2ME(4:6,:), 2)'];
            fit.Cx = [];
            fit.Cy = [];
            fit.Cz = [];

            % DIFFERENTIAL CORRECTIONS %
            if diff
                % get effective Keplerian elements
                kepmsg = [t0 0 0 0 t0 0 0 0 0 fit.t_oe ...
                          fit.a fit.e fit.i fit.RAAN fit.w fit.M0 fit.A];
    
                xbase = zeros(6,n);
                for k=1:n
                    temp = RadiometricObsSim.geteph(teval(k), 0, kepmsg, 1, obj.pri.GM);
                    xbase(:,k) = temp(1:6);
                    % T_ME2J = [T_J2ME(1:3,1:3)' zeros(3)
                    %           T_J2ME(4:6,1:3)' T_J2ME(4:6,4:6)'];
                end
                % state error
                dx = xeval - xbase;
        
                % compute coefficients for basis and generate model function
                % % polynomial basis
                % phi = (span').^(0:n-1);
                % chebyshev basis
                phi = chebyshev(0:n-1, span);
                G = pinv(phi) * dx(1:3,:)';
                fit.Cx = G(:,1)';
                fit.Cy = G(:,2)';
                fit.Cz = G(:,3)';
    
                % how to compute values: (basis(2*(tau)/dt - 1) * H)
            end

            
        end

        function plot(obj,traj,frame)
            %PLOT Generates a plot of the provided satellite trajectories in
            %the specified frame.
            %   Input:
            %    - traj; Trajectory instance(s)
            %    - frame; reference frame to plot trajectories in
            arguments
                obj     (1,1)   LunarPropagator
                traj    (1,:)   Trajectory
                frame   (1,:)   char
            end

            ts = traj(1).ts;
            nsats = length(traj);
            data = zeros(length(ts),3,nsats);

            % convert data to new frame if required
            for i=1:nsats
                x = traj(i).get(ts, outframe=frame);
                data(:,:,i) = x(1:3,:)';
            end

            % % not supported utility
            figure();
            plotformat("APA", 1);
            % Display moon in trajectory plot
            R = obj.pri.R;
            [Imoon, ~] = imread("Moon_HermesCelestiaMotherlode.jpg");
            [xx, yy, zz] = ellipsoid(0, 0, 0, R, R, R);
            
            % Rotate moon from MOON_ME frame to plot_frame
            T = cspice_pxform('MOON_ME', frame, ts(end));
            for j=1:size(xx,1)
                for k=1:size(xx,2)
                    % -z to flip image
                    tmp = T * [xx(j,k); yy(j,k); -zz(j,k)];
                    xx(j,k) = tmp(1); yy(j,k) = tmp(2); zz(j,k) = tmp(3);
                end
            end
            
            globe = surf(xx, yy, zz);
            set(globe, 'FaceColor', 'texturemap', 'CData', Imoon, 'FaceAlpha', 1, ...
                'EdgeColor', 'none');
            hold on;
            
            % put a little star on the south pole
            lsp = T * [0;0;-R];
            scatter3(lsp(1), lsp(2), lsp(3), 100, "red", "filled", "pentagram");
            
            % plot user trajectory for same time frame
            styles = {'-', '--', '-.', ':'};
            for k=1:size(data,3)
                plot3(data(:,1,k), data(:,2,k), data(:,3,k), "LineWidth", 1.5, ...
                    "Marker", "diamond", "MarkerIndices", length(ts), ...
                    LineStyle=styles{mod(k-1,4)+1});
            end
            
            grid on; axis equal;
            if strcmp(frame, 'J2000'), frame = 'ICRF'; end
            SUB = strsplit(frame,"_");
            SUB = SUB(end);

            % make sure units are appropriately tagged in plots
            if obj.unit == 1e3
                UNIT = "m";
            elseif obj.unit == 1
                UNIT = "km";
            else
                error("plot:unknownUnit", ...
                    "Unit multiplier %d not known.", obj.unit);
            end
            xlabel("x_{"+SUB+"} ("+UNIT+")");
            ylabel("y_{"+SUB+"} ("+UNIT+")");
            zlabel("z_{"+SUB+"} ("+UNIT+")");
            title("Satellite trajectories");
        end
    end
end

