classdef EarthPropagator < OrbitPropagator
    %EARTHPROPAGATOR Generic propagation class for Earth-centered satellites.
    %   Can propagate orbits for various lengths of time and starting
    %   conditions. Mainly created to reduce repetition / verbosity of
    %   scripts.
    
    properties
        % Add any additional properties here
    end
    
    methods
        function obj = EarthPropagator(ord,nbods,options)
            %EARTHPROPAGATOR Construct an EarthPropagator instance.
            %   Inputs:
            %    - t0; character string, 'DD-MMM-YYYY XX:XX:XX'
            %    - x0; starting states -- either array of OE structs, (6,n)
            %          or array of starting states (J2000 frame)
            %    - ord; maximum degree and order of gravity model to use
            %    - nbods; what secondary bodies to include (1:+moon,
            %            2:+sun, 3:+jupiter)
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
            [R,C,S,norms] = cofloader("JGM3.cof", false);
            
            % planetary info
            bods = getplanets('EARTH', obj.unit, "EARTH", "MOON", "SUN", "JUPITER");
            bods(1).R = R * 1e-3 * obj.unit;           % convert from m to km
            bods(1).C = C;                  % store in earth struct for orbitaldynamics
            bods(1).S = S;                  % store in earth struct for orbitaldynamics
            bods(1).norms = norms;          % store in earth struct for orbitaldynamics
            bods(1).frame = 'ITRF93';       % body-fixed frame of coefficients
            obj.pri = bods(1);              % primary body
            obj.sec = bods(2:nbods+1);      % secondary bodies
        end

        function [comp,exp] = computedriftrates(obj)
            %COMPUTEDRIFTRATES Computes the drift of right ascension for
            %each orbit over the propagation period.
            %   Output:
            %    - comp; (1,nsats) drift rates computed from frozen orbit eqs
            %    - exp; (1,nsats) drift rates calculated from propagation
            arguments
                obj (1,1)   EarthPropagator
            end
            
            % throw error if there hasn't been a propagation yet
            if isempty(obj.frame)
                error("plotlastorbits:noData", ...
                    "No data has been generated yet!");
            end

            comp = zeros(1,obj.nsats);
            exp = zeros(1,obj.nsats);
            for j=1:obj.nsats
                xo = cspice_sxform(obj.frame, 'J2000', obj.ts(1)) * obj.xs(:,1,j);
                xf = cspice_sxform(obj.frame, 'J2000', obj.ts(end)) * obj.xs(:,end,j);
                [a,e,i,r0,~,~] = rv2oe(xo(1:3), xo(4:6), obj.pri.GM);
                [~,~,~,rf,~,~] = rv2oe(xf(1:3), xf(4:6), obj.pri.GM);

                comp(j) = ascendingnodedrift(a,e,i);
                exp(j) = (rf - r0) / (obj.ts(end) - obj.ts(1));
            end
        end

        function plot(obj,traj,frame)
            %PLOT Generates a plot of the provided satellite trajectories in
            %the specified frame.
            %   Input:
            %    - traj; Trajectory instance(s)
            %    - frame; reference frame to plot trajectories in
            arguments
                obj     (1,1)   EarthPropagator
                traj    (1,:)   Trajectory
                frame   (1,:)   char
            end

            ts = traj(1).ts;
            nsats = length(traj);
            data = zeros(length(ts),3,nsats);

            % convert data to new frame if required
            for i=1:nsats
                x = traj(i).get(ts, frame);
                data(:,:,i) = x(1:3,:)';
            end

            figure();
            plotformat("APA", 1);
            % Display Earth in trajectory plot
            R = obj.pri.R;
            [Iearth, ~] = imread("ModifiedBlueMarble.jpg");
            [xx, yy, zz] = ellipsoid(0, 0, 0, R, R, R);
            
            % Rotate Earth from IAU frame to plot_frame
            T = cspice_pxform('IAU_EARTH', frame, ts(end));
            for j=1:size(xx,1)
                for k=1:size(xx,2)
                    % -z to flip image
                    tmp = T * [xx(j,k); yy(j,k); -zz(j,k)];
                    xx(j,k) = tmp(1); yy(j,k) = tmp(2); zz(j,k) = tmp(3);
                end
            end
            
            globe = surf(xx, yy, zz);
            set(globe, 'FaceColor', 'texturemap', 'CData', Iearth, 'FaceAlpha', 1, ...
                'EdgeColor', 'none');
            hold on;
            
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
            xlabel("x_{"+SUB+"} (km)");
            ylabel("y_{"+SUB+"} (km)");
            zlabel("z_{"+SUB+"} (km)");
            title("Satellite trajectories");
        end

        function fit = AFSfit(obj,traj,n)
            %MODELFIT Fit the AFS navigation message ephemeris format to the
            %provided trajectory (not well defined atm lol).
            arguments
                obj     (1,1)   EarthPropagator
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
            xeval = traj.get(teval, 'J2000');
            
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
                kepmsg = [t0 0 0 0 0 0 0 0 0 fit.t_oe ...
                          fit.a fit.e fit.i fit.RAAN fit.w fit.M0 fit.A];
    
                xbase = zeros(6,n);
                for k=1:n
                    temp = RadiometricObsSim.geteph(teval(k), 0, ...
                        kepmsg, 1, obj.pri.GM);
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
    end
end

