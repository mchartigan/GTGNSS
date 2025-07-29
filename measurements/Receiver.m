classdef Receiver < handle
    %RECEIVER Defines all the attributes of a RNSS receiver, including code
    %and carrier tracking loops.
    
    properties
        % antenna object
        ant     (1,1)   ReceiveAntenna
        % receiver clock, units in m %
        clock       (1,1)   Clock = Clock("none", zeros(4,1))

        % signal characteristics %
        % Hz, receiving center frequency (default 2492.028 MHz, LunaNet AFS)
        freq        (1,1)   double {mustBePositive} = 2492.028e6
        % Hz (or chips/s), (default AFS Q-channel BPSK(5) speading code rate)
        Rc          (1,1)   double {mustBePositive} = 5.115e6

        % code tracking loop characteristics %
        % code tracking loop order (first-, second-, or third-order)
        codeorder   (1,1)   {mustBePositive,mustBeInteger,mustBeLessThan(codeorder,4)} = 2
        % is data present on signal? (i.e. does it have a nav msg or is it a 
        % pilot channel?) 1 - yes, 0 - no
        data        (1,1)   {mustBeNonnegative,mustBeInteger,mustBeLessThan(data,2)} = 1
        % Hz, front-end bandwidth of receiver
        % bandwidth of AFS is 2500-2483.5=16.5 MHz, so Bfe >= 16.5 MHz
        % or, Bfe >= 2*Rc = 10.23 MHz
        Bfe         (1,1)   double {mustBePositive} = 10.23e6
        % Hz, DLL noise bandwidth (default 3 Hz for unaided code tracking)
        Bn          (1,1)   double {mustBePositive} = 3
        % s, code predetection integration time; w/o a pilot channel, T can only be 
        % as long as the nav message bit transitions (20ms for GPS or 1-4ms for AFS)
        T           (1,1)   double {mustBePositive} = 0.004
        % chips, early-late correlator spacing; no use making < Rc/Bfe (default that)
        D           (1,1)   double {mustBePositive} = 0.5

        % carrier tracking loop characteristics %
        % carrier tracking loop type and order
        carrier     (1,1)   string {mustBeCarrierLoop} = "none"
        % second- or third-order for PLL, first- or second-order for FLL
        carrierorder(1,1)   {mustBePositive,mustBeInteger,mustBeLessThan(carrierorder,4)} = 2
        % 1 - high C/N0, 2 - near threshold (default high)
        F           (1,1)   {mustBePositive,mustBeInteger,mustBeLessThan(F,3)} = 1
        % Hz, carrier loop noise bandwidth (default moderate)
        Bn_c        (1,1)   double {mustBePositive} = 2
        % s, carrier predetection integration time; must be half of data
        % bit transition time so it can bet two samples to form the
        % discriminator; doesn't matter if data = 0 (default AFS-complaint)
        T_c         (1,1)   double {mustBePositive} = 0.002
        % s, phase measurement spacing to derive rate (default 0s). If Tm <
        % T_c, then T_c is used to derive rate (like in the case of a FLL).
        % If Tm > T_c, Tm is used; this is mainly for PLLs when you want a
        % more accurate measurement because you have the full phase counts
        % for that whole time. Similar techniques may be implementable for
        % a FLL?
        Tm          (1,1)   double {mustBeNonnegative} = 0
    end

    properties (Access = private)
        % code tracking loop type, can only be DLL
        code        (1,1)   string = "DLL"
    end

    properties (Constant)
        % speed of light (m/s)
        c           (1,1)   double = 299792458
    end

    methods
        function obj = Receiver(antenna,clock,carrierloop)
            %RECEIVER Creates a Receiver instance. Specify at least the
            %carrier tracking loop type (or "none").
            %   Input:
            %    - clock; receiver clock
            %    - carrierloop; carrier tracking loop type, "PLL", "FLL",
            %                   or "none"

            if nargin ~= 0
                obj.ant = antenna;
                obj.clock = clock;
                obj.carrier = carrierloop;
            end
        end

        function [y,err,var] = tracksat(obj,ts,T,dT,AP)
            %NOISE Returns the range and range-rate error (and variance) of
            %the receiver measurements.
            %   Input:
            %    - ts; eval time steps, seconds past J2000
            %    - T; transmitter-receiver delays (s) to compute error for
            %    - dT; transmitter-receiver Doppler (s/s)
            %    - AP; power at the user antenna, dBW
            %   Output:
            %    - y; returned measurements of signal. For DLL, y(1,:) is
            %         the measured transmission delay in s. If PLL, y(2,:) is
            %         the carrier phase (w/ integer ambiguity) in s. If FLL,
            %         y(3,:) is the Doppler shift in s/s.
            %    - err; random zero-mean variables with variance vrec,
            %           first row is range (m) and second is range-rate (m/s)
            %    - var; variance of range (m^2) and range-rate
            %           measurements (m^2/s^2)
            arguments
                obj         (1,1)   Receiver
                ts          (1,:)   double
                T           (1,:)   double
                dT          (1,:)   double
                AP          (1,:)   double
            end

            CN0 = obj.rxlinkbudget(AP);
            CN0 = 10.^(CN0/10);                 % Hz, converted from dB-Hz for below equations
            Tc = 1/obj.Rc;                      % s (or s/chip), chip period

            % compute additional line-of-sight dynamics
            ddT = gradient(dT, ts);             % s/s^2, Doppler rate (accel)
            dddT = gradient(ddT, ts);           % s/s^3, Doppler acceleration (jerk)

            % receiver tracking constraints
            n = length(T);                      % number of data points

            % create output measurement vector
            % [code offset (s); carrier offset (s); Doppler (s/s)]
            y = [T; T; dT];
            % make carrier offset from first valid measurement in
            % contiguous series of valid measurements
            phi0 = y(2,1);
            for i=1:n-1
                if isnan(y(2,i)) && ~isnan(y(2,i+1))
                    phi0 = y(i+1);
                elseif ~isnan(y(2,i))
                    y(2,i) = y(2,i) - phi0;
                end
            end
            y(2,end) = y(2,end) - phi0;

            % initialize variances and error
            err = nan(3,n);
            track = false(3,n);
            var.thermal = zeros(3,n);
            var.clk = zeros(3,n);
            var.dyn = zeros(3,n);
            var.total = zeros(3,n);

            % assign variance based on code tracking loop design
            if strcmpi(obj.code, "DLL")         % delay lock loop
                % assign thermal noise based on E-L correlator spacing
                if obj.D >= pi*obj.Rc/obj.Bfe
                    var.thermal(1,:) = obj.Bn./(2*CN0) .* obj.D .* ...
                          (1 + obj.data*2./(obj.T*CN0*(2-obj.D)));
                elseif obj.D > obj.Rc/obj.Bfe
                    var.thermal(1,:) = obj.Bn./(2*CN0) .* ...
                          (1/(obj.Bfe*Tc) + obj.Bfe*Tc/(pi-1)*(obj.D - (1/(obj.Bfe*Tc)))^2) .* ...
                          (1 + obj.data*2./(obj.T*CN0*(2-obj.D)));
                else
                    var.thermal(1,:) = obj.Bn./(2*CN0) .* (1/(obj.Bfe*Tc)) .* ...
                          (1 + obj.data*1./(obj.T*CN0));
                end

                % uncertainty from clock phase noise is zero
                % add dynamic stress error if not carrier-aided
                if strcmpi(obj.carrier, "none")
                    switch obj.codeorder
                        case 1
                            w0 = 4 * obj.Bn;
                            dyn = dT * obj.Rc;      % chips/s
                        case 2
                            w0 = obj.Bn / 0.53;
                            dyn = ddT * obj.Rc;     % chips/s^2            
                        case 3
                            w0 = obj.Bn / 0.7845;
                            dyn = dddT * obj.Rc;    % chips/s^3
                        otherwise
                            error("noise:invalidCodeOrder", ...
                                "Supported code tracking loop orders are: 1, 2, 3.");
                    end

                    % DLL noise is thermal noise + dynamic stress, in chips^2 
                    var.dyn(1,:) = (abs(dyn) / (3*w0^obj.codeorder)).^2;            
                end

                var.total(1,:) = var.thermal(1,:) + var.clk(1,:) + var.dyn(1,:);

                % compute validity
                track(1,:) = 3*sqrt(var.total(1,:)) <= obj.D / 2;
                
                % convert from chips^2 to s^2
                var.thermal(1,:) = var.thermal(1,:) * Tc^2;
                var.clk(1,:) = var.clk(1,:) * Tc^2;
                var.dyn(1,:) = var.dyn(1,:) * Tc^2;
                var.total(1,:) = var.total(1,:) * Tc^2;

            else
                error("noise:invalidCodeLoop", ...
                    "Supported code tracking loops are: 'DLL'.");
            end

            % assign variance based on carrier tracking loop design
            if strcmpi(obj.carrier, "PLL")      % phase lock loop
                % thermal noise, rad
                var.thermal(2,:) = obj.Bn_c./CN0 .* (1 + obj.data*1./(obj.T_c*CN0));

                % add dynamic stress
                switch obj.carrierorder
                    case 2
                        w0 = obj.Bn_c / 0.53;
                        dyn = ddT*obj.freq*2*pi;    % rad/s^2
                        coef = 2.5;

                    case 3
                        w0 = obj.Bn_c / 0.7845;
                        dyn = dddT*obj.freq*2*pi;   % rad/s^3
                        coef = 2.25;
                        
                    otherwise
                        error("noise:invalidCarrierOrder", ...
                            "Supported PLL loop orders are: 2, 3.");
                end

                % add oscillator phase noise, rad^2
                var.clk(2,:) = (2*pi/coef/obj.Bn_c * obj.freq)^2 * ...
                           obj.clock.stability(1/obj.Bn_c) * ones(1,size(var.thermal,2));
                % add dynamic stress error
                var.dyn(2,:) = (abs(dyn) / (3*w0^obj.carrierorder)).^2;
                var.total(2,:) = var.thermal(2,:) + var.clk(2,:) + var.dyn(2,:);

                % compute validity (assume ATAN2 discriminator)
                track(2,:) = 3*sqrt(var.total(2,:)) <= 2*pi / (4*(1+obj.data));

                % convert from rad^2 to s^2
                var.thermal(2,:) = var.thermal(2,:) * (2*pi*obj.freq)^(-2);
                var.clk(2,:) = var.clk(2,:) * (2*pi*obj.freq)^(-2);
                var.dyn(2,:) = var.dyn(2,:) * (2*pi*obj.freq)^(-2);
                var.total(2,:) = var.total(2,:) * (2*pi*obj.freq)^(-2);
                
            elseif strcmpi(obj.carrier, "FLL")  % frequency lock loop
                % thermal noise, Hz
                var.thermal(3,:) = 1/(2*pi*obj.T_c)^2 * ...
                       (4*obj.F*obj.Bn_c./CN0 .* (1 + obj.data./(obj.T_c*CN0)));

                % add dynamic stress, Hz
                switch obj.carrierorder
                    case 1
                        w0 = 4 * obj.Bn_c;
                        dyn = ddT*obj.freq;     % cycles/s^2
                    case 2
                        w0 = obj.Bn_c / 0.53;
                        dyn = dddT*obj.freq;    % cycles/s^3
                    otherwise
                        error("noise:invalidCarrierOrder", ...
                            "Supported FLL loop orders are: 1, 2.");
                end

                % uncertainty from clock phase noise is zero
                % FLL has one more integrator, so dynamic stress is proportional
                % to d^(n+1)R/dt^(n+1)
                var.dyn(3,:) = (abs(dyn) / (3*w0^obj.carrierorder)).^2;
                var.total(3,:) = var.thermal(3,:) + var.clk(3,:) + var.dyn(3,:);

                % compute validity (assume ATAN2 discriminator)
                track(3,:) = 3*sqrt(var.total(3,:)) <= 1 / (4*obj.T_c);

                % convert from Hz^2 to (s/s)^2
                var.thermal(3,:) = var.thermal(3,:) * obj.freq^(-2);
                var.clk(3,:) = var.clk(3,:) * obj.freq^(-2);
                var.dyn(3,:) = var.dyn(3,:) * obj.freq^(-2);
                var.total(3,:) = var.total(3,:) * obj.freq^(-2);

            elseif ~strcmpi(obj.carrier, "none")
                % loop isn't none (no carrier tracking)
                error("linkbudget:invalidCarrierLoop", ...
                    "Supported carrier tracking loops are: 'PLL', 'FLL', 'none'.");
            end

            % generate noise based on var
            for i=1:n
                ind = ~isnan(var.total(:,i));
                err(ind,i) = mvnrnd(zeros(1,sum(ind)), diag(var.total(ind,i)))';
            end
            % apply noise to measurements
            y = y + err;
            % mask out invalid measurements
            y(~track) = NaN;
        end

        function CN0 = rxlinkbudget(obj,AP)
            %RXLINKBUDGET Computes the carrier to noise density ratio at
            %the user receiver based on received power and antenna
            %parameters.
            %   Input:
            %    - AP; power at the user antenna, dBW
            %   Output:
            %    - CN0; carrier-to-noise density ratio, dB-Hz
            arguments (Input)
                obj (1,1)   Receiver
                AP  (1,:)   double
            end
            
            RP = AP + obj.ant.gain + obj.ant.As;        % dBW, gain before amps
            k = 1.3803e-23;                             % J/K, Boltzmann's constant
            N0 = 10*log10(k * obj.ant.Ts);              % dBW/Hz, noise power spectral density
            CN0 = RP + obj.ant.Nf + obj.ant.L - N0;     % dB*Hz, carrier to noise density ratio
        end
    end
end

% custom validation
function mustBeCarrierLoop(txt)
%MUSTBEOCARRIERLOOP Tests that string option is a valid carrier tracking
%loop type. Options are "PLL", "FLL", or "none".

if ~sum(strcmpi(txt,["PLL","FLL","none"]))
    error("carrier:notOption", "Carrier tracking loop type must be 'PLL', 'FLL', or 'none'.");
end
end

