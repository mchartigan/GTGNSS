classdef Receiver < handle
    %RECEIVER Defines all the attributes of a RNSS receiver, including code
    %and carrier tracking loops.
    
    properties
        % speed of light (m/s)
        c           (1,1)   double = 299792458

        % signal characteristics %
        % Hz, receiving center frequency (default 2492.028 MHz, LunaNet AFS)
        freq        (1,1)   double {mustBePositive} = 2492.028e6
        % Hz (or chips/s), (default AFS Q-channel BPSK(5) speading code rate)
        Rc          (1,1)   double {mustBePositive} = 5.115e6

        % code tracking loop characteristics %
        % code tracking loop order
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
        carrierorder(1,1)   {mustBePositive,mustBeInteger,mustBeLessThan(carrierorder,4)} = 2
        % 1 - high C/N0, 2 - near threshold (default high)
        F           (1,1)   {mustBePositive,mustBeInteger,mustBeLessThan(F,3)} = 1
        % Hz, carrier loop noise bandwidth (default moderate)
        Bn_c        (1,1)   double {mustBePositive} = 2
        % s, carrier predetection integration time; must be half of data
        % bit transition time so it can bet two samples to form the
        % discriminator; doesn't matter if data = 0 (default AFS-complaint)
        T_c         (1,1)   double {mustBePositive} = 0.002
        % s, delta pseudorange measurement spacing to derive rate (default 1s)
        Tm          (1,1)   double {mustBePositive} = 1
    end

    properties (Access = private)
        % code tracking loop type, can only be DLL
        code        (1,1)   string = "DLL"
    end

    methods
        function obj = Receiver(carrierloop)
            %RECEIVER Creates a Receiver instance. Specify at least the
            %carrier tracking loop type (or "none").
            %   Input:
            %    - carrierloop; carrier tracking loop type, "PLL", "FLL",
            %                   or "none"

            obj.carrier = carrierloop;
        end

        function [err,var] = noise(obj,CN0,ts,r,dr)
            %NOISE Returns the range and range-rate error (and variance) of
            %the receiver measurements.
            %   Input:
            %    - CN0; signal to noise density ratio (dB-Hz)
            %    - ts; eval time steps, seconds past J2000
            %    - r; transmitter-receiver ranges (m) to compute error for
            %    - dr; transmitter-receiver range-rates (mm/s)
            %   Output:
            %    - err; random zero-mean variables with variance vrec,
            %           first row is range (m) and second is range-rate (mm/s)
            %    - var; variance of range (m^2) and range-rate
            %           measurements (mm^2/s^2)
            arguments
                obj     (1,1)   Receiver
                CN0     (1,:)   double
                ts      (1,:)   double
                r       (1,:)   double {mustBePositive}
                dr      (1,:)   double
            end


            CN0 = 10.^(CN0/10);                 % Hz, converted from dB-Hz for below equations
            lambda = obj.c / obj.freq * 1e3;    % mm, wavelength of carrier
            Tc = 1/obj.Rc;                      % s (or s/chip), chip period

            % compute additional line-of-sight dynamics
            ddr = gradient(dr, ts);             % mm/s^2, acceleration
            dddr = gradient(ddr, ts);           % mm/s^3, jerk

            % assign variance based on code tracking loop design
            if strcmpi(obj.code, "DLL")         % delay lock loop
                % assign thermal noise based on E-L correlator spacing
                if obj.D >= pi*obj.Rc/obj.Bfe
                    var = obj.Bn./(2*CN0) .* obj.D .* ...
                          (1 + obj.data*2./(obj.T*CN0*(2-obj.D)));
                elseif obj.D > obj.Rc/obj.Bfe
                    var = obj.Bn./(2*CN0) .* ...
                          (1/(obj.Bfe*Tc) + obj.Bfe*Tc/(pi-1)*(obj.D - (1/(obj.Bfe*Tc)))^2) .* ...
                          (1 + obj.data*2./(obj.T*CN0*(2-obj.D)));
                else
                    var = obj.Bn./(2*CN0) .* (1/(obj.Bfe*Tc)) .* ...
                          (1 + obj.data*1./(obj.T*CN0));
                end

                % add dynamic stress error if not carrier-aided
                if strcmpi(obj.carrier, "none")
                    switch obj.codeorder
                        case 1
                            w0 = 4 * obj.Bn;
                            dyn = dr;
                        case 2
                            w0 = obj.Bn / 0.53;
                            dyn = ddr;
                        case 3
                            w0 = obj.Bn / 0.7845;
                            dyn = dddr;
                        otherwise
                            error("noise:invalidCodeOrder", ...
                                "Supported code tracking loop orders are: 1, 2, 3.");
                    end

                    % DLL noise is thermal noise + dynamic stress error
                    var = (sqrt(var) + abs(dyn)*1e-3 / (3*w0^obj.codeorder)).^2;
                end

                % convert from chips^2 to meters^2
                var = var * (obj.c * Tc)^2;
            else
                error("noise:invalidCodeLoop", ...
                    "Supported code tracking loops are: 'DLL'.");
            end

            % assign variance based on carrier tracking loop design
            if strcmpi(obj.carrier, "PLL")      % phase lock loop
                % thermal noise
                vdop = (lambda/(2*pi))^2 * (obj.Bn_c./CN0 .* ...
                       (1 + obj.data*1./(obj.T_c*CN0)));

                % add dynamic stress
                switch obj.carrierorder
                    case 1
                        w0 = 4 * obj.Bn_c;
                        dyn = dr;
                    case 2
                        w0 = obj.Bn_c / 0.53;
                        dyn = ddr;
                    case 3
                        w0 = obj.Bn_c / 0.7845;
                        dyn = dddr;
                    otherwise
                        error("noise:invalidCarrierOrder", ...
                            "Supported PLL loop orders are: 1, 2, 3.");
                end

                % add dynamic stress error
                vdop = (sqrt(vdop) + abs(dyn) / (3*w0^obj.carrierorder)).^2;
                % vdop in mm^2, multiply by 2/Tm^2 to convert to mm^2/s^2
                % (difference of two random variables)
                var = [var; vdop * 2/obj.Tm^2];
                
            elseif strcmpi(obj.carrier, "FLL")  % frequency lock loop
                % thermal noise
                vdop = (lambda/(2*pi*obj.T_c))^2 * ...
                       (4*obj.F*obj.Bn_c./CN0 .* (1 + obj.data*1./(obj.T_c*CN0)));

                % add dynamic stress
                switch obj.carrierorder
                    case 1
                        w0 = 4 * obj.Bn_c;
                        dyn = ddr;
                    case 2
                        w0 = obj.Bn_c / 0.53;
                        dyn = dddr;
                    otherwise
                        error("noise:invalidCarrierOrder", ...
                            "Supported FLL loop orders are: 1, 2.");
                end

                % FLL has one more integrator, so dynamic stress is proportional
                % to d^(n+1)R/dt^(n+1)
                vdop = (sqrt(vdop) + abs(dyn) / (3*w0^obj.carrierorder)).^2;
                var = [var; vdop];
            elseif ~strcmpi(obj.carrier, "none")
                % loop isn't none (no carrier tracking)
                error("linkbudget:invalidCarrierLoop", ...
                    "Supported carrier tracking loops are: 'PLL', 'FLL', 'none'.");
            end

            % generate noise based on var
            m = size(var,1);        % range or range and Doppler?
            n = length(r);          % number of data points
            err = zeros(m,n);
            for i=1:n
                err(:,i) = mvnrnd(zeros(1,m), diag(var(:,i)))';
            end
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

