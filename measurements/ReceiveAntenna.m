classdef ReceiveAntenna < Antenna
    %RECEIVEANTENNA Properties specific to receiving antenna
    
    properties
        % dB, system losses in receiver (default no losses)
        As      (1,1)   double = 0
        % K, noise temperature of antenna [default 150K, ballpark temp of
        % moon (100-400K) + pointing at space]
        Ts      (1,1)   double = 150
        % dB, noise figure of receiver [default -3 dB from ODTBX gpsmeas()]
        Nf      (1,1)   double = -3
        % dB, receiver conversion losses [default -1.5 dB from ODTBX gpsmeas()]
        L       (1,1)   double = -1.5
    end
    properties (Constant)
        k = 1.3803e-23      % J/K, Boltzmann's constant
    end

    methods
        function CN0 = getCN0(obj,AP,beta)
            %GETGAIN Computes the carrier to noise density ratio at
            %the user receiver based on received power and antenna
            %parameters.
            %   Input:
            %    - AP; power at the user antenna, dBW
            %    - beta; angle between direction to central body of satellite
            %       and to satellite, rad. A lunar satellite would be
            %       Earth-pointing for GNSS and nadir pointing for LunaNet
            %   Output:
            %    - CN0; carrier-to-noise density ratio, dB-Hz
            arguments (Input)
                obj     (1,1)   ReceiveAntenna
                AP      (1,:)   double
                beta    (1,:)   double
            end

            if numel(obj.gain) > 1
                % convert to deg, wrap to [0,180] (ignore directionality)
                beta = abs(wrapToPi(beta)) * 180/pi;
                G = interp1(obj.gain(:,1), obj.gain(:,2), beta);
            else
                G = ones(size(beta)) * obj.gain;
            end
            % if outside of defined antenna pattern, -300 dB gain
            G(isnan(G)) = -300;

            RP = AP + G + obj.As;           % dBW, gain before amps
            N0 = 10*log10(obj.k * obj.Ts);  % dBW/Hz, noise power spectral density
            CN0 = RP + obj.Nf + obj.L - N0; % dB*Hz, carrier to noise density ratio
        end
    end
end

% NOTES %
%   P_sv    double  1       spacecraft transmit power [dBW]
%   Ts      double  1       System noise temp [K]
%   Ae      double  1       attenuation due to atmosphere (should be negative) [dB]
%   Nf      double  1       dB, Noise figure of receiver/LNA
%   L       double  1       Receiver implementation, A/D conversion losses [dB]
%   Ar      double  Nx1     receive antenna gain (dBi)
%   At      double  Nx1     transmit antenna gain (dBi)
%   As      double  1       dB, System losses, in front of LNA
%   Ad      double  Nx1     Attenuation from R^2 losses (dB)
%   AP      double  Nx1     budget gain before receiver antenna (dBW)
%   RP      double  Nx1     budget gain before receiver amplifiers and conversion (dBW)
% The following are the link budget equations:
%   Ad  = 20.*log10((C/freq)./(4*pi*los_mag));
%   AP  = P_sv + At + Ad + Ae;
%   RP  = AP + Ar + As;
%   CN0 = RP + Nf + L - (10*log10(k*Ts));
