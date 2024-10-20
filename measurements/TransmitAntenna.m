classdef TransmitAntenna < Antenna
    %TRANSMITANTENNA Properties specific to transmitting antenna
    
    properties
        % dBW, transmit power of antenna (default 25W)
        P       (1,1)   double = 10*log10(25)
        % Hz, transmit center frequency (default 2492.028 MHz, LunaNet AFS)
        freq    (1,1)   double {mustBePositive} = 2492.028e6
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
