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
