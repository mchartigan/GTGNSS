function var = getclockjitter(clock,Bn)
%GETCLOCKJITTER Returns the jitter noise of a clock at a specific noise
%bandwidth, based on the phase noise statistics provided in the datasheets.
%   Inputs:
%    - clock; string of clock type, must be from list: "MicrochipCSAC",
%             "MicrochipMAC", "SafranMAC", "SafranMiniRAFS",
%             "RakonMiniUSO", "AccubeatUSO", "SafranRAFS", "ExcelitasRAFS"
%    - Bn; carrier loop noise bandwidth
%   Outputs:
%    - var; variance of clock jitter, radians^2
%
%   Ref: Zucca, C. and Tavella, P.; doi.org/10.1109/TUFFC.2005.1406554
arguments
    clock   (1,:)   string
    Bn      (1,1)   double {mustBePositive}
end

c = 299792458;                          % m/s, speed of light

switch clock
    case "MicrochipCSAC"
        0

    case "MicrochipMAC"
        freq = [1e0 1e1 1e2 1e3 1e4];       % frequency offset, Hz
        noise = [-70 -90 -114 -135 -140];   % phase noise, dBc/Hz

    case "SafranMAC"
        0

    case "SafranMiniRAFS"
        0

    case "RakonMiniUSO"
        0

    case "AccubeatUSO"
        0

    case "SafranRAFS"
        0

    case "ExcelitasRAFS"
        0

    otherwise
        error("getclockjitter:invalidClock", ...
            "%s not found, valid clock options found in " + ...
            "documentation. Strings must match exactly.", clock);
end

Bn = Bn / 2;        % noise bandwidth presumed two-sided, so get one side
noise = 10.^(noise/10);
n_Bn = interp1(freq, noise, Bn);
ii = find(freq < Bn);
f_int = [freq(ii) Bn];
n_int = [noise(ii) n_Bn];
A = trapz(f_int, n_int);
A = 10*log10(A);

var = 2*10^(A/10);
end

