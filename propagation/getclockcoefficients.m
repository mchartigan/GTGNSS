function [coeff,a,dev] = getclockcoefficients(clock,fit,plot)
%GETCLOCKCOEFFICIENTS Returns the diffusion coefficients and aging rate of
%a clock, specified by name.
%   Diffusion coefficients are found by fitting short-term stability data
%   from product datasheets to classical Allan or Hadamard variance models,
%   described in (Zucca and Tavella, 2005).
%
%   Inputs:
%    - clock; string of clock type, must be from list: "MicrochipCSAC",
%             "MicrochipMAC", "SafranMAC", "SafranMiniRAFS",
%             "RakonMiniUSO", "AccubeatUSO", "SafranRAFS", "ExcelitasRAFS"
%    - fit; optional, model to fit. Either "Allan" (default) or "Hadamard"
%    - plot; optional, boolean (default false) to plot model fit
%   Outputs:
%    - coeff; row vector of diffusion coefficients
%    - a; aging rate of clock, in s/s/s
%    - dev; 3-sigma deviation @10s (mm/s)
%
%   Ref: Zucca, C. and Tavella, P.; doi.org/10.1109/TUFFC.2005.1406554
arguments
    clock   (1,:)   string
    fit     (1,:)   string = "Allan"
    plot    (1,1)   double {mustBeNonnegative} = false
end

c = 299792458;                          % m/s, speed of light

switch clock
    case "MicrochipCSAC"
        a = 9e-10 / (30 * 86400);           % s/s/s, aging rate a
        taus = [1 10 100 1000]';            % intervals for short-term stability
        stds = [3e-10 1e-10 3e-11 1e-11]';  % deviations for taus
        % 3-sigma deviation @ 10s (mm/s)
        dev = 1e3 * c * 3 * stds(2);

    case "MicrochipMAC"
        a = 5e-11 / (30 * 86400);           % s/s/s, aging rate a
        taus = [1 10 100 1000]';            % intervals for short-term stability
        stds = [1.5e-11 5e-12 1.5e-12 5e-13]';  % deviations for taus
        % 3-sigma deviation @ 10s (mm/s)
        dev = 1e3 * c * 3 * stds(2);

    case "SafranMAC"
        a = 5e-12 / (86400);                % s/s/s, aging rate a
        taus = [1 10 100]';                 % intervals for short-term stability
        stds = [4e-11 1.3e-11 4e-12]';      % deviations for taus
        % 3-sigma deviation @ 10s (mm/s)
        dev = 1e3 * c * 3 * stds(2);

    case "SafranMiniRAFS"
        a = 1e-10 / (365 * 86400);          % s/s/s, aging rate a
        taus = [1 10 100 1000]';            % intervals for short-term stability
        stds = [1e-11 3e-12 1e-12 3e-13]';  % deviations for taus
        % 3-sigma deviation @ 10s (mm/s)
        dev = 1e3 * c * 3 * stds(2);

    case "RakonMiniUSO"
        a = 1e-8 / (365 * 86400);           % s/s/s, aging rate a
        taus = [1 10 100 1000 1e4]';        % intervals for short-term stability
        stds = [2e-13 2e-13 4e-13 7e-13 7e-12]';  % deviations for taus
        % 3-sigma deviation @ 10s (mm/s)
        dev = 1e3 * c * 3 * stds(2);

    case "AccubeatUSO"
        a = 7e-11 / (86400);                % s/s/s, aging rate a
        taus = [1 10 100 1000]';            % intervals for short-term stability
        stds = [5e-13 5e-13 5e-13 6e-13]';  % deviations for taus
        % 3-sigma deviation @ 10s (mm/s)
        dev = 1e3 * c * 3 * stds(2);

    case "SafranRAFS"
        a = 1e-10 / (365 * 86400);          % s/s/s, aging rate a
        taus = [1 10 100 1000]';            % intervals for short-term stability
        stds = [3e-12 1e-12 3e-13 6e-14]';  % deviations for taus
        % 3-sigma deviation @ 10s (mm/s)
        dev = 1e3 * c * 3 * stds(2);

    case "ExcelitasRAFS"
        a = 5e-14 / (86400);                % s/s/s, aging rate a
        taus = [1 10 100 1e3 1e4 1e5]';     % intervals for short-term stability
        stds = 2e-12./sqrt(taus) + 2e-14;   % deviations for taus
        % 3-sigma deviation @ 10s (mm/s)
        dev = 1e3 * c * 3 * stds(2);

    otherwise
        error("getclockcoefficients:invalidClock", ...
            "%s not found, valid clock options found in " + ...
            "documentation. Strings must match exactly.", clock);
end

if strcmpi(fit, "Allan")
    B = stability2allan(taus, stds, a);
    coeff = allan2diffcoeff(B);
elseif strcmpi(fit, "Hadamard")
    B = stability2hadamard(taus, stds);
    coeff = hadamard2diffcoeff(B);
else
    error("getclockcoefficients:invalidFit", ...
        "%s not found, valid fit options found in documentation.", fit);
end

if plot
    % plot short-term stability fit to check quality
    ta = taus(1):taus(end);
    if strcmpi(fit, "Allan")
        sa = sqrt(B(1) ./ ta + B(2) .* ta + a^2/2 * ta .^ 2);
    else
        sa = sqrt(B(1) ./ ta + B(2) * ta + B(3) * ta .^ 3);
    end

    figure();
    loglog(ta, sa);
    hold on;
    scatter(taus, stds, 100, "rx");
    xlabel("Interval (s)");
    ylabel("\sigma_Y(\tau)");
end
end

