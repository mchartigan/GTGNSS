function [oes,t0] = gal2oes(file,rollover,tol)
%GAL2OES Converts data from a Galileo Almanac (XML) into a list of SV OE 
%structs. Requires SPICE time kernel to be loaded.
%
%   https://www.gsc-europa.eu/gsc-products/almanac#parameters
%
%   OE structs will contain the following information:
%   {  SV, space vehicle #
%      PRN, PRN code #
%      t0, time information was collected at (in s past J2000, UTC/TDB)
%      a, semi-major axis (in km)
%      e, eccentricity
%      i, orbit inclination (in rad)
%      RAAN, right ascension of the ascending node (in rad)
%      w, argument of perigee (in rad)
%      f, true anomaly at starting time (in rad)
%      af0, clock bias (in s)
%      af1, clock drift rate (in s/s) }
%
%   Inputs:
%    - file; path to Yuma Almanac
%    - rollover; number of GPS rollover weeks (default 2 since 2019)
%    - tol; newton iteration tolerance for finding eccentric anomaly
%   Output:
%    - oes; list of OE structs
%    - t0; start time, seconds past J2000
arguments
    file     (1,:) {mustBeText}
    rollover (1,1) {mustBeInteger,mustBeNonnegative} = 2
    tol      (1,1) double = 1e-9
end

% farm out actual XML reading
data = galalmanacread(file);
% convert datetime to UTC string for SPICE
START = convertStringsToChars(sprintf("%s", data(1,1).Time) + " UTC");
t0 = cspice_str2et(START);

n = size(data, 1);      % number of SVs recorded
blank = cell(1,n);
oes = struct('SV',blank,'PRN',blank,'t0',blank,'a',blank,'e',blank,'i',blank, ...
             'RAAN',blank,'w',blank,'f',blank,'af0',blank,'af1',blank);

for j=1:n
    el = data(j,:);
    oes(j).t0 = t0;
    oes(j).SV = el.SVID;
    tk = el.t0a;

    % aSqRoot is offset w.r.t. sqrt of nominal a (29,600 km)
    oes(j).a = (sqrt(29600000) + el.aSqRoot)^2 / 1000;
    oes(j).e = el.ecc;
    % deltai is offset in semicircles (1 semicirc. = pi rad) from i (56 deg)
    oes(j).i = (56/180 + el.deltai) * pi;
    % RAAN is RAAN0 (at weekly epoch) + dRAAN * seconds of week
    oes(j).RAAN = (el.omega0 + el.omegaDot*tk) * pi;
    oes(j).w = el.w * pi;

    % find true anomaly
    M = el.m0 * pi;
    % use Newton iteration to solve for eccentric anomaly corresp. to M
    E = M;
    for k=1:100
        E = E - (E - oes(j).e*sin(E) - M) / (1 - oes(j).e*cos(E));
        if abs(E - oes(j).e*sin(E) - M) < tol
            break;
        elseif k >= 100
            error('gal2oes:ESolver', ...
                'Failed to converge to E in %d iterations.', k);
        end
    end
    % find the true anomaly corresp. to eccentric anomaly
    oes(j).f = 2 * atan2(sqrt(1 + oes(j).e) * tan(E / 2), sqrt(1 - oes(j).e));

    oes(j).af0 = el.af0;
    oes(j).af1 = el.af1;
end
end