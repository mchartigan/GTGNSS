function [oes,t0] = yuma2oes(file,rollover,tol)
%YUMA2OES Converts data from a Yuma Almanac into a list of SV OE structs.
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

% convert GPS week and seconds of week to seconds past J2000, TDB
temp = split(file, "week");
temp = split(temp(2), ".");
week = double(temp(1)) + 1024*rollover;
secs = double(temp(2));
t0 = gpstow2j2000(week, secs);

% gather and parse lines
lines = readlines(file);
m = length(lines);          % lines in file
n = floor(m / 15);          % number of SVs recorded (15 lines ea.)
blank = cell(1,n);
oes = struct('SV',blank,'PRN',blank,'t0',blank,'a',blank,'e',blank,'i',blank, ...
             'RAAN',blank,'w',blank,'f',blank,'af0',blank,'af1',blank, ...
             'gcfL1',blank,'gcfL2',blank);

% iterate over each bank of data in file
for j=1:n
    oes(j).t0 = t0;
    PRN = strip(split(lines((j-1)*15 + 2), ":"));
    oes(j).PRN = double(PRN(2));        % satellite PRN #
    e = strip(split(lines((j-1)*15 + 4), ":"));
    oes(j).e = double(e(2));            % orbit eccentricity
    tk = strip(split(lines((j-1)*15 + 5), ":"));
    tk = double(tk(2));                 % time past GPS week info was collected
    i = strip(split(lines((j-1)*15 + 6), ":"));
    oes(j).i = double(i(2));            % orbit inclination
    dRAAN = strip(split(lines((j-1)*15 + 7), ":"));
    dRAAN = double(dRAAN(2));           % drift rate of RAAN
    a = strip(split(lines((j-1)*15 + 8), ":"));
    oes(j).a = double(a(2))^2 / 1000;   % orbit semi-major axis
    % right ascension from file is given at start of GPS week, so must
    % modify with drift rate and time past start of week
    RAAN0 = strip(split(lines((j-1)*15 + 9), ":"));
    oes(j).RAAN = double(RAAN0(2)) + dRAAN * tk;
    w = strip(split(lines((j-1)*15 + 10), ":"));
    oes(j).w = double(w(2));            % orbit argument of perigee
    M = strip(split(lines((j-1)*15 + 10), ":"));
    M = double(M(2));                   % mean anomaly

    % use Newton iteration to solve for eccentric anomaly corresp. to M
    E = M;
    for k=1:100
        E = E - (E - oes(j).e*sin(E) - M) / (1 - oes(j).e*cos(E));
        if abs(E - oes(j).e*sin(E) - M) < tol
            break;
        elseif k >= 100
            error('yuma2oes:ESolver', ...
                'Failed to converge to E in %d iterations.', k);
        end
    end
    % find the true anomaly corresp. to eccentric anomaly
    oes(j).f = 2 * atan2(sqrt(1 + oes(j).e) * tan(E / 2), sqrt(1 - oes(j).e));

    % get SV clock parameters
    af0 = strip(split(lines((j-1)*15 + 12), ":"));
    oes(j).af0 = double(af0(2));
    af1 = strip(split(lines((j-1)*15 + 12), ":"));
    oes(j).af1 = double(af1(2));
end
end

