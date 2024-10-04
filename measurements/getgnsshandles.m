function [fGnss, namesGnss] = getgnsshandles(includeGnss)
%GETGNSSHANDLES Returns function handles for each satellite in includeGnss
%constellation; handle accepts time(s) in seconds past J2000 and returns
%[pos (km); vel (km/s)] state(s).
%   Input:
%    - includeGnss; "GPS", "GALILEO", or "BOTH"
%   Output:
%    - fGnss; cell array of fcn handles corresponding to namesGnss sats
%    - namesGnss; SPICE names of GNSS satellites included
arguments
    includeGnss (1,1) string
end

% GNSS SPICE kernel dataset (SKD) provided by the ESA SPICE Service:
% https://spiftp.esac.esa.int/data/SPICE/GNSS/misc/gnss.html
%
% Data provided from 13 Nov 2022 to 18 Apr 2024
% cspice_furnsh(strcat(userpath,'/kernels/GNSS/mk/gnss.tm'));

% load GNSS trajectories
namesGnss = {};                 % SPICE names of GNSS satellites
if strcmpi(includeGnss, "GPS") || strcmpi(includeGnss, "BOTH")
    namesGnss = [namesGnss; importdata(strcat(userpath,'/kernels/GNSS/GPS_sats.txt'))];
end
if strcmpi(includeGnss, "GALILEO") || strcmpi(includeGnss, "BOTH")
    namesGnss = [namesGnss; importdata(strcat(userpath,'/kernels/GNSS/Galileo_sats.txt'))];
end
nGnss = length(namesGnss);      % # of GNSS satellites to include
if ~nGnss                       % throw error if none selected
    error("main:noGNSS", "No GNSS satellites included. Ensure include_GNSS is set properly.");
end

fGnss = cell(nGnss, 1);
for i=1:nGnss
    fGnss{i} = @(tau) cspice_spkezr(namesGnss{i}, tau, 'J2000', 'NONE', 'MOON');
end
end

