function epoch = gpstow2j2000(week,secs)
%GPSTOW2J2000 Converts a GPS week and seconds of week to seconds past J2000
%(NAIF SPICE's central time epoch).
%   Input:
%    - week; GPS week (rollover must be added)
%    - secs; seconds past start of GPS week

epoch = (week * 7 - 7300.5) * 86400 + secs;
epoch = cspice_unitim(epoch, 'GPS', 'TDB');
end

