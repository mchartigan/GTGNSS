% formatAntennaData.m
% Author: Mark Hartigan
% Date  : September 18, 2024
% Description:
%    Organize the GPS block IIR and IIR-M antenna pattern data into a
%    format more machine-readable. Also average the azimuth patterns,
%    making it only dependent on elevation, and subtract the Gain
%    Correction Factor (GCF).
%
% TODO:
%  - Create average antenna patterns for generic use, call SVNX_L1.txt and
%    SVNX_L2.txt; this should be made using the IIR-M data

%% reset
clc, clear, close all;
addpath(genpath(pwd));

%% file operations
sheets = dir(fullfile("data/gps/patterns", "*.xlsx"));
GCFdata = importdata("data/gps/patterns/GCF_from_SVN.txt");
GCFdata = GCFdata.data;

for i=1:length(sheets)
    % get sheet name and extract info
    sheet = sheets(i);
    name = sheet.name;
    temp = split(name, {'-','.'});
    SVN = temp(1);                      % name in format 'SVNXX'
    L = temp(2);                        % transmit band, 'L1' or 'L2'
    newname = SVN + "-" + L + ".txt";   % filename for formatted data

    % data formatting; average phi angles and +/- theta
    data = importdata(name);
    data = data.data;
    rows = ~isnan(data(:,1));
    data = data(rows,:);                % remove header

    % average directivity over power, not dB, then convert back to dB
    % avgphi_old = mean(data(:,2:end), 2);
    avgphi = mean(10.^(data(:,2:end)/10), 2);   % convert to power/isotropic
    mid = ceil(length(avgphi)/2);
    avgtht = mean([avgphi(mid:end) flipud(avgphi(1:mid))], 2);
    avgtht = 10 * log10(avgtht);                % convert back to dBi
    thetas = data(mid:end,1);

    % subtract gain correction factor to arrive at antenna gain
    SVN = split(SVN, 'SVN');
    SVN = str2double(SVN{2});
    j = find(GCFdata(:,1) == SVN, 1);
    if strcmpi(L, 'L1'), k = 2; else, k = 3; end

    if isempty(j) || isempty(k)
        error("formatAntennaData:GCFnotFound", ...
            "Couldn't find GCF for SVN %d, %s.", SVN, L{1});
    else
        avgtht = avgtht - GCFdata(j,k);
    end

    % open and write to file
    id = fopen(sheet.folder + "/" + newname, 'w');
    fprintf(id, "theta (deg),gain (dBi)\n");
    for j=1:mid
        fprintf(id, "%f,%f\n", thetas(j), avgtht(j));
    end
    fclose(id);
end
    