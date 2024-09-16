% PlotFromYuma.m
% Author: Mark Hartigan
% Date  : August 26, 2024
% Description:
%    Plot GPS satellite orbits over time, given a Yuma file

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% initialize
% load SPICE kernels
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_earth.tm'));
% get OE data
oes = yuma2oes("data/gps/almanac.yuma.week0278.061440.txt");
t0 = oes(1).t0;
% create propagator object (8x8 gravity)
prop = EarthPropagator(t0, oes, 8, 1, "opts", odeset("RelTol", 1));

%% run propagator
prop.run(3600 * 12, 1000, 'J2000');
prop.plotlastorbits('J2000');
funcs = prop.statetotrajectory();

% %% get PRN -> SV info
% % assign rough corresponding datenum for use with ODTBX
% date = datetime(t0, 'ConvertFrom', 'epochtime', 'Epoch', '2000-01-01 12:00:00');
% % use ODTBX to obtain PRN -> SV conversions
% odtbxopts = setOdtbxOptions('epoch', datenum(date));
% txinfo = get_gps_block(odtbxopts);
% 
% % assign SV # based on PRN #
% for i=1:length(oes)
%     oes(i).SV = txinfo{oes(i).PRN}.GPS_ID;
% end