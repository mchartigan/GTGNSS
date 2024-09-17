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
oes = yuma2oes("almanac.yuma.week0278.061440.txt");
t0 = oes(1).t0;
% create propagator object (8x8 gravity)
prop = EarthPropagator(t0, oes, 8, 1, "opts", odeset("RelTol", 1));

%% run propagator
prop.run(3600 * 12, 1000, 'J2000');
prop.plotlastorbits('J2000');
funcs = prop.statetotrajectory();
