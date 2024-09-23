% GNSSFilteringDemo.m
% Author: Mark Hartigan
% Date  : August 26, 2024
% Description:
%    A sample script (nonfunctional) that demonstrates how to propagate
%    orbits and construct a navigation filter.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% initialize
% load SPICE kernels
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_earth.tm'));

%% create GNSS satellites (8x8 gravity)
% get OE data
oes = yuma2oes("almanac.yuma.week0278.061440.txt");
t0 = oes(1).t0;
gpsprop = EarthPropagator(t0, oes, 8, 1);

% run propagator
gpsprop.run(3600 * 12, 1000, 'J2000');
gpsprop.plotlastorbits('J2000');
funcs = gpsprop.statetotrajectory();

%% create LEO satellite
% sample orbital elements
leo.a = 6880;
leo.e = 0;
leo.i = 45 * pi/180;
leo.RAAN = 0;
leo.w = 0;
leo.f = 0;

% run propagator
leoprop = EarthPropagator(t0, leo, 32, 1);
leoprop.run(3600 * 12, 1000, 'J2000');
leoprop.plotlastorbits('J2000');

%% gather measurements
% this is where GNSSmeasurements() would be used, no time for
% implementation rn
h = @(t,x) 0;       % measurement model
y = [];             % dataset of measurements
t_meas = [];
dydx = @(t,x) 0;    % jacobian of measurements w.r.t. state
R = 0;              % measurement noise matrix

%% example filter demo
filter = LUMVE("hybrid", 6, @leoprop.dynamics, @leoprop.partials, ...
    h, y, dydx, R, 1, t_meas, "opts", leoprop.opts, "t_sim", leoprop.ts);
filter.run(leoprop.x0, eye(6));

diff = filter.x - leoprop.xs;

figure();
plot((filter.t - filter.t(1))/3600, sum(diff(1:3,:), 1));
grid on;
xlabel("Time (hrs)"); ylabel("error (km)");
title("Position error between filter and truth data");
