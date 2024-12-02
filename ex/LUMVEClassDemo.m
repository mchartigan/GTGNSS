% LUMVEClassDemo.m
% Author: Mark Hartigan
% Date  : November 7, 2024
% Description:
%    A sample script (nonfunctional) that demonstrates how to propagate
%    orbits and construct a navigation filter.

%% reset
clc, clear, close all;

%% initialize
% load SPICE kernels
spkpath = fileread('SPKPATH');
cspice_furnsh(strcat(spkpath,'/kernels/generic/mk/generic_earth.tm'));

% create simulation starting time
START = '2024 Nov 7 00:00:00';
t0 = cspice_str2et(START);

%% create LEO satellite from sample orbital elements
leo.a = 6880;
leo.e = 0;
leo.i = 45 * pi/180;
leo.RAAN = 0;
leo.w = 0;
leo.f = 0;

% run propagator
leoprop = EarthPropagator(t0, leo, 32, 1);
[~, xs] = leoprop.run(3600 * 12, 1000, 'J2000');
leoprop.plotlastorbits('J2000');

%% gather measurements
% this is where measurements are implemented, assuming none at the moment
h = @(t,x) 0;       % measurement model
y = [];             % dataset of measurements
t_meas = [];
dydx = @(t,x) 0;    % jacobian of measurements w.r.t. state
R = 0;              % measurement noise matrix

%% example filter demo
% this is sort of a nothing burger of a navigation simulation as we have no
% measurements, have perfect knowledge of the starting state, and perfect
% knowledge of the dynamics. Therefore, any error we see is just
% propagation errors between the "truth" (when we run leoprop) and the
% "estimate" (when we run the filter). The point here is to illustrate how
% to construct the filter and run it.

% initialize batch filter (weighted least squares)
filter = LUMVE("hybrid", 6, @leoprop.dynamics, @leoprop.partials, ...
    h, y, dydx, R, 1, t_meas, "opts", leoprop.opts, "t_sim", leoprop.ts);
% run the filter w/ a priori info
filter.run(leoprop.x0, eye(6));

% diff the two runs and plot the errors
diff = filter.x - xs;

figure();
plot((filter.t - filter.t(1))/3600, sum(diff(1:3,:), 1));
grid on;
xlabel("Time (hrs)");
ylabel("error (km)");
title("Position error between filter and truth data");
