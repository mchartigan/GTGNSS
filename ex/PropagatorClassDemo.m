% PropagatorClassDemo.m
% Author: Mark Hartigan
% Date  : November 7, 2024
% Description:
%    Demonstrates how to use the OrbitPropagator classes (EarthPropagator / 
%    LunarPropagator) in various ways.

%% reset
clc, clear, close all;

%% load SPICE kernels
spkpath = fileread('SPKPATH');
cspice_furnsh(strcat(spkpath,'/kernels/generic/mk/generic_earth.tm'));
cspice_furnsh(strcat(spkpath,'/kernels/generic/mk/generic_lunar.tm'));

%% initialize GPS propagator from Yuma almanac
% get GPS OE data
oes = yuma2oes("almanac.yuma.week0278.061440.txt");
t0 = oes(1).t0;                         % grab starting time from OE struct
% create propagator for GPS satellites from Keplerian orbital elements (assumed 
% in ECI) with 8x8 gravity and lunar 3rd body
gpsprop = EarthPropagator(t0, oes, 8, 1);

%% run and plot GPS orbits
% propagate orbits for 12 hours, 1000 data points, in ECI
[t_gps, x_gps] = gpsprop.run(3600 * 12, 1000, 'J2000');
gpsprop.plotlastorbits('J2000');        % plot orbits in ECI

%% initialize lunar satellite propagator from initial state
% generate orbit in frozen orbit frame from OEs (MOON_OP in SPICE)
a = 12000;
e = 0.8;
i = 60 * pi/180;
RAAN = 0;
w = 90 * pi/180;
f = 0;
% convert orbital elements to pos/vel vectors in MOON_OP frame
[r_op, v_op] = oe2rv(a, e, i, RAAN, w, f, cspice_bodvrd('MOON', 'GM', 1));
% convert state in MOON_OP frame to Moon-centered ICRF
x_mci = cspice_sxform('MOON_OP', 'J2000', t0) * [r_op; v_op];
% create propagator for lunar satellite from state in ICRF with 32x32
% gravity and Earth + sun 3rd body
lunarprop = LunarPropagator(t0, x_mci, 32, 2);

%% run and plot lunar orbit
% propagate and get data once/min for 24 hours
ts = 0:60:86400;
[~,x_lunar] = lunarprop.runat(ts, 'J2000');
lunarprop.plotlastorbits('MOON_OP');    % plot in frozen orbit frame