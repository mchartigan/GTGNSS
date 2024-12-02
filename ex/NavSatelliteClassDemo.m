% NavSatelliteClassDemo.m
% Author: Mark Hartigan
% Date  : November 7, 2024
% Description:
%    Demonstrate the functionality of the NavSatellite class by simulating
%    the error budget to a user from a lunar NSNS node.

%% reset
clc, clear, close all;

%% initialize
% load SPICE kernels
spkpath = fileread('SPKPATH');
cspice_furnsh(strcat(spkpath,'/kernels/generic/mk/generic_lunar.tm'));

% get timing info
START = '2027 Feb 2 00:00:00';
END = 5 * 3600;
t0 = cspice_str2et(START);

% get lunar satellite OE data
oe.a = 9850; oe.e = 0.714; oe.i = 123 * pi/180; oe.RAAN = pi; oe.w = pi/2; oe.f = -pi/2;

% orbital propagator
satprop = LunarPropagator(t0, oe, 64, 3);
% clock info
a = 5e-14 / 86400;
x0Clk = [0 0 a];
satclock = Clock(t0, x0Clk, "RAFS");    % oscillator propagator

%% filter info
% assume s/c is getting ephemeris and time from the ground every X minutes
meas.h = @(~,x) x(1:8);
meas.dhdx = @(~,x) [eye(8) zeros(8,1); zeros(1,9)];
% measurement uncertainty
meas.R = diag([1e-3 1e-3 1e-3 1e-6 1e-6 1e-6 1e-8 1e-11]);
% measurements times (past start), every 5 min for 12 hours
tm = 0:60*5:END;
[~,xs] = satprop.runat(tm, 'J2000');
[~,xc] = satclock.runat(tm);
meas.y = [xs; xc(1:2,:)];
meas.t_meas = tm + t0;
meas.varargin = {"opts", odeset("RelTol", 1e-9, "AbsTol", 1e-11)};

%% antenna info
% transmitting antenna on spacecraft
transant = TransmitAntenna();
% dBi, gain of transmit antenna (based on IIR-M antenna at 0deg elev)
transant.gain = 14.75;
% receiving antenna on user
recant = ReceiveAntenna();
% dBi, receive antenna gain (based on ODTBX sensysmeas_ant.txt peak gain, 157deg half-beamwidth)
recant.gain = 4;
% receiver on user
rec = Receiver("PLL");
rec.codeorder = 3;          % increase DLL loop order
rec.carrierorder = 3;       % increase PLL loop order
rec.data = 0;               % tracking pilot channel
rec.Bn = 0.1;               % tighten noise bandwidth cuz DLL is carrier-aided
rec.T = 1;                  % increase T cuz channel is dataless
rec.T_c = 0.01;
rec.Tm = 10;                % measurements @10s, per LNSP specifications

%% create satellite
sat = NavSatellite(satprop, satclock, "EKF", meas, transant, false);
% simulation time steps; start at 10 so times aren't < 0 when accounting for TOF
% and delta-pseudorange computation
ts = (11:15:END) + t0;
n = length(ts);         % number of time steps

% build user
pos = [0 0 -1736 0 0 0]';       % on lunar south pole
user = User(@(t) repmat(pos,1,length(t)), 'MOON_ME', rec, recant);

% build options structure
navopts.x0 = [];
navopts.P0 = diag([1e-3 1e-3 1e-3 1e-6 1e-6 1e-6 1e-10 1e-13 1e-19].^2);
navopts.cadence = [0.5 0.5] * 3600;

%% get measurements to user
% plotformat("IEEE", 1);
[meas,var,true] = sat.getmeasurements(ts, user, navopts);

psr_sandpile = [var.psr.prop; var.psr.mdl; var.psr.clk; var.psr.rec];
psrr_sandpile = [var.psrr.prop; var.psrr.mdl; var.psrr.clk; var.psrr.rec];
labels = ["Propagation", "Model", "Clock", "Receiver"];
tplot = (ts - ts(1)) / 3600;
psr_3s = 3*sqrt(var.psr.total);
psrr_3s = 3*sqrt(var.psrr.total);

%% analysis
colors = [0   0   0
          0.3 0.3 0.3
          0.6 0.6 0.6];                 % greyscale for plotting

% plot the true measurement error, including the user-computed 3-sigma bound
figure();
subplot(2,1,1);
plot(tplot, true.err.psr.total, "Color", colors(1,:));
hold on;
patch([tplot flip(tplot)], [-psr_3s flip(psr_3s)], colors(2,:), ...
    "FaceAlpha", 0.6, "EdgeColor", "none");
hold off; grid on;
ylabel("Error (m)");
subplot(2,1,2);
plot(tplot, true.err.psrr.total, "Color", colors(1,:));
hold on;
patch([tplot flip(tplot)], [-psrr_3s flip(psrr_3s)], colors(2,:), ...
    "FaceAlpha", 0.6, "EdgeColor", "none");
hold off;grid on;
xlabel("Time (hrs)");
ylabel("Error (mm/s)");
legend(["Error", "3\sigma bound"], "location", "best");
sgtitle("User measurement error and uncertainty over time");

% plot a covariance sandpile for all the error contributors
figure();
subplot(2,1,1);
NavSatellite.sandpile(ts, psr_sandpile, {}, "\sigma^2 (m)^2", false);
subplot(2,1,2);
NavSatellite.sandpile(ts, psrr_sandpile, labels, "\sigma^2 (mm/s)^2");
sgtitle("Covariance sandpile for measurement uncertainty");
