% ReceiverClassDemo.m
% Author: Mark Hartigan
% Date  : November 6, 2024
% Description:
%    Demonstrates how to use the Receiver class to evaluate receiver noise
%    on pseudorange and Doppler measurements.

%% reset
clc, clear, close all;

%% initialize receiver info
% receiver on user
rec = Receiver("PLL");
% change parameters from default as necessary
rec.freq = 1575.42e6;       % Hz, GPS L1 frequency
rec.Rc = 1.023e6;           % chips/s, L1 C/A chipping rate
rec.data = 1;               % tracking a channel w/ nav message
rec.codeorder = 3;          % DLL loop order
rec.Bfe = 2*rec.Rc;         % DLL front-end bandwidth
rec.Bn = 0.1;               % tighten noise bandwidth, DLL is carrier-aided
rec.T = 0.020;              % PIT as wide as nav msg bit transitions
rec.carrierorder = 3;       % PLL loop order
rec.T_c = 0.010;            % PIT half or bit transitions
rec.Tm = 10;                % delta-pseudorange @10s

%% generate receiver noise from selected data
% all inputs to Receiver.noise() must be the same size, as they should all
% be time-tagged by the input times ts
ts = 0:10;                  % usually in seconds past J2000 (SPICE time system)
CN0 = 15:3:45;              % assume C/N0 sweeps through range over ts
r = 4e8 * ones(1,11);       % range (in m) is static, roughly lunar distance
dr = zeros(1,11);           % zero line-of-sight speed (mm/s)
[err, var] = rec.noise(CN0, ts, r, dr);

%% plot info to look at it
figure();
yyaxis left;
plot(CN0, 3*sqrt(var(1,:)));
grid on;
xlabel("C/N0 (dB-Hz)");
ylabel("3\sigma pseudorange uncertainty (m)");
title("Test case receiver noise vs. signal strength");
yyaxis right;
plot(CN0, 3*sqrt(var(2,:)));
ylabel("3\sigma Doppler uncertainty (mm/s)");
