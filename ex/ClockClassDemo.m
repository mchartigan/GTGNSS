% ClockClassDemo.m
% Author: Mark Hartigan
% Date  : November 7, 2024
% Description:
%    Demonstrate use of the Clock class and generate various plots. This
%    script in particular is adapted from one to demonstrate why the
%    Microchip CSAC won't work as the onboard oscillator for NSNS
%    satellites -- its short-term stability is too poor to meet the SISE
%    velocity requirement of 1.2 mm/s 3-sigma @10s measurement intervals.

%% reset
clc, clear, close all;

%% initialize clock info
c = 299792458;                              % m/s, speed of light

% create Clock instance
clk = Clock(0, [0 0 0], "SafranMiniRAFS");  % oscillator propagator
% set starting state again since now we've gotten the aging rate
clk.x0 = [0 0 clk.a];
a = clk.a;
% compute the 3-sigma allan deviation (in mm/s) @10s
adev = c * 1e3 * 3 * clk.s_allan(2);

%% generate data
ts = 0:1:3600*2;                        % data once/sec for 2 hours
n = length(ts);
[~, xs, vs] = satclock.runat(ts);       % run to get data
% separate data into phase and frequency offsets, minus initial model
% propagation (just aging) to get error
b = c * (xs(1,:) - 0.5 * a * ts.^2);
f = c * 1e3 * (xs(2,:) - a * ts);
% 3-sigma uncertainty of offsets
v_b = c * 3 * sqrt(reshape(vs(1,1,:), 1, n));
v_f = c * 1e3 * 3 * sqrt(reshape(vs(2,2,:), 1, n));

% can't directly measure frequency -- rather, it's two time-differenced
% phase measurements -- so add frequency stability noise (Allan variance)
% from measurement period (@10s per Lunar Relay SRD)
v_f = 3 * sqrt((v_f/3).^2 + (adev/3)^2);
f = f + mvnrnd(0, (adev/3)^2, n)';

%% plot data
% plotformat("IEEE", 1, "coloring", "greyscale");
tp = ts ./ 60;                          % min, time for plotting
colors = [0   0   0
          0.3 0.3 0.3
          0.6 0.6 0.6];                 % greyscale for plotting

% plot phase measurement (no receiver error), uncertainty, and NSNS requirement
figure();
subplot(2,1,1);
patch([tp flip(tp)], [-v_b flip(v_b)], colors(3,:), ...
      "FaceAlpha", 0.6, "EdgeColor", "none");
hold on;
% NSNS requirement is 13.43m 3-sigma for ALL SISE errors; so in reality the
% clock will need to perform better than this.
patch([tp(1) tp(end) tp(end) tp(1)], [13.43 13.43 -13.43 -13.43], ...
      colors(2,:), "FaceAlpha", 0.6, "EdgeColor", "none");
plot(tp, b, "Color", colors(1,:));
hold off; grid on;
ylabel("Phase error (m)");
legend(["3\sigma Bound", "NSNS Req.", "Error"], "location", "best");

% plot freq. measurement (no receiver error), uncertainty, and NSNS requirement
subplot(2,1,2);
patch([tp flip(tp)], [-v_f flip(v_f)], colors(3,:), ...
      "FaceAlpha", 0.6, "EdgeColor", "none");
hold on;
% NSNS requirement is 1.2 mm/s 3-sigma @10s for ALL SISE errors; so in reality 
% the clock will need to perform better than this.
patch([tp(1) tp(end) tp(end) tp(1)], [1.2 1.2 -1.2 -1.2], colors(2,:), ...
      "FaceAlpha", 0.6, "EdgeColor", "none");
plot(tp, f, "Color", colors(1,:));
hold off;
grid on;
axis([-inf inf -125 125]);
xlabel("Time (min)");
ylabel("Frequency error (mm/s)");
sgtitle("Microchip Space CSAC dynamics, trial simulation");