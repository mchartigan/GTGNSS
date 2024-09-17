% errorBudgetSimulation.m
% Author: Mark Hartigan
% Date  : June 10, 2024
% Description:
%    Chaining together simulations to develop a demonstration of error
%    budgets for lunar navigation systems.
%
%    NOTE: THIS SCRIPT WILL NOT EXECUTE PROPERLY WITH THE BASE INSTALLATION
%    OF GTGNSS. IT SHOULD BE USED FOR REFERENCE ONLY ON COMMAND
%    IMPLEMENTATION.

%% reset
clc, clear, close all;
addpath(genpath(pwd));
format long g;          % display long numbers, no scientific notation

%% init
% starting epoch
START = '2024 May 1 00:00:00';
% frame of data
FRAME = 'J2000';
% number of days in simulation
DAYS = 0.5;
% number of spherical harmonics to use in EKF
N_SPH = 16;
% should filter performance statistics and plots be displayed?
FILTER_PERFORMANCE = true;

% load generic lunar information
cspice_furnsh(strcat(userpath,'/kernels/generic/mk/generic_lunar.tm'));
% load ELFO orbit 
cspice_furnsh('data/gmat-to-spk/ELFO_20240501-20240531.bsp');
% load GNSS satellite info
cspice_furnsh(strcat(userpath,'/kernels/GNSS/mk/gnss.tm'));

ELFO = '-909';                  % SPICE ID of ELFO s/c
t0 = cspice_str2et(START);

% planetary info
bods = getplanets('MOON', "MOON", "EARTH", "SUN", "JUPITER");

% store spherical harmonic coefficients for the moon from Lunar Prospector
[~,C,S] = cofloader("LP165P.cof");
bods(1).C = C; bods(1).S = S;
bods(1).frame = 'MOON_ME';      % body-fixed frame of coefficients
moon = bods(1);                 % primary body
sec = bods(2:end);              % secondary bodies

%% generate truth orbit
ns = 1440;
ts = linspace(t0, t0+86400*DAYS, ns);
sat = zeros(9,ns);

for i=1:ns
    sat(1:6,i) = cspice_spkezr(ELFO, ts(i), FRAME, 'NONE', 'MOON');
end

x0_OP = cspice_sxform(FRAME, 'MOON_OP', t0) * sat(1:6,1);
[oe.a,oe.e,oe.i,oe.RAAN,oe.w,oe.f] = rv2oe(x0_OP(1:3), x0_OP(4:6), moon.GM);
% plotLunarOrbit(ts, sat', FRAME, "ELFO");

a = 2e-10 / (86400 * 30);       % Hz/Hz/s, upper end of aging rate a
x0Clk = [-0.8e-7, 0.5e-10, a];  % starting state of clock
[s1,s2,s3] = DiffCoeffCSAC();
xClk = clockStateOverTime(ts, x0Clk, "CSAC");
sat(7:9,:) = xClk;
x0 = sat(:,1);

%% get GNSS satellite setup
[f_GNSS, names_GNSS] = getgnsshandles("BOTH");
n_GNSS = length(f_GNSS);

% generate measurement model
R = zeros(2*n_GNSS,2*n_GNSS);       % measurement covariance matrix
% per GPS specifications
R(1:n_GNSS,1:n_GNSS) = diag(repmat((9.7/1.96*1e-3)^2, 1, n_GNSS));
% per GPS specifications
R(n_GNSS+1:end,n_GNSS+1:end) = diag(repmat((.006/1.96*1e-3)^2, 1, n_GNSS));

measurements = GNSSmeasurements("BOTH",ts,f_GNSS,R,sec(1),moon,sat);

%% create filter
% process noise
% no uncertainty about velocity, accel ~4e-7 from experimental_process_noise.m
% stds = 3.98642855800357e-07 3.34295767996475e-07 3.81824183418184e-07
Qorb = [0 0 0 1 1 1] * (2e-6)^2;  
Qclk = [s1 s2 s3].^2;
Q = diag([Qorb Qclk]);

opts = odeset("RelTol", 1e-9, "AbsTol", 1e-11);
filter = EKF("hybrid", @(t,x) lnss_elfodyn(t,x,moon,N_SPH,sec), ...
             @(t,x) lnss_elfopartials(t,x,moon,sec(1)), ...
             Q, ...
             @(t,x) measurements.compute(t,x), ...
             measurements.ymeas, ...
             @(t,x) measurements.partials(t,x), ...
             R, ts, "opts", opts);

% run filter
P0hat = diag([.03 .03 .03 3e-10 3e-10 3e-10 7e-13 1e-21 1e-32]);
x0hat = mvnrnd(x0', P0hat)';
filter.run(x0hat, P0hat);

%% convert error data to RTN
x_err = filter.x - sat;
x_err_rtn = x_err;
P_rtn = filter.P;
for i=1:ns
    ti = ts(i);
    r_r = sat(1:3,i);
    u_r = r_r / norm(r_r);          % radial direction
    r_n = cross(r_r, sat(4:6,i));
    u_n = r_n / norm(r_n);          % normal direction
    u_t = cross(u_n, u_r);          % tangential direction
    T_J2RTN = [u_r u_t u_n]';
    T = [T_J2RTN zeros(3,3); zeros(3,3) T_J2RTN];
    x_err_rtn(1:6,i) = T * x_err(1:6,i);
    P_rtn(1:6,1:6,i) = T * filter.P(1:6,1:6,i) * T';
end

%% performance statistics
if FILTER_PERFORMANCE
    filterperformance("EKF", ts, x_err, filter.P, ...
                      x_err_rtn, P_rtn, measurements.visible);
end

%% ephemeris fit
t02 = ts(end);              % starting point of ephemeris propagation
x02 = filter.x(:,end);      % starting state of ephemeris propagation
hrs = 4;
tf2 = 3600*hrs;             % final bound of ephemeris propagation
[a,e,i,RAAN,w,~] = rv2oe(x02(1:3), x02(4:6), moon.GM);
[r0,v0] = oe2rv(6540, 0.6, 56.2*pi/180, pi, pi/2, -pi/4, moon.GM);
x02 = [r0; v0];

% generate ephemeris approximation
x0_eph = cspice_sxform('J2000', 'MOON_OP', t02) * x02;
N_APPX = 14;
ephprop = LunarPropagator(t02, x0_eph, 50, 3);
f_eph1 = ephprop.ephemerisfit("polynomial", tf2, N_APPX);
f_eph2 = ephprop.ephemerisfit("Kepler", tf2, N_APPX - 2);

% truth orbit
[ts2,x_true2] = ephprop.run(tf2, 1600, 'J2000');
ephprop.plotlastorbits('MOON_OP');
% error between truth and approximations
err1_x = f_eph1(ts2) - x_true2;
err1_pos = sqrt(sum(err1_x(1:3,:).^2, 1)) * 1e3;
err1_vel = sqrt(sum(err1_x(4:6,:).^2, 1)) * 1e6;
err2_x = f_eph2(ts2) - x_true2;
err2_pos = sqrt(sum(err2_x(1:3,:).^2, 1)) * 1e3;
err2_vel = sqrt(sum(err2_x(4:6,:).^2, 1)) * 1e6;

fprintf("  int. |  Poly (3s) |  Kep  (3s)\n");
fprintf("--------------------------------\n");
fprintf(" %d hrs | %8.2f m | %8.2f m \n", hrs, 3*std(err1_pos), 3*std(err2_pos));
fprintf("       | %5.2f mm/s | %5.2f mm/s \n", 3*std(err1_vel), 3*std(err2_vel));

plotformat("IEEE", 1, "scaling", 2);
figure();
subplot(2,1,1);
plot((ts2 - t02)/3600, err1_pos);
hold on;
plot((ts2 - t02)/3600, err2_pos);
hold off; grid on;
linestyleorder("mixedstyles");
ylabel("Error (m)");
title("Approximation error of different surrogate models");

subplot(2,1,2);
plot((ts2 - t02)/3600, err1_vel);
hold on;
plot((ts2 - t02)/3600, err2_vel);
hold off; grid on;
linestyleorder("mixedstyles");
xlabel("Time (hrs)"); ylabel("Error (mm/s)");
legend(["Polynomial (N=14)", "Keplerian + Polynomial (N=12)"], ...
    "location", "northwest");


% TODO:
% - compare / contrast with max approximation intervals bar chart from
%   cortinovis
% - add additional observable (ranging to moon simulating opnav?) for nav
%   solution improvement
