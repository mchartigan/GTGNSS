function handle = statetotrajectory(t0,tf,x0,dyn)
%STATETOTRAJECTORY Converts initial position and velocity of spacecraft to
%a spline-interpolated trajectory handle (from t0 to tf).
%   Inputs:
%    - t0; starting time, seconds past J2000
%    - tf; ending time, seconds past J2000
%    - x0; starting state [km, km/s]
%    - dyn; dynamics function @(t,x)

ts = linspace(t0, tf, 1000);
opts = odeset("RelTol", 1e-10, "AbsTol", 1e-12);
[~,X] = ode45(dyn, ts, x0, opts);
plotLunarOrbit(ts, X, 'J2000', "test");
pp = spline(ts, X');
handle = @(tau) ppval(pp, tau);
end

