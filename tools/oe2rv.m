function [r_, v_] = oe2rv(a, e, i, RAAN, w, f, mu)
%OE2RV computes the ECI position and velocity vectors associated with
%Keplerian orbital elements.
% Input:
%  - a; semimajor axis [km]
%  - e; eccentricity
%  - i; inclination [rad]
%  - RAAN; right ascension of the ascending node [rad]
%  - w; argument of periapsis [rad]
%  - f; true anomaly [rad]
%  - mu; gravitational parameter of body being orbited, [km^3 / s^2]
%
% Alternate input:
%  - oe; struct of orbital elements {a,e,i,RAAN,w,f}
%  - mu; same as above

% parse optional struct input
if nargin == 2
    mu = e;
    e = a.e;
    i = a.i;
    RAAN = a.RAAN;
    w = a.w;
    f = a.f;
    a = a.a;
end

p = a*(1 - e^2);                    % semi-latus rectum [km]
r = p / (1 + e*cos(f));             % radial position [km]
dr = sqrt(mu/p) * e * sin(f);       % dr/dt [km / s]
rdf = sqrt(mu/p) * (1 + e*cos(f));  % r * df/dt [km / s]

r_p = [r*cos(f) r*sin(f) 0]';
v_p = [dr * cos(f) - rdf * sin(f); dr * sin(f) + rdf * cos(f); 0];

R_P2ECI = [cos(-RAAN) sin(-RAAN) 0; -sin(-RAAN) cos(-RAAN) 0; 0 0 1] * ...
          [1 0 0; 0 cos(-i) sin(-i); 0 -sin(-i) cos(-i)] * ...
          [cos(-w) sin(-w) 0; -sin(-w) cos(-w) 0; 0 0 1];

r_ = R_P2ECI * r_p;
v_ = R_P2ECI * v_p;

% accomodate outputting as a single state vector
if nargout == 1
    r_ = [r_; v_];
end
end

