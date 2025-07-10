function [a, e, i, RAAN, w, f] = rv2oe(r_, v_, mu)
%RV2OE computes the Keplerian orbital elements associated with the position
%and velocity vectors provided.
% Input:
%  - r_; position (3x1), [km]
%  - v_; velocity (3x1), [km/s]
%  - mu; gravitational parameter of body being orbited, [km^3 / s^2]

% if only 2 args get passed in, then it's (x, mu) and reassign things
if nargin == 2
    mu = v_;
    v_ = r_(4:6);
    r_ = r_(1:3);
end

% a_ is a vector, a is a magnitude
r_ = reshape(r_, 3, 1); v_ = reshape(v_, 3, 1); % ensure proper notation
r = norm(r_);
v = norm(v_);

h_ = cross(r_,v_);                              % angular momentum
h = norm(h_);
e_ = cross(v_, h_) / mu - r_ / r;               % eccentricity vector
e = norm(e_);
n_ = cross([0 0 1]', h_);                       % line of nodes
n = norm(n_);

a = 1 / (2 / r - v^2 / mu);                     % semimajor axis
i = acos(dot([0 0 1]', h_) / h);                % inclination

% if orbit is equatorial (RAAN ambiguity), change n to x unit vector
if n < 1e-8, n_ = [1 0 0]'; n = 1; end
% if orbit is circular (w ambiguity), change e to x unit vector
if e < 1e-8, e_ = [1 0 0]'; end

RAAN = acos(dot(n_, [1 0 0]') / n);             % right ascension
if dot(n_, [0 1 0]') < 0, RAAN = 2*pi - RAAN; end
w = acos(dot(e_, n_) / (norm(e_) * n));         % argument of periapsis
if dot(e_, [0 0 1]') < 0, w = 2*pi - w; end
f = acos(dot(e_, r_) / (norm(e_) * r));         % true anomaly
if dot(r_, v_) < 0, f = 2*pi - f; end
end
