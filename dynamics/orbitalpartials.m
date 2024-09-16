function A = orbitalpartials(t,x,pri,sec)
%ORBITALPARTIALS Returns the Jacobian of orbitaldynamics w.r.t. x 
%(spherical harmonics approximated to J2).
%   Input:
%    - t; time, seconds past J2000
%    - x_; state [pos (km); vel (km/s)] of spacecraft
%    - pri; struct, {GM: gravitational parameter, x: @(t) position [km]}
%           for primary body
%    - sec; optional struct, {GM: gravitational parameter, x: @(t) position [km]}
%           for secondary body
arguments
    t   (1,1) double
    x  (6,1) double
    pri (1,:) struct
    sec (1,:) struct = []
end

x_1s = x(1:3);
r_1s = norm(x_1s);
J2 = -pri.C(3,1);
A = zeros(6,6);

A(1:3,4:6) = eye(3);
% point mass gravity gradient
A(4:6,1:3) = pri.GM * (3*(x_1s*x_1s')/r_1s^5 - 1/r_1s^3 * eye(3));

% third-body perturbations
for i=1:length(sec)
    x_si = sec(i).x(t) - x_1s;
    r_si = norm(x_si);
    A(4:6,1:3) = A(4:6,1:3) + sec(i).GM * (3*(x_si*x_si')/r_si^5 - 1/r_si^3 * eye(3));
end

% J2 gravity gradient
% x is in J2000, transform to Moon ME to apply this partial
T = cspice_pxform('MOON_ME', 'J2000', t);
x_1s = T' * x_1s;
z = x_1s(3);
S = diag([1 1 3]);
Q = [1 1 3; 1 1 3; 3 3 3];
V1 = (5*z^2-r_1s^2)/r_1s^7 * S - 35*z^2*(x_1s*x_1s')/r_1s^9 + Q .* (5*(x_1s*x_1s')/r_1s^7);
% transform the partial back to J2000
A(4:6,1:3) = A(4:6,1:3) + T * 3*pri.GM*pri.R^2*J2/2 * V1 * T';
end

