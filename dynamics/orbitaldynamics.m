function dxdt = orbitaldynamics(t,x,pri,N,bodies)
%ORBITALDYNAMICS Inertial dynamics for a satellite orbiting a primary
%body with gravitational influence from a secondary (and tertiary);
%includes nonspherical gravity.
%   Input:
%    - t; simulation time
%    - x; state [pos (km); vel (km/s)] of spacecraft
%    - pri; struct, {GM: gravitational parameter, x: @(t) position [km]}
%           for primary body
%    - N; maximum degree and order of harmonics to compute
%    - sec; array of struct, {GM: gravitational parameter, x: @(t) position [km]}
%           for secondary bodies
arguments
    t       (1,1) double
    x       (6,1) double
    pri     (1,1) struct
    N       (1,1) {mustBePositive, mustBeInteger}
    bodies  (1,:) struct = []
end

dxdt = zeros(6,1);
x_1s = x(1:3);
r_1s = norm(x_1s);

% get lunar nonspherical gravity effects
T = cspice_pxform('J2000', pri.frame, t);
x_me = T * x_1s;
f_ns = fast_harmonics(x_me, N, pri.GM, pri.R, pri.C, pri.S);
f_ns = T' * f_ns;

f_sec = zeros(3,1);
% get influence of secondary bodies
for i=1:length(bodies)
    x_13 = bodies(i).x(t);
    r_13 = norm(x_13);
    x_s3 = x_13 - x_1s;
    r_s3 = norm(x_s3);
    % Formulation according to Geyling and Westerman for numerical improvements
    % %    Fundamentals of Astrodynamics, Vallado, pp.575 8-36
    % f_sec = f_sec - bodies(i).GM * (x_1s - 3*x_13*(x_1s'*x_13)/r_13^2 - 15/2*((x_1s'*x_13)/r_13^2)^2*x_13)/r_13^3;
    % Standard formulation
    f_sec = f_sec + bodies(i).GM * (x_s3 / r_s3^3 - x_13/r_13^3);
    % Formulation according to Roy for numerical improvements
    % Q = (r_1s^2 + 2*x_1s'*x_s3) * (r_13^2 + r_13*r_s3 + r_s3^2) / (r_13^3 * r_s3^3 * (r_13 + r_s3));
    % f_sec = f_sec + bodies(i).GM * (x_s3 * Q - x_1s/r_s3^3);
end 

% compute drag
f_d = 0;
% compute SRP
f_SRP = 0;

dxdt(1:3) = x(4:6);
dxdt(4:6) = f_ns + f_sec + f_d + f_SRP;
end