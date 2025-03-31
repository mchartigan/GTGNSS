function [dxdt] = halodyn(X, mu)
%HALODYN Normalized dynamics of the circular restricted 3-body problem 
%(CRTBP).
%   Inputs: (dims),[units]
%    - X ; (6x1),[D,D/T] state vector
%    - mu; (1x1),[N/A] normalized m2

% sort state vector
x  = X(1); y  = X(2); z  = X(3);    % position
vx = X(4); vy = X(5); vz = X(6);    % velocity
phi = reshape(X(7:end), 6, 6);      % STM

% compute radial distances
r1 = sqrt((x + mu)^2     + y^2 + z^2);  % from Earth to point
r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);  % from moon to point

% compute acceleration
ax =  2*vy + x - (1 - mu) * (x + mu)/r1^3 - mu * (x - 1 + mu)/r2^3;
ay = -2*vx + y - (1 - mu) * y/r1^3 - mu * y/r2^3;
az =           - (1 - mu) * z/r1^3 - mu * z/r2^3;

% compute jacobian F of [dU/dx dU/dy dU/dz], where U is the potential fn
% See https://adsabs.harvard.edu/full/1984CeMec..32...53H (Howell, 1983)
dx2  = 1 - (1-mu)*(r1^2 - 3*(x+mu)^2)/r1^5 - mu*(r2^2 - 3*(x-1+mu)^2)/r2^5;
dy2  = 1 - (1-mu)*(r1^2 - 3*y^2)     /r1^5 - mu*(r2^2 - 3*y^2)       /r2^5;
dz2  =   - (1-mu)*(r1^2 - 3*z^2)     /r1^5 - mu*(r2^2 - 3*z^2)       /r2^5;
dxdy = (1-mu)*y*3*(x+mu)/r1^5 + mu*y*3*(x-1+mu)/r2^5;
dxdz = (1-mu)*z*3*(x+mu)/r1^5 + mu*z*3*(x-1+mu)/r2^5;
dydz = (1-mu)*z*3*y     /r1^5 + mu*z*3*y       /r2^5;
Uxx = [dx2  dxdy dxdz
       dxdy dy2  dydz
       dxdz dydz dz2 ];
F = [zeros(3) eye(3); Uxx [0 2 0; -2 0 0; 0 0 0]];

dphi = F * phi;                                     % derivative of STM
dxdt = [vx vy vz ax ay az reshape(dphi, 1, 36)]';   % state derivative

end

