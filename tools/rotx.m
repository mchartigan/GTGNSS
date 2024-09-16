function R = rotx(a)
%ROTX Returns a rotation matrix about the x-axis by a radians (given in the
%passive notation).
%   Input: (dims), [units]
%    - a; (1x1),[rad] angle to rotate about
%   Output:
%    - R; (3x3),[N/A] x rotation matrix corresponding to a angle

R = [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)];
end

