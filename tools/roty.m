function R = roty(a)
%ROTY Returns a rotation matrix about the y-axis by a radians (given in the
%passive notation).
%   Input: (dims), [units]
%    - a; (1x1),[rad] angle to rotate about
%   Output:
%    - R; (3x3),[N/A] y rotation matrix corresponding to a angle

R = [cos(a) 0 -sin(a); 0 1 0; sin(a) 0 cos(a)];
end

