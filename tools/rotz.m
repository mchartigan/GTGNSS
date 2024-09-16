function R = rotz(a)
%ROTZ Returns a rotation matrix about the z-axis by a radians (given in the
%passive notation).
%   Input: (dims), [units]
%    - a; (1x1),[rad] angle to rotate about
%   Output:
%    - R; (3x3),[N/A] z rotation matrix corresponding to a angle

R = [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];
end

