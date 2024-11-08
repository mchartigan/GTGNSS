function [s1,s2,s3] = HadamardToDiffCoeff(B)
%HADAMARDTODIFFCOEFF Converts Hadamard variance coefficients to diffusion
%coefficients.
%   The Hadamard variance formula is given in (Zucca and Tavella, 2005) as
%    s_y^2 (t) = s_1^2/t + s_2^2/6 * t + 11*s_3^2/120 * t^3 + mu_3^2/6 t^4
%   or
%    s_y^2 (t) = B(1)/t + B(2) * t + B(3) * t^3 + B(4) t^4
%
%   Inputs: (dims),[units]
%    - B; (1x3),[?] coefficients of Allan variance formula above

s1 = sqrt(B(1));                % diffusion coefficient of white noise
s2 = sqrt(B(2) * 6);            % ^ of random walk frequency noise
s3 = sqrt(120 * B(3) / 11);     % set to zero (aging rate is constant)
end

