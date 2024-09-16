function [s1,s2,s3] = AllanToDiffCoeff(B)
%ALLANTODIFFCOEFF Converts Allan variance coefficients to diffusion
%coefficients.
%   The Allan variance formula is given in (Zucca and Tavella, 2005) as
%    s_y^2 (t) = s_1^2/t + s_2^2/3 * t + c_3^2 / 2 * t^2
%   or
%    s_y^2 (t) = B(1)/t + B(2) * t + B(3) * t^2
%
%   Inputs: (dims),[units]
%    - B; (1x3),[?] coefficients of Allan variance formula above
%   Output:
%    - s; (1x3),[?] diffusion coefficients

s1 = sqrt(B(1));            % diffusion coefficient of white noise
s2 = sqrt(3 * B(2));        % ^ of random walk frequency noise
s3 = 0;                     % set to zero (aging rate is constant)
end

