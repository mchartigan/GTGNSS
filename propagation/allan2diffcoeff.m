function coeff = allan2diffcoeff(B)
%ALLAN2DIFFCOEFF Converts Allan variance coefficients to diffusion
%coefficients.
%   The Allan variance formula is given in (Zucca and Tavella, 2005) as
%    s_y^2 (t) = s_1^2/t + s_2^2/3 * t + c_3^2 / 2 * t^2
%   or
%    s_y^2 (t) = B(1)/t + B(2) * t + B(3) * t^2
%
%   Inputs:
%    - B; coefficients of Allan variance formula above
%   Output:
%    - coeff; diffusion coefficients
arguments
    B (1,3) double {mustBeNonnegative}
end

s1 = sqrt(B(1));            % diffusion coefficient of white noise
s2 = sqrt(3 * B(2));        % ^ of random walk frequency noise
s3 = 0;                     % set to zero (aging rate is constant)
coeff = [s1 s2 s3];
end

