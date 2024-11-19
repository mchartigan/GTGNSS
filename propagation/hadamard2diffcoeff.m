function coeff = hadamard2diffcoeff(B)
%HADAMARD2DIFFCOEFF Converts Hadamard variance coefficients to diffusion
%coefficients.
%   The Hadamard variance formula is given in (Zucca and Tavella, 2005) as
%    s_y^2 (t) = s_1^2/t + s_2^2/6 * t + 11*s_3^2/120 * t^3 + mu_3^2/6 t^4
%   or
%    s_y^2 (t) = B(1)/t + B(2) * t + B(3) * t^3 + B(4) t^4
%
%   Inputs:
%    - B; coefficients of Allan variance formula above
%   Output:
%    - coeff; diffusion coefficients
arguments
    B (1,3) double {mustBeNonnegative}
end

s1 = sqrt(B(1));                % diffusion coefficient of white noise
s2 = sqrt(B(2) * 6);            % ^ of random walk frequency noise
s3 = sqrt(120 * B(3) / 11);     % set to zero (aging rate is constant)
coeff = [s1 s2 s3];
end

