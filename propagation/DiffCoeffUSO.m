function [s1,s2,s3] = DiffCoeffUSO()
%DIFFCOEFFUSO Returns the diffusion coefficients for an Ultra-Stable 
%Oscillator (USO), obtained from Small et al., 2022, and Zucca and
%Tavella, 2005.
%   Small et al., 2022:
%   https://www.doi.org/10.33012/2022.18221
%   Zucca and Tavella, 2005: 
%   https://www.doi.org/10.1109/TUFFC.2005.1406554
%
%   Outputs: (dims),[units]
%    - s1; (1x1),[?] diffusion coefficient of white noise
%    - s2; (1x1),[?] diffusion coefficient of random walk frequency noise
%    - s3; (1x1),[?] diffusion coefficient assoc. w/ aging rate (0)

coeff = [2.53e-23 4.22e-24/6 11*1.00e-38/120];
[s1, s2, s3] = HadamardToDiffCoeff(coeff);
end

