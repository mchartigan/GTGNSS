function [s1,s2,s3] = DiffCoeffRAFS()
%DIFFCOEFFRAFS Returns the diffusion coefficients for a Rubidium Atomic 
%Frequency Standard (RAFS), obtained from Small et al., 2022, and Zucca and
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

coeff = [3.70e-24 1.87e-33/6 11*7.56e-59/120];
[s1, s2, s3] = HadamardToDiffCoeff(coeff);
end

