function B = StabilityToHadamard(taus,stds)
%STABILITYTOHADAMARD Takes in short-term stability data points and fits 
%them to the Hadamard deviation model, providing the model coefficients.
%These coefficients must be passed through HadamardToDiffCoeff() to get 
%the resultant diffusion coefficients, though.
%   Input:
%    - taus; time spans, (1xn)
%    - stds; corresponding deviations, in s/s (1xn)
arguments
    taus (:,1) double {mustBePositive}
    stds (:,1) double {mustBePositive}
end
nrm = stds(end)^2;
b = stds.^2 / nrm;

opts = fitoptions('Method', 'LinearLeastSquares', ...
    'Lower', [0 0 0], 'Robust', 'off');
type = fittype({'1/x', 'x', 'x^3'}, 'options', opts);
[curve, ~] = fit(taus, b, type);
B = coeffvalues(curve)' * nrm;
end

