function B = stability2allan(taus,stds,a)
%STABILITY2ALLAN Takes in short-term stability data points and fits them
%to the Allan deviation model, providing the model coefficients. These
%coefficients must be passed through AllanToDiffCoeff() to get the
%resultant diffusion coefficients, though.
%   Input:
%    - taus; time spans, (1xn)
%    - stds; corresponding deviations, in s/s (1xn)
%    - a; aging rate of oscillator, in s/s^2
arguments
    taus (:,1) double {mustBePositive}
    stds (:,1) double {mustBePositive}
    a    (1,1) double
end
nrm = stds(end)^2;
b = (stds.^2 - a^2 / 2 * taus.^2) / nrm;

opts = fitoptions('Method', 'LinearLeastSquares', ...
    'Lower', [0 0], 'Robust', 'off');
type = fittype({'1/x', 'x'}, 'options', opts);
[curve, gof] = fit(taus, b, type);
B = [coeffvalues(curve) * nrm a^2/2]';
end

