function [s1,s2,s3] = DiffCoeffCSAC()
%DIFFCOEFFCSAC Returns the diffusion coefficients for a Microsemi Space
%CSAC, obtained from information in their datasheet and Zucca and Tavella's
%2005 paper.
%   Microsemi Space CSAC datasheet:
%   https://www.microsemi.com/document-portal/doc_download/1243238-space-csac-datasheet
%   Zucca and Tavella, 2005: 
%   https://www.doi.org/10.1109/TUFFC.2005.1406554
%
%   Outputs: (dims),[units]
%    - s1; (1x1),[?] diffusion coefficient of white noise
%    - s2; (1x1),[?] diffusion coefficient of random walk frequency noise
%    - s3; (1x1),[?] diffusion coefficient assoc. w/ aging rate (0)

taus = [1 10 100 1000]';            % intervals for allan deviation
stds = [3e-10 1e-10 3e-11 1e-11]';  % deviations for different intervals
a_hi = 9e-10 / (86400 * 30);        % Hz/Hz/s, upper end of aging rate a
a_lo = 6e-10 / (86400 * 30);        % Hz/Hz/s, lower end of aging rate a

% Assume frequency drift a is constant, c_3 = a. Zucca and Tavella, 2005, 
% give the Allan variance formula as
%    s_y^2 (t) = s_1^2/t + s_2^2/3 * t + c_3^2 / 2 * t^2
% Since c_3 is known, the formula can become
%    s_y^2 (t) - c_3^2 / 2 * t^2 = s_1^2/t + s_2^2/3 * t
% We can then do a polynomial fit to this with the provided data points in
% the MicroSemi CSAC datasheet. The solution to Ax = b takes
% the form (A' * A)^(-1) * A' * b.

A = [1./taus taus];
b = stds.^2 - a_hi^2 / 2 * taus.^2;
coeff = [(A' * A) \ A' * b; a_hi^2 / 2];

% % plot model to confirm it's working
% t = linspace(taus(1), taus(end), 1000)';
% mdl = sqrt([1./t t t.^2] * coeff);
% 
% h1 = figure();
% fsz = 12;
% loglog(t, mdl, 'LineWidth', 1.5);
% hold on;
% scatter(taus, stds, 'rx');
% hold off; grid on;
% xlabel('\tau (s)', 'FontSize', fsz);
% ylabel('\sigma_y', 'FontSize', fsz);
% title('Microsemi Space CSAC Allan Deviation Fit', 'FontSize', fsz);
% legend(["Model fit", "From datasheet"], "location", "northeast", ...
%     'FontSize', fsz);
% fontname(h1, "Times New Roman");

[s1, s2, s3] = AllanToDiffCoeff(coeff);
end

