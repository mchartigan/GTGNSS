function [R,C,S,norms] = cofloader(fname,normalized)
%COFLOADER Loads nonspherical gravity coefficients from .cof files.
%   Input:
%    - fname; path to .cof file
%    - normalized; should coefficients be normalized? (default no)
%   Output:
%    - R; reference radius of body
%    - C; cosine coefficients of spherical harmonics
%    - S; sine coefficients of spherical harmonics
%    - norms; normalization factors for coefficients
arguments
    fname       (1,:) {mustBeText}
    normalized  (1,1) = false
end

lines = readlines(fname);               % get lines in file

i = 1;                                  % get row where numbers start
while ~startsWith(lines(i,1), "POTFIELD"), i = i + 1; end

% all data in these .cof files are fixed-width, so it should be accessed by
% character indexing.

% First, the row starting with POTFIELD contains the max degree and order;
% it also contains the reference radius.
meta = lines{i};
n_max = str2double(meta(9:11));         % max degree
R = str2double(meta(40:59));            % reference radius
lines = lines(i+1:end);

C = zeros(n_max+1, n_max+1);
C(1,1) = 1;                             % value not provided in .cof files
S = C;
norms = ones(n_max+1, n_max+1);         % normalization factors

i = 1;
while ~startsWith(lines{i}, "END")
    data = lines{i};
    n = str2double(data(9:11));         % extract degree
    m = str2double(data(12:14));        % extract order
    C(n+1,m+1) = str2double(data(18:38));
    % if order is zero the row will be empty so don't try to grab it
    if m ~= 0, S(n+1,m+1) = str2double(data(39:59)); end

    i = i + 1;
end

% if coefficients are requested to be denormalized, do that. Otherwise,
% return as-is. Also enter this loop if the normalization factors were
% requested.
if ~normalized || nargout > 3
    for n = 0:n_max
        for m = 0:n
            if m == 0, k = 1; else, k = 2; end

            % assign normalization factor to matrix
            norms(n+1,m+1) = sqrt(factorial(n+m)/(k*(2*n+1)*factorial(n-m)));

            % if coefficients shouldn't be normalized, divide by factor
            if ~normalized
                C(n+1,m+1) = C(n+1,m+1) / norms(n+1,m+1);
                S(n+1,m+1) = S(n+1,m+1) / norms(n+1,m+1);
            end
        end
    end
end
end

