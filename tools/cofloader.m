function [R,C,S] = cofloader(fname)
%COFLOADER Loads nonspherical gravity coefficients from .cof files.
%   Input:
%    - filename; path to .cof file
%   Output:
%    - 
arguments
    fname (1,:) {mustBeText}
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

C_norm = zeros(n_max+1, n_max+1);
C_norm(1,1) = 1;                        % value not provided in .cof files
S_norm = C_norm; C = C_norm; S = S_norm;

i = 1;
while ~startsWith(lines{i}, "END")
    data = lines{i};
    n = str2double(data(9:11));         % extract degree
    m = str2double(data(12:14));        % extract order
    C_norm(n+1,m+1) = str2double(data(18:38));
    % if order is zero the row will be empty so don't try to grab it
    if m ~= 0, S_norm(n+1,m+1) = str2double(data(39:59)); end

    i = i + 1;
end

for n = 0:n_max
    for m = 0:n
        if m == 0
            k = 1;
        else
            k = 2;
        end
        C(n+1,m+1) = C_norm(n+1,m+1) / sqrt(factorial(n+m)/(k*(2*n+1)*factorial(n-m)));
        S(n+1,m+1) = S_norm(n+1,m+1) / sqrt(factorial(n+m)/(k*(2*n+1)*factorial(n-m)));
    end
end
end

