function phi = chebyshev(n,x,type)
%CHEBYSHEV Provides Chebyshev polynomials of the [type] kind of order n at
%point(s) x. Overrides native MATLAB solution since it uses symbolic
%variables and is trash for this purpose.
%   Input:
%    - n; degree(s) of polynomial, increasing
%    - x; point(s) to evaluate at
%    - type; optional, default 1 (type T). Could be 2 (type U).
%   Output:
%    - phi (length(x), length(n); coefficient matrix

if nargin < 3, type = 1; end

nmax = n(end)+1;
m = length(x);
x = reshape(x,m,1);
temp = zeros(m, nmax);

for i=1:nmax
    if i == 1
        temp(:,i) = ones(m,1);
    elseif i == 2
        temp(:,i) = type * x;
    else
        temp(:,i) = 2*x.*temp(:,i-1) - temp(:,i-2);
    end
end

phi = temp(:,n+1);
end

