function [T,E,df] = avar(t,x)
%AVAR Compute the traditional Allan variance of provided phase data
%   Input:
%    - t; time of phase measurements
%    - phase; measurements of phase offset (in s) of DUT vs reference
%   Output:
%    - T; time interval corresponding to E
%    - E; Allan variance for corresponding time interval
%    - df; degrees of freedom for estimator (ffm model)
arguments (Input)
    t   (1,:)   double
    x   (1,:)   double
end
arguments (Output)
    T   (1,:)   double
    E   (1,:)   double {mustBePositive}
    df  (1,:)   double
end

t = t - t(1);       % set start of interval to 0
N = length(t);      % number of samples
tau = t(2) - t(1);  % time between adjacent samples
% log-scale intervals to reduce computational load
nf = floor((N-1)/2);                    % max timestep
ints = 10.^(log10(tau):0.1:log10(nf));  % log-spaced intervals
n = unique(round(ints/tau)*tau);        % unique intervals at multiples of tau
E = zeros(1,length(n));

% wpm = (N+1)*(N-2*n)./(2*(N-n));
% fpm = exp(1./sqrt(log((N-1)./(2*n)).*log((2*n+1)*(N-1)/4)));
% wfm = (3*(N-1)./(2*n) - 2*(N-2)/N).*(4*n.^2./(4*n.^2+5));
% rwfm = (N-2)./n.*((N-1)^2-3*n*(N-1)+4*n.^2)/(N-3)^2;

% use flicker FM model to compute degrees of freedom since it averages all
% the effects (wpm and fpm don't play a part)
df = 5*N^2./(4*n.*(N+3*n));

% figure();
% plotformat("APA", 0.5);
% loglog(n, wpm);
% hold on;
% loglog(n, fpm);
% loglog(n, wfm);
% loglog(n, ffm);
% loglog(n, rwfm);
% hold off; grid on;
% legend(["wpm","fpm","wfm","ffm","rwfm"]);

% compute overlapping Allan variance
for i=1:length(n)
    k = 1:N-2*n(i);
    x_sum = sum((x(k+2*n(i)) - 2*x(k+n(i)) + x(k)).^2);
    E(i) = x_sum / (2*n(i)^2*tau^2*(N - 2*n(i)));
end

T = n*tau;          % time intervals
end