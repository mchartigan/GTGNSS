function f = fast_harmonics(x, nn, mu, R, C, S)
%F_HARMONICS compute gravitational acceleration from spherical harmonics 
%of body
%   Input:
%    - x; position vector of point (3,) [km]
%    - nn; max degree (and order m) of harmonics
%    - mu; gravitational parameter of body [km^3/s^2]
%    - R; reference radius of body [km]
%    - C; C coefficients of spherical harmonics (n,n)
%    - S; S coefficients of spherical harmonics (n,n)

    [V, W] = VandW(x, nn + 1, R);

    n = 1:nn+1; m = 1:nn+1;
    A = -mu/R^2;
    f= A * [C(1:nn+1,1)' * V(2:nn+2,2)
            C(1:nn+1,1)' * W(2:nn+2,2)
            trace((n'-m+1) * (C(n,m).*V(n+1,m) + S(n,m).*W(n+1,m))')];

    m = m(2:end);
    fact = factorial(abs(n'-m+2))./factorial(abs(n'-m));
    f(1) = f(1) + A/2 * trace((C(n,m)*V(n+1,m+1)' + S(n,m)*W(n+1,m+1)') ...
           - fact * (C(n,m).*V(n+1,m-1) + S(n,m).*W(n+1,m-1))');
    f(2) = f(2) + A/2 * trace((C(n,m)*W(n+1,m+1)' - S(n,m)*V(n+1,m+1)') ...
           + fact * (C(n,m).*W(n+1,m-1) - S(n,m).*V(n+1,m-1))');
end

function [V, W] = VandW(x, nn, R)
%VANDW Computes recursive Legendre polynomials for spherical harmonics
%   Input
%    - x; position vector of point (3,) [km]
%    - nn; max degree and order of harmonics
%    - R; reference radius of body [km]

    r = norm(x);
    V = zeros(nn+1,nn+1); W = V;
    V(1,1) = R / r; W(1,1) = 0;         % initial conditions
    A = R / r^2;

    for m=0:nn          % iterate over m cols (add 1 for indexing)
        if m ~= 0       % if not V_00, W_00 (already defined)
            V(m+1,m+1) = (2*m-1)*A*(x(1)*V(m,m) - x(2)*W(m,m));
            W(m+1,m+1) = (2*m-1)*A*(x(1)*W(m,m) + x(2)*V(m,m));
        end

        for n=m+1:nn    % iterate over n rows (add 1 for indexing)  
            V(n+1,m+1) = (2*n-1)/(n-m)*x(3)*A*V(n,m+1);
            W(n+1,m+1) = (2*n-1)/(n-m)*x(3)*A*W(n,m+1);
            if n ~= m+1
                V(n+1,m+1) = V(n+1,m+1) - (n+m-1)/(n-m)*A*R * V(n-1,m+1);
                W(n+1,m+1) = W(n+1,m+1) - (n+m-1)/(n-m)*A*R * W(n-1,m+1);
            end
        end
    end
end
