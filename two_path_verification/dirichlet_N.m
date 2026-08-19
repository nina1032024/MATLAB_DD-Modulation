function D = dirichlet_N(x, N)
%DIRICHLET_N  D_N(x) = sum_{p=0}^{N-1} exp(j*2*pi*x*p/N)
%
%   Closed form: D_N(x) = sin(pi*x)/sin(pi*x/N) * exp(j*pi*x*(N-1)/N),
%   which holds for every real x (finite geometric series, no
%   approximation). At x integer the closed form is 0/0; the exact
%   limiting value is D_N(x) = N when x is a multiple of N, and 0
%   otherwise -- note the "mod N", not a bare x==0 test, since e.g.
%   x = -2 with N = 8 is also on-grid.

D = zeros(size(x));
isInt = abs(x - round(x)) < 1e-12;
D(isInt)  = N * (mod(round(x(isInt)), N) == 0);
xi = x(~isInt);
D(~isInt) = sin(pi*xi)./sin(pi*xi/N) .* exp(1j*pi*xi*(N-1)/N);
end
