function D = periodic_sinc(x, N)
% =========================================================================
% Program : periodic_sinc.m
% Description :
%   Periodic sinc (Dirichlet) kernel used to describe fractional-Doppler
%   leakage in the delay-Doppler domain (chap4_list.md Sec. C.9, eq. (4.79)):
%
%       D_N(x) = sin(pi*x) / sin(pi*x/N) * exp(j*pi*x*(N-1)/N)
%
%   Derivation (used in Sim-14, matches theory_dd_relation.m's style of
%   re-deriving formulas from first principles rather than trusting the
%   compressed summary blindly): for a single path with integer delay l and
%   possibly FRACTIONAL Doppler tap k_i,
%       nu_tilde_{m,l}[n] = g_i * z^{k_i(m-l)} * e^{j*2*pi*k_i*n/N}
%   so its forward DFT (4.77) is a geometric sum
%       nu_{m,l}[k] = g_i*z^{k_i(m-l)} * (1/N) * sum_{n=0}^{N-1} e^{j2pi(k_i-k)n/N}
%   and the closed form of that geometric sum is exactly (1/N)*D_N(k_i-k)
%   using the standard Dirichlet-kernel identity. When k_i is an integer,
%   D_N(k_i-k) collapses to N*delta[k-k_i] (see the removable-singularity
%   branch below), recovering the plain on-grid result of Stage 2.
%
%   Removable singularity: D_N(x) has 0/0 at every integer multiple of N
%   (both sin(pi*x) and sin(pi*x/N) vanish there). Taking the limit
%   (L'Hopital) gives exactly D_N(mN) = N for every integer m — implemented
%   directly below rather than relying on floating-point cancellation.
%
% Input  : x - array of any size (real), typically (k_i - k) for k=0...N-1
%          N - Doppler grid size
% Output : D - D_N(x), same size as x (complex in general)
% =========================================================================

sz = size(x);
x = x(:);
tol = 1e-9;
xmodN = mod(x, N);
isSing = (xmodN < tol) | (xmodN > N - tol);

D = zeros(size(x));
xs = x(~isSing);
D(~isSing) = sin(pi*xs) ./ sin(pi*xs/N) .* exp(1j*pi*xs*(N-1)/N);
D(isSing) = N;

D = reshape(D, sz);

end
