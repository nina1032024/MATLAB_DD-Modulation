function Y = ideal_pulse_relation(p, X, g, l, k)
% =========================================================================
% Program : ideal_pulse_relation.m
% Description :
%   Symbol-level DD I/O relation for the IDEAL (bi-orthogonal) pulse,
%   chap4_list.md Sim-11, eq. (4.15): a pure, TIME-INVARIANT 2D circular
%   convolution with NO extra phase term and no m<l_i boundary correction
%   (unlike the rectangular-pulse case in theory_dd_relation.m):
%
%       Y_ideal[m,n] = sum_i g_i * X[[m-l_i]_M, [n-k_i]_N]
%
%   The gap between this and the rectangular-pulse result (theory_dd_relation.m)
%   is exactly the "loss of bi-orthogonality" effect discussed in
%   chap4_list.md Sec. C.8: the rectangular pulse's Heisenberg/Wigner pair is
%   not perfectly bi-orthogonal, which shows up as the extra phase
%   z^{k_i(m-l_i)} and, for m < l_i, the extra e^{-j2*pi*n/N} boundary term.
%
% Input  : p - parameter struct (needs p.M, p.N)
%          X - M x N DD-domain transmit symbol matrix
%          g - 1xP complex path gains
%          l - 1xP delay taps (integers)
%          k - 1xP Doppler taps (signed integers, only used mod N here)
% Output : Y - M x N DD-domain ideal-pulse receive symbol matrix
% =========================================================================

M = p.M;  N = p.N;
Y = zeros(M, N);

for m = 0:M-1
    for n = 0:N-1
        acc = 0;
        for i = 1:numel(g)
            mp = mod(m - l(i), M);
            np = mod(n - k(i), N);
            acc = acc + g(i) * X(mp+1, np+1);
        end
        Y(m+1, n+1) = acc;
    end
end

end
