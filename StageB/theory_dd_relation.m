function Y = theory_dd_relation(p, X, g, l, k)
% =========================================================================
% Program : theory_dd_relation.m
% Description :
%   Pure for-loop "theoretical value generator" for the symbol-level
%   delay-Doppler input-output relation under a RECTANGULAR pulse, integer
%   (on-grid) taps, and the RCP channel model (chap4_list.md Sim-10,
%   Sec. E.2, eqs. (4.118)-(4.122)):
%
%       Y[m,n] = sum_i g_i * z^{k_i(m-l_i)}
%                       * X[[m-l_i]_M, [n-k_i]_N]                     , m >= l_i
%       Y[m,n] = sum_i g_i * z^{k_i([m-l_i]_M)} * e^{-j2*pi*n/N}
%                       * X[[m-l_i]_M, [n-k_i]_N]                     , m <  l_i
%
%   Derivation sketch (also in the Ch.4 progress report): under RCP, a
%   sample r[q] with q = m+nM and m < l_i must borrow its input sample from
%   time-block (n-1) instead of block n (mod N, because of the RCP
%   frame-periodic wraparound) — that one-block time shift is exactly what
%   turns into the extra e^{-j2*pi*n/N} phase after the DFT-in-n
%   demodulation step (a discrete-time shift becomes a linear phase ramp in
%   the DFT domain). For m >= l_i no such borrowing happens and the
%   relation is the plain "shift + phase" form of the first branch.
%
%   *** k_i is used RAW (signed) in the z^{...} phase, but WRAPPED
%   mod N when it indexes into X (delay index is wrapped mod M
%   analogously) — see gen_channel_taps.m. ***
%
% Input  : p - parameter struct (needs p.M, p.N, p.z)
%          X - M x N DD-domain transmit symbol matrix
%          g - 1xP complex path gains
%          l - 1xP delay taps (integers)
%          k - 1xP Doppler taps (signed integers)
% Output : Y - M x N DD-domain theoretical receive symbol matrix
% =========================================================================

M = p.M;  N = p.N;  z = p.z;
Y = zeros(M, N);

for m = 0:M-1
    for n = 0:N-1
        acc = 0;
        for i = 1:numel(g)
            li = l(i);  ki = k(i);
            mp = mod(m - li, M);
            np = mod(n - ki, N);
            if m >= li
                acc = acc + g(i) * z^(ki*(m-li)) * X(mp+1, np+1);
            else
                acc = acc + g(i) * z^(ki*mp) * exp(-1j*2*pi*n/N) * X(mp+1, np+1);
            end
        end
        Y(m+1, n+1) = acc;
    end
end

end
