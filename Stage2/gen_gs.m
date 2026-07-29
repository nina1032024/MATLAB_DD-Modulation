function gs = gen_gs(p, g, l, k)
% =========================================================================
% Program : gen_gs.m
% Description :
%   Build the integer-tap time-domain delay-Doppler spreading function
%   g^s[l,q]  (chap4_list.md Sec. A.3, Appendix-C code 8):
%
%       g^s[l,q] = sum_{i : l_i = l}  g_i * z^{k_i (q-l)},   z = e^{j2pi/(NM)}
%                                                                       (4.8)
%
%   Stored as an (l_max+1) x NM matrix, row index l = 0...l_max,
%   column index q = 0...NM-1.  Paths that happen to share the same delay
%   tap l_i are summed together (rare with the default taps, but handled
%   correctly in general).
%
% Input  : p - parameter struct (needs p.z, p.NM)
%          g - 1xP complex path gains
%          l - 1xP delay taps (integers, 0-indexed, l_i < M)
%          k - 1xP Doppler taps (signed integers, RAW — not mod N, see
%              gen_channel_taps.m)
% Output : gs - (l_max+1) x NM matrix, gs(row=l+1, col=q+1) = g^s[l,q]
% =========================================================================

l_max = max(l);
NM = p.NM;
q  = 0:NM-1;                        % 1 x NM

gs = zeros(l_max+1, NM);
for i = 1:numel(g)
    li = l(i);
    phase = p.z .^ (k(i) * (q - li));   % 1 x NM
    gs(li+1, :) = gs(li+1, :) + g(i) * phase;
end

end
