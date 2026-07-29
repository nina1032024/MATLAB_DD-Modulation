function gs = gen_gs_frac(p, g, l, k, Lrange)
% =========================================================================
% Program : gen_gs_frac.m
% Description :
%   Generalization of gen_gs.m (Stage 2) to FRACTIONAL delay taps l_i,
%   chap4_list.md Sec. A.3 (4.6): the exact delta function delta[l-l_i]
%   used for on-grid taps (4.8) is replaced by a sinc INTERPOLATION kernel,
%   spreading each path's energy across a range of integer delay indices l:
%
%       g^s[l,q] = sum_i g_i * z^{k_i(q-l)} * sinc(l - l_i),   l = 0...Lrange-1
%
%   This reduces EXACTLY to gen_gs.m's delta-based (4.8) when every l_i is
%   an integer (sinc(integer)=0 except sinc(0)=1) — verified as a
%   regression check in Sim-15. Lrange must extend well past the nominal
%   l_max, or the sinc tail gets truncated and leaked energy is lost
%   (chap4_list.md Sec. H point 6 / Sim-15: "l 的範圍要擴大到 l_max 以外
%   (例如 0…l_max+8)").
%
% Input  : p      - parameter struct (needs p.z, p.NM)
%          g      - 1xP complex path gains
%          l      - 1xP delay taps (may be fractional, >= 0)
%          k      - 1xP Doppler taps (signed, integer or fractional)
%          Lrange - number of delay-index rows to compute, l = 0...Lrange-1
%                   (default: ceil(max(l)) + 9, i.e. nominal l_max + 8)
% Output : gs     - Lrange x NM matrix, gs(row=l+1, col=q+1) = g^s[l,q]
% =========================================================================

if nargin < 5 || isempty(Lrange)
    Lrange = ceil(max(l)) + 9;
end

NM = p.NM;
q  = 0:NM-1;
gs = zeros(Lrange, NM);

for i = 1:numel(g)
    li = l(i);
    for ll = 0:Lrange-1
        w = sinc_fn(ll - li);
        if abs(w) < 1e-15
            continue;
        end
        phase = p.z .^ (k(i) * (q - ll));
        gs(ll+1, :) = gs(ll+1, :) + g(i) * w * phase;
    end
end

end
