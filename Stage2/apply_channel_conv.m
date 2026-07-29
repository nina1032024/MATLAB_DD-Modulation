function r = apply_channel_conv(p, gs, s, variant)
% =========================================================================
% Program : apply_channel_conv.m
% Description :
%   Direct time-varying convolution implementation of (4.36),
%   COMPLETELY INDEPENDENT of the matrix form G (gen_G.m), so that the two
%   can be cross-checked against each other (chap4_list.md Sim-5):
%
%       r[q] = sum_{l=0}^{l_max} g^s[l,q] * s[q-l] ,   q = 0...NM-1
%
%   with the same RZP / RCP boundary convention as gen_G.m:
%     'RZP' : s[q-l] treated as 0 when q-l < 0  (zero padding)
%     'RCP' : s[q-l] read circularly, index mod(q-l, NM)
%
%   Deliberately written as a plain double for-loop (no matrix products)
%   per chap4_list.md's instruction "寫一個 for-loop 的時變摺積函式" — the
%   whole point of Sim-5 is that this and G*s are independent derivations
%   of the same equation.
%
% Input  : p       - parameter struct (needs p.NM)
%          gs      - (l_max+1) x NM spreading matrix from gen_gs.m
%          s       - NM x 1 time-domain transmit vector
%          variant - 'RZP' (default) | 'RCP'
% Output : r       - NM x 1 time-domain receive vector
% =========================================================================

if nargin < 4 || isempty(variant), variant = 'RZP'; end

NM = p.NM;
l_max = size(gs, 1) - 1;
s = s(:);
r = zeros(NM, 1);

for q = 0:NM-1
    acc = 0;
    for l = 0:l_max
        ql = q - l;
        switch upper(variant)
            case 'RZP'
                if ql >= 0
                    acc = acc + gs(l+1, q+1) * s(ql+1);
                end
            case 'RCP'
                qlw = mod(ql, NM);
                acc = acc + gs(l+1, q+1) * s(qlw+1);
            otherwise
                error('apply_channel_conv:variant', 'Unknown variant ''%s''. Use RZP | RCP.', variant);
        end
    end
    r(q+1) = acc;
end

end
