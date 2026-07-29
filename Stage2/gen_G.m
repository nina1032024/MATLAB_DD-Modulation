function G = gen_G(p, gs, variant)
% =========================================================================
% Program : gen_G.m
% Description :
%   Build the NM x NM time-domain channel matrix G from the TDL spreading
%   function g^s[l,q]  (chap4_list.md Sec. A.3 (4.38), Appendix-C code 8/15/16):
%
%       G[q, q-l] = g^s[l,q],     l = 0...l_max
%
%   Two variants of how the "q-l < 0" boundary (the first l_max samples of
%   the frame, which would need samples from before the frame) is handled:
%
%   'RZP' (Reduced Zero Padding) - DEFAULT for Sec. 4.4
%       q-l < 0 entries are simply left at zero (no wraparound).  This is
%       the "single isolated frame" assumption of Sec. 4.4 itself
%       (chap4_list.md Sec. H point 7) -> G is purely banded/lower-triangular
%       with l_max+1 sub-diagonals (reproduces Fig. 4.8).
%       Used for Sim-5 ... Sim-9.
%
%   'RCP' (Reduced Cyclic Prefix)
%       q-l < 0 entries wrap around using mod(q-l, NM), landing in the
%       upper-right corner of G. This is the idealization needed to derive
%       the clean symbol-level closed forms (4.118)-(4.122), and is the
%       variant used for Sim-10 (per chap4_list.md Sec. E.1, code 16 -> Sim-10).
%
% Input  : p       - parameter struct (needs p.NM)
%          gs      - (l_max+1) x NM spreading matrix from gen_gs.m
%          variant - 'RZP' (default) | 'RCP'
% Output : G       - NM x NM time-domain channel matrix
% =========================================================================

if nargin < 3 || isempty(variant), variant = 'RZP'; end

NM = p.NM;
l_max = size(gs, 1) - 1;
G = zeros(NM, NM);
q = (0:NM-1).';

for l = 0:l_max
    qml = q - l;                          % NM x 1
    switch upper(variant)
        case 'RZP'
            valid = qml >= 0;
            rows  = q(valid) + 1;
            cols  = qml(valid) + 1;
            vals  = gs(l+1, q(valid)+1).';
            idx   = sub2ind([NM, NM], rows, cols);
            G(idx) = G(idx) + vals;

        case 'RCP'
            colsW = mod(qml, NM) + 1;
            rows  = q + 1;
            vals  = gs(l+1, :).';
            idx   = sub2ind([NM, NM], rows, colsW);
            G(idx) = G(idx) + vals;

        otherwise
            error('gen_G:variant', 'Unknown variant ''%s''. Use RZP | RCP.', variant);
    end
end

end
