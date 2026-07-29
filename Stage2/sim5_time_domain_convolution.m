function results = sim5_time_domain_convolution(p, doPlot)
% =========================================================================
% Program : sim5_time_domain_convolution.m
% Description :
%   Sim-5  "時域通道兩種實作一致 (4.4.1)"  (chap4_list.md Sec. D, Stage 2)
%
%   Purpose  : verify that the MATRIX form r = G*s and the direct
%              time-varying CONVOLUTION form r[q] = sum_l g^s[l,q]s[q-l]
%              are two independent derivations of the same equation (4.36).
%   Theory   : r = G*s  ==  r[q] = sum_{l=0}^{l_max} g^s[l,q]*s[q-l]  (RZP:
%              q-l<0 treated as zero, i.e. the isolated single-frame
%              assumption of Sec. 4.4 — see chap4_list.md Sec. H point 7).
%   Method   : build (g_i,l_i,k_i) -> g^s[l,q] (gen_gs.m) -> G (gen_G.m,
%              RZP) and independently apply_channel_conv.m (pure for-loop).
%   Criterion: norm(G*s - r_conv)/norm(r_conv) < 1e-12.
%
% Input  : p      - parameter struct from otfs_params()  (optional)
%          doPlot - true to draw imagesc(abs(G)), Fig. 4.8 style (default false)
% Output : results - struct with .name .err_rel .tol .pass, plus .G .gs .taps
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Stage1'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;  N = p.N;

fprintf('\n=== Sim-5 : time-domain channel, matrix G vs. direct convolution (M=%d, N=%d) ===\n', M, N);

[g, l, k, taps_info] = gen_channel_taps(p);
gs = gen_gs(p, g, l, k);
G  = gen_G(p, gs, 'RZP');

X  = gen_dd_symbols(p);
Pm = gen_perm_matrix(M, N);
s  = otfs_modulate(X, 'idzt', Pm);

r_mat  = G * s;
r_conv = apply_channel_conv(p, gs, s, 'RZP');

err = relerr(r_mat, r_conv);

fprintf('    paths: l_i = [%s],  k_i = [%s],  l_max = %d, k_max = %d\n', ...
        num2str(l), num2str(k), taps_info.l_max, taps_info.k_max);
fprintf('    ||G*s - r_conv|| / ||r_conv|| = %.3e\n', err);

results.name    = 'Sim-5 time-domain channel: G*s vs. direct convolution';
results.err_rel = err;
results.tol     = 1e-12;
results.pass    = err < results.tol;
results.G       = G;
results.gs      = gs;
results.taps    = struct('g', g, 'l', l, 'k', k);

fprintf('    -> %s\n', ternary(results.pass, 'PASS', 'FAIL'));

if doPlot
    figure('Name', 'Sim-5 |G| (RZP)', 'Color', 'w');
    imagesc(abs(G)); axis square; colorbar;
    xlabel('column index q'''); ylabel('row index q');
    title(sprintf('|G|  (RZP, Fig. 4.8 style),  l_{max}=%d -> %d sub-diagonals', ...
                  taps_info.l_max, taps_info.l_max+1));
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
