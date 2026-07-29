function results = sim10_symbol_level_theory(p, doPlot)
% =========================================================================
% Program : sim10_symbol_level_theory.m
% Description :
%   Sim-10  "符號級 DD 域理論式驗證（最終驗收）"  (chap4_list.md Sec. D, Stage 2)
%
%   Purpose  : THE formal validation of Stage 2 — compare the pure symbol-
%              level "theory generator" (4.118)-(4.122) directly against
%              the full matrix-pipeline simulation (channel matrix G -> DZT
%              demodulation), with NO shared code path between the two.
%   Theory   : rectangular pulse, integer taps, RCP channel (chap4_list.md
%              Sec. E.1, code 16 -> Sim-10):
%                Y[m,n] = sum_i g_i z^{k_i(m-l_i)} X[[m-l_i]_M,[n-k_i]_N]      , m>=l_i
%                Y[m,n] = sum_i g_i z^{k_i([m-l_i]_M)} e^{-j2pi n/N} X[...]    , m< l_i
%              implemented in theory_dd_relation.m.
%   Method   : X -> s -> r = G_RCP*s -> Y_sim (via otfs_demodulate 'dzt');
%              independently Y_theory = theory_dd_relation(p,X,g,l,k);
%              compare.
%   Criterion: norm(Y_sim - Y_theory)/norm(Y_theory) < 1e-10.
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw |Y_sim|,|Y_theory|,|diff| heatmaps
% Output : results - struct with .name .err_rel .tol .pass
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Stage1'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;  N = p.N;

fprintf('\n=== Sim-10 : symbol-level DD-domain theory vs. simulation (RCP)  (M=%d, N=%d) ===\n', M, N);

[g, l, k] = gen_channel_taps(p);
gs = gen_gs(p, g, l, k);
G  = gen_G(p, gs, 'RCP');          % RCP variant required to match (4.118)-(4.122)

X  = gen_dd_symbols(p);
Pm = gen_perm_matrix(M, N);
s  = otfs_modulate(X, 'idzt', Pm);
r  = G * s;
Y_sim = otfs_demodulate(r, M, N, 'dzt', Pm);

Y_theory = theory_dd_relation(p, X, g, l, k);

err = relerr(Y_sim, Y_theory);
fprintf('    ||Y_sim - Y_theory|| / ||Y_theory|| = %.3e\n', err);

results.name    = 'Sim-10 symbol-level DD-domain theory vs. simulation (final Stage-2 validation)';
results.err_rel = err;
results.tol     = 1e-10;
results.pass    = err < results.tol;
results.Y_sim    = Y_sim;
results.Y_theory = Y_theory;

fprintf('    -> %s\n', ternary(results.pass, 'PASS', 'FAIL'));

if doPlot
    figure('Name', 'Sim-10 Y_sim vs Y_theory', 'Color', 'w');
    subplot(1,3,1); imagesc(0:N-1, 0:M-1, abs(Y_sim));    axis square; colorbar;
    xlabel('n'); ylabel('m'); title('|Y_{sim}|');
    subplot(1,3,2); imagesc(0:N-1, 0:M-1, abs(Y_theory)); axis square; colorbar;
    xlabel('n'); ylabel('m'); title('|Y_{theory}|  (4.120)/(4.121)');
    subplot(1,3,3); imagesc(0:N-1, 0:M-1, abs(Y_sim-Y_theory)); axis square; colorbar;
    xlabel('n'); ylabel('m'); title('|Y_{sim} - Y_{theory}|  (~1e-15)');
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
