function results = sim_appendixC_channel_models(p, doPlot)
% =========================================================================
% Program : sim_appendixC_channel_models.m
% Description :
%   Supplementary Sim (not numbered Sim-1~16): validate the 3GPP standard
%   channel models (Appendix C, MATLAB code 5 + code 6 — EPA/EVA/ETU),
%   which no Sim in Stage 1-3 exercised before this file (every earlier
%   Sim used only synthetic on-grid/off-grid taps from gen_channel_taps.m).
%
%   Four checks:
%     (a) PDP normalization: sum(pdp_linear) = 1 exactly (deterministic,
%         before random Rayleigh fading is applied).
%     (b) Monte Carlo validation of the Rayleigh gain model itself:
%         E|g_i|^2 should converge to pdp_linear(i) (g_i = sqrt(pdp_i)*CN(0,1)).
%     (c) THE delay-resolution collapse this predicts (chap4_list.md Sec.
%         H point 5 / Sec. E.3.5): "M=64, Δf=15kHz 時，EVA 的 9 條路徑經
%         round(delays/delay_resolution) 後會塌成 l_i=[0 0 0 0 0 1 1 2 2]，
%         只剩 3 個相異 delay tap" — reproduced numerically here (not just
%         asserted) across a sweep of M, confirming the collapse and its
%         recovery at larger M.
%     (d) integration test: feed a realistic EVA channel (at M=64, where
%         the collapse above is severe) through the ALREADY-VALIDATED
%         Stage-2 pipeline (gen_gs/gen_G/build_channel_matrices) and rerun
%         the Sim-5 (self-consistency) and Sim-9 (circulant, "safe" m>=l
%         blocks) style checks — confirming that machinery generalizes
%         beyond the toy synthetic taps used in Sim-5~13.
%
% Input  : p      - parameter struct from otfs_params()   (optional; only
%                    df, fc, seed are taken from it — M/N are swept/reset
%                    internally for parts (c)/(d))
%          doPlot - true to draw the distinct-tap-count vs M curve
% Output : results - struct with .name .err_rel .tol .pass, plus
%              .normErr, .mcErr, .Mlist, .distinctCount,
%              .err_selfconsistency, .circErrSafe
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..', 'Stage1'));
addpath(fullfile(thisDir, '..'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;      end

fprintf('\n=== Appendix C supplementary Sim : 3GPP standard channels (code 5+6) ===\n');

models = {'EPA', 'EVA', 'ETU'};

% ---- (a) PDP normalization (deterministic) ------------------------------
normErr = zeros(1, 3);
for i = 1:3
    [~, ~, ~, info] = gen_3gpp_channel(p, models{i});
    normErr(i) = abs(sum(info.pdp_linear) - 1);
    fprintf('  (a) %-4s: sum(pdp_linear) = %.15f  (taps=%d)\n', models{i}, sum(info.pdp_linear), info.taps);
end

% ---- (b) Monte Carlo Rayleigh gain validation (EVA) ---------------------
nMC = 5000;
[~, ~, ~, infoEVA] = gen_3gpp_channel(p, 'EVA');
taps = infoEVA.taps;
gPow = zeros(nMC, taps);
for tr = 1:nMC
    [g, ~, ~, ~] = gen_3gpp_channel(p, 'EVA', [], tr);
    gPow(tr, :) = abs(g).^2;
end
empMean = mean(gPow, 1);
mcErr = max(abs(empMean - infoEVA.pdp_linear) ./ infoEVA.pdp_linear);
fprintf('  (b) Monte Carlo (EVA, %d trials): max rel. err, E|g_i|^2 vs pdp_linear(i) = %.3e\n', nMC, mcErr);

% ---- (c) distinct delay taps vs M (delay-resolution collapse) ----------
Mlist = [8 32 64 128 256 512];
distinctCount = zeros(numel(models), numel(Mlist));
for mi = 1:numel(models)
    for j = 1:numel(Mlist)
        pj = otfs_params('M', Mlist(j), 'N', p.N, 'df', p.df, 'fc', p.fc, 'seed', p.seed);
        [~, ~, ~, infoj] = gen_3gpp_channel(pj, models{mi});
        distinctCount(mi, j) = infoj.distinct_l;
    end
    fprintf('  (c) %-4s: M=%-24s -> distinct l_i = %-24s (of %d taps)\n', ...
            models{mi}, mat2str(Mlist), mat2str(distinctCount(mi, :)), infoj.taps);
end

p64 = otfs_params('M', 64, 'N', p.N, 'df', p.df, 'fc', p.fc, 'seed', p.seed);
[~, l64, ~, ~] = gen_3gpp_channel(p64, 'EVA');
expected_l64 = [0 0 0 0 0 1 1 2 2];
matches_reference = isequal(l64, expected_l64);
fprintf('  (c) EVA @ M=64: l_i = %s  (chap4_list.md Sec.E.3.5 states [0 0 0 0 0 1 1 2 2]) -> %s\n', ...
        mat2str(l64), ternary(matches_reference, 'MATCHES', 'differs'));

% ---- (d) integration test: realistic channel through Stage-2 pipeline --
p_int = p64;
[g, l, k] = gen_3gpp_channel(p_int, 'EVA');
gs = gen_gs(p_int, g, l, k);
G  = gen_G(p_int, gs, 'RZP');

X  = gen_dd_symbols(p_int);
Pm = gen_perm_matrix(p_int.M, p_int.N);
s  = otfs_modulate(X, 'idzt', Pm);
r_mat  = G * s;
r_conv = apply_channel_conv(p_int, gs, s, 'RZP');
err_selfconsistency = relerr(r_mat, r_conv);

Cm = build_channel_matrices(p_int, G, Pm);
x = reshape(X.', p_int.NM, 1);
l_max = max(l);
circErrSafe = 0;
for m = 0:p_int.M-1
    for ll = 0:l_max
        if m >= ll
            mp = mod(m - ll, p_int.M);
            K = get_block(Cm.H, p_int.N, m, mp);
            c = K(:, 1);
            Kcirc = zeros(p_int.N, p_int.N);
            for col = 0:p_int.N-1
                Kcirc(:, col+1) = circshift(c, col);
            end
            circErrSafe = max(circErrSafe, norm(K - Kcirc, 'fro') / max(norm(K, 'fro'), eps));
        end
    end
end

fprintf('  (d) integration test (EVA @ M=%d, N=%d, real channel through Stage-2 pipeline)\n', p_int.M, p_int.N);
fprintf('      self-consistency ||G*s - r_conv|| / ||r_conv||  = %.3e\n', err_selfconsistency);
fprintf('      circulant check ("safe" m>=l blocks)             = %.3e\n', circErrSafe);

results.name             = 'Appendix C 3GPP channel models (code 5+6)';
results.normErr          = normErr;
results.mcErr            = mcErr;
results.Mlist            = Mlist;
results.distinctCount    = distinctCount;
results.err_selfconsistency = err_selfconsistency;
results.circErrSafe      = circErrSafe;
results.err_rel          = max([normErr, err_selfconsistency, circErrSafe]);
results.tol              = 1e-9;
results.pass             = (max(normErr) < 1e-9) && (err_selfconsistency < 1e-9) && (circErrSafe < 1e-9);

fprintf('  -> %s  (gated on (a)/(d); (b) Monte Carlo and (c) distinct-tap counts are informational)\n', ...
        ternary(results.pass, 'PASS', 'FAIL'));

if doPlot
    figure('Name', 'Appendix C 3GPP channel models', 'Color', 'w');
    semilogx(Mlist, distinctCount(1,:), '-o', 'LineWidth', 1.5); hold on;
    semilogx(Mlist, distinctCount(2,:), '-s', 'LineWidth', 1.5);
    semilogx(Mlist, distinctCount(3,:), '-^', 'LineWidth', 1.5);
    grid on; xlabel('M (delay bins)'); ylabel('distinct delay taps after quantization');
    legend(models, 'Location', 'southeast');
    title('Distinct l_i vs. M — delay-resolution collapse (Sec. H.5 / E.3.5)');
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
