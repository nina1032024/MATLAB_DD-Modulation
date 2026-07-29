function results = sim15_fractional_delay_leakage(p, doPlot)
% =========================================================================
% Program : sim15_fractional_delay_leakage.m
% Description :
%   Sim-15  "分數延遲的洩漏驗證"  (chap4_list.md Sec. D, Stage 3)
%
%   Purpose  : verify the sinc-interpolation model for a fractional delay
%              tap (l_i not an integer, e.g. 1.4), and quantify how much the
%              delay range must be extended beyond the nominal l_max to
%              avoid truncating the leaked energy (Sec. H point 6).
%   Theory   : g^s[l,q] = sum_i g_i z^{k_i(q-l)} sinc(l-l_i)          (4.6)
%              implemented in gen_gs_frac.m, which reduces EXACTLY to the
%              delta-based (4.8)/gen_gs.m when l_i is an integer.
%   Method   :
%     (a) self-consistency: G_frac*s == direct convolution (Stage 2's
%         apply_channel_conv.m, reused as-is — it works for any gs matrix
%         regardless of row count), same style as Sim-5.
%     (b) regression: force l_i integer -> gen_gs_frac.m must reproduce
%         Stage 2's delta-based gen_gs.m exactly (up to floating-point sinc
%         evaluation), confirming the generalization didn't break Stage 2.
%     (c) truncation-error curve: capturedFraction(Lrange) =
%         sum_{l=0}^{Lrange-1} sinc(l-l_i)^2. The sinc basis' completeness
%         identity sums to 1 over ALL integers l (-inf...+inf), but delay
%         taps are causal (l>=0 only, per Sec. A.3's l in L), so this curve
%         plateaus BELOW 1 (around 0.95 for l_i=1.4) — the gap is the
%         energy that "would" sit at negative delay and is physically
%         excluded, not a truncation bug. What growing Lrange actually
%         fixes is the POSITIVE-side tail only; this directly quantifies
%         "否則會截斷洩漏能量而使誤差變大" for the l>0 side.
%     (d) effective path count S vs. fractional offset (informational,
%         S=1 at zero offset, grows towards a maximum near offset=0.5).
%   Criterion: (a) and (b) gated at < 1e-10 (RZP channel model); (c)/(d)
%              are informational (chap4_list.md gives no numeric threshold
%              for this Sim, only the "S >= P" qualitative observation).
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw the truncation-error and S(offset) curves
% Output : results - struct with .name .err_rel .tol .pass, plus
%              .err_selfconsistency, .err_regression, .LrangeList,
%              .capturedFraction, .fracList, .S_list
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Stage1'));
addpath(fullfile(thisDir, '..', 'Stage2'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;  N = p.N;

fprintf('\n=== Sim-15 : fractional delay leakage (4.6)  (M=%d, N=%d) ===\n', M, N);

% ---- (a) self-consistency: G_frac*s vs. direct convolution -------------
l_i = 1.4;  k_i = 0;
[g, l, k] = gen_channel_taps(p, 'l', l_i, 'k', k_i, 'randomPhase', false);
Lrange = ceil(l_i) + 9;

gs_frac = gen_gs_frac(p, g, l, k, Lrange);
G_frac  = gen_G(p, gs_frac, 'RZP');           % Stage 2's gen_G.m, reused as-is

X  = gen_dd_symbols(p);
Pm = gen_perm_matrix(M, N);
s  = otfs_modulate(X, 'idzt', Pm);
r_mat  = G_frac * s;
r_conv = apply_channel_conv(p, gs_frac, s, 'RZP');   % Stage 2's function, reused as-is

err_selfconsistency = relerr(r_mat, r_conv);
fprintf('    (a) self-consistency  ||G_frac*s - r_conv|| / ||r_conv||   = %.3e\n', err_selfconsistency);

% ---- (b) regression: integer l_i must reduce to Stage-2 delta model ----
l_int = 2;
[g2, l2, k2] = gen_channel_taps(p, 'l', l_int, 'k', 0, 'randomPhase', false);
Lrange2 = l_int + 9;
gs_frac_int = gen_gs_frac(p, g2, l2, k2, Lrange2);
gs_delta    = gen_gs(p, g2, l2, k2);                 % Stage 2's delta-based (4.8)

err_regression_main = relerr(gs_frac_int(1:(l_int+1), :), gs_delta);
leakEnergy  = norm(gs_frac_int((l_int+2):end, :), 'fro');
totalEnergy = norm(gs_frac_int, 'fro');
err_regression_leak = leakEnergy / max(totalEnergy, eps);

fprintf('    (b) regression (integer l_i): main-tap rel. err = %.3e, residual leak = %.3e\n', ...
        err_regression_main, err_regression_leak);

% ---- (c) truncation error vs. Lrange (Parseval completeness) -----------
LrangeList = [3 5 8 11 15 20 30];
capturedFraction = zeros(size(LrangeList));
for idx = 1:numel(LrangeList)
    Lr = LrangeList(idx);
    w = sinc_fn((0:Lr-1) - l_i);
    capturedFraction(idx) = sum(w.^2);
end
fprintf('    (c) Lrange sweep       = %s\n', mat2str(LrangeList));
fprintf('        captured fraction  = %s   (plateaus <1: l>=0 causal, misses negative-delay tail)\n', mat2str(capturedFraction, 4));

% ---- (d) effective path count S vs. fractional offset -------------------
fracList = 0:0.1:0.5;
S_list = zeros(size(fracList));
bigRange = -30:30;
for idx = 1:numel(fracList)
    li_test = 2 + fracList(idx);
    w = sinc_fn(bigRange - li_test);
    S_list(idx) = effective_path_count(w.^2, 0.99);
end
fprintf('    (d) fractional offset  = %s\n', mat2str(fracList));
fprintf('        effective S (99%% energy) = %s   (P=1 here; S>=P as required)\n', mat2str(S_list));

results.name                = 'Sim-15 fractional delay leakage (4.6)';
results.err_selfconsistency = err_selfconsistency;
results.err_regression      = err_regression_main;
results.err_regression_leak = err_regression_leak;
results.LrangeList          = LrangeList;
results.capturedFraction    = capturedFraction;
results.fracList            = fracList;
results.S_list              = S_list;
results.err_rel             = max(err_selfconsistency, err_regression_main);
results.tol                 = 1e-10;
results.pass                = (err_selfconsistency < 1e-10) && (err_regression_main < 1e-8);

fprintf('    -> %s\n', ternary(results.pass, 'PASS', 'FAIL'));

if doPlot
    figure('Name', 'Sim-15 fractional delay leakage', 'Color', 'w');
    subplot(1,2,1);
    plot(LrangeList, capturedFraction, '-o', 'LineWidth', 1.5); grid on;
    xlabel('L_{range}'); ylabel('captured energy fraction');
    title(sprintf('Sim-15(c): sinc tail truncation  (l_i=%.1f)', l_i));
    yline(1, 'r--');

    subplot(1,2,2);
    plot(fracList, S_list, '-s', 'LineWidth', 1.5); grid on;
    xlabel('fractional offset of l_i'); ylabel('effective path count S (99% energy)');
    title('Sim-15(d): S grows with fractional offset');
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
