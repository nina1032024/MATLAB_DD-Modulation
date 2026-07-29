function results = sim11_ideal_vs_rectangular_pulse(p, doPlot)
% =========================================================================
% Program : sim11_ideal_vs_rectangular_pulse.m
% Description :
%   Sim-11  "理想脈波 vs. 矩形脈波比較 (4.2 節)"  (chap4_list.md Sec. D, Stage 2)
%
%   Purpose  : contrast the ideal (bi-orthogonal) pulse relation (4.15) —
%              a pure, time-invariant 2D circular convolution — against the
%              rectangular pulse's relation (4.120)/(4.121), which carries
%              an extra phase z^{k_i(m-l_i)} and, for m<l_i, an extra
%              e^{-j2*pi*n/N} boundary term. This is the numerical evidence
%              behind chap4_list.md item C.8 ("雙正交性喪失").
%
%   NOTE ON THE PASS CRITERION: chap4_list.md gives Sim-11 only qualitative
%   "觀察重點" bullets, no numeric threshold. This script therefore reports
%   three DIFFERENT comparisons rather than forcing one pass/fail number:
%     (1) general case (default 3-path channel)       -> reported, informational
%     (2) all k_i = 0, l_max > 0 (delay only, no Doppler)
%                                                      -> reported, informational
%     (3) l_i = 0 AND k_i = 0 (trivial identity channel)
%                                                      -> must match EXACTLY
%
%   *** Finding (do not "fix" this away) ***: a naive reading of
%   "當所有 k_i=0 時兩者應完全相同" is NOT correct in general. The
%   rectangular-pulse formula's phase term z^{k_i(m-l_i)} vanishes when
%   k_i=0 (good), but its m<l_i branch STILL carries the extra
%   e^{-j2*pi*n/N} boundary phase regardless of k_i — that phase comes
%   from the frame-boundary time-shift itself (see theory_dd_relation.m),
%   not from Doppler. So case (2) generally shows a NONZERO residual
%   whenever l_max>0. The two formulas coincide exactly only when there is
%   truly no channel distortion at all: l_i=0 (no delay, so the m<l_i
%   branch never triggers) AND k_i=0 (no phase) — case (3). This was
%   confirmed by first-principles re-derivation of (4.120)/(4.121) from
%   the RCP channel model (see theory_dd_relation.m header / progress
%   report), and matches Sim-10's independently-validated formula.
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw the |Y_rect - Y_ideal| heatmap (default false)
% Output : results - struct with .name .tol .pass (gated on case 3 only),
%              .diff_general, .diff_zero_doppler, .diff_trivial
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Stage1'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;  N = p.N;

fprintf('\n=== Sim-11 : ideal pulse (4.15) vs. rectangular pulse (4.120/4.121) (M=%d, N=%d) ===\n', M, N);

X = gen_dd_symbols(p);

% ---- (1) general case: default 3-path channel ---------------------------
[g, l, k] = gen_channel_taps(p);
Y_rect  = theory_dd_relation(p, X, g, l, k);
Y_ideal = ideal_pulse_relation(p, X, g, l, k);
diff_general = norm(Y_rect - Y_ideal, 'fro') / norm(Y_ideal, 'fro');

% ---- (2) all k_i = 0, but l_max > 0 (delay wraparound still present) ----
[g0, l0, k0] = gen_channel_taps(p, 'l', l, 'k', zeros(size(k)), 'randomPhase', false);
Y_rect0  = theory_dd_relation(p, X, g0, l0, k0);
Y_ideal0 = ideal_pulse_relation(p, X, g0, l0, k0);
diff_zero_doppler = norm(Y_rect0 - Y_ideal0, 'fro') / norm(Y_ideal0, 'fro');

% ---- (3) l_i = 0 AND k_i = 0 : trivial identity channel -----------------
[g1, l1, k1] = gen_channel_taps(p, 'l', 0, 'k', 0, 'randomPhase', false);
Y_rect1  = theory_dd_relation(p, X, g1, l1, k1);
Y_ideal1 = ideal_pulse_relation(p, X, g1, l1, k1);
diff_trivial = norm(Y_rect1 - Y_ideal1, 'fro') / norm(Y_ideal1, 'fro');

fprintf('    (1) general 3-path channel              : rel. diff = %.3e   (informational)\n', diff_general);
fprintf('    (2) all k_i=0, l_max=%d (still wraps)     : rel. diff = %.3e   (informational — NOT expected to vanish)\n', ...
        max(l), diff_zero_doppler);
fprintf('    (3) l_i=0 AND k_i=0 (trivial channel)   : rel. diff = %.3e   (must be ~0)\n', diff_trivial);

results.name              = 'Sim-11 ideal vs. rectangular pulse (4.2)';
results.diff_general      = diff_general;
results.diff_zero_doppler = diff_zero_doppler;
results.diff_trivial      = diff_trivial;
results.err_rel           = diff_trivial;     % the one gate-able number
results.tol               = 1e-10;
results.pass              = diff_trivial < results.tol;

fprintf('    -> %s  (gated on case (3) only; (1)/(2) are observations, per chap4_list.md having no numeric threshold here)\n', ...
        ternary(results.pass, 'PASS', 'FAIL'));

if doPlot
    figure('Name', 'Sim-11 ideal vs rectangular pulse', 'Color', 'w');
    subplot(1,3,1); imagesc(0:N-1, 0:M-1, abs(Y_rect));  axis square; colorbar;
    xlabel('n'); ylabel('m'); title('|Y_{rect}|  (4.120/4.121)');
    subplot(1,3,2); imagesc(0:N-1, 0:M-1, abs(Y_ideal)); axis square; colorbar;
    xlabel('n'); ylabel('m'); title('|Y_{ideal}|  (4.15, circular conv.)');
    subplot(1,3,3); imagesc(0:N-1, 0:M-1, abs(Y_rect-Y_ideal)); axis square; colorbar;
    xlabel('n'); ylabel('m'); title('|Y_{rect} - Y_{ideal}|');
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
