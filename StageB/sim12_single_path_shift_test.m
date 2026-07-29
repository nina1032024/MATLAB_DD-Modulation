function results = sim12_single_path_shift_test(p, doPlot)
% =========================================================================
% Program : sim12_single_path_shift_test.m
% Description :
%   Sim-12  "單路徑位移測試（強烈建議先做）"  (chap4_list.md Sec. D, Stage 2)
%
%   Purpose  : THE most direct debugging tool for axis/sign errors. With a
%              single path (P=1, g_1) and a delta symbol X = delta at
%              (m0,n0), the DD-domain output must be a single delta at
%              ([m0+l1]_M, [n0+k1]_N) with magnitude |g1| and phase
%              z^{k1(m-l1)} — any swapped axis or flipped sign shows up
%              immediately as the peak landing in the wrong place.
%   Method   : three cases — pure delay (k1=0), pure Doppler (l1=0), and
%              both — each checked TWO independent ways:
%                (a) theory_dd_relation.m (symbol-level formula, RCP)
%                (b) full pipeline: X -> s -> r=G_RCP*s -> Y (otfs_demodulate)
%              For each: verify peak position, peak magnitude = |g1|,
%              peak phase (for m0+l1 < M, i.e. no wraparound, so the plain
%              m>=l1 branch applies: phase = z^{k1(m-l1)}), and that all
%              other DD bins carry (numerically) zero energy.
%   Criterion: position exact; magnitude & phase relative error < 1e-10;
%              energy elsewhere / total energy < 1e-10; theory vs. full
%              pipeline agree to < 1e-10.
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw |X| / |Y| for each case   (default false)
% Output : results - struct with .name .err_rel .tol .pass, plus .cases
%              (1x3 struct array: delay-only, Doppler-only, both)
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Stage1'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;  N = p.N;
Pm = gen_perm_matrix(M, N);

m0 = 2;  n0 = 1;                      % interior point (avoids m0=0 edge cases)
X = zeros(M, N);  X(m0+1, n0+1) = 1;

fprintf('\n=== Sim-12 : single-path delta shift test  (M=%d, N=%d, delta at m0=%d,n0=%d) ===\n', M, N, m0, n0);

caseDefs = struct('name', {'pure delay', 'pure Doppler', 'delay + Doppler'}, ...
                   'l1',   {3, 0, 2}, ...
                   'k1',   {0, 2, -1});

cases = struct('name', {}, 'err', {}, 'expected_pos', {}, 'peak_pos', {});
maxErr = 0;

for c = 1:numel(caseDefs)
    l1 = caseDefs(c).l1;  k1 = caseDefs(c).k1;
    [g1, l1v, k1v] = gen_channel_taps(p, 'l', l1, 'k', k1, 'randomPhase', false);

    mExp = mod(m0 + l1, M);
    nExp = mod(n0 + k1, N);

    % (a) theory generator
    Y_theory = theory_dd_relation(p, X, g1, l1v, k1v);

    % (b) full channel pipeline (RCP, matches theory_dd_relation's model)
    gs = gen_gs(p, g1, l1v, k1v);
    G  = gen_G(p, gs, 'RCP');
    s  = otfs_modulate(X, 'idzt', Pm);
    r  = G * s;
    Y_pipe = otfs_demodulate(r, M, N, 'dzt', Pm);

    err_theory_vs_pipe = relerr(Y_theory, Y_pipe);

    [peakMag, peakIdx] = max(abs(Y_theory(:)));
    [mPeak, nPeak] = ind2sub([M, N], peakIdx);
    mPeak = mPeak - 1;  nPeak = nPeak - 1;
    posOK = (mPeak == mExp) && (nPeak == nExp);

    magErr = abs(peakMag - abs(g1)) / abs(g1);

    % expected phase: since m0+l1 < M here (no wraparound), m=mExp >= l1
    expPhase = p.z ^ (k1v * (mExp - l1v));
    phaseErr = abs(Y_theory(mExp+1, nExp+1) - g1*expPhase) / abs(g1);

    energyTotal    = norm(Y_theory, 'fro')^2;
    energyElsewhere = energyTotal - abs(Y_theory(mExp+1, nExp+1))^2;
    energyRatio = energyElsewhere / energyTotal;

    err_c = max([err_theory_vs_pipe, double(~posOK), magErr, phaseErr, energyRatio]);
    maxErr = max(maxErr, err_c);

    fprintf('  [%s]  l1=%d, k1=%d\n', caseDefs(c).name, l1, k1);
    fprintf('      expected peak at (m,n) = (%d,%d), got (%d,%d)  -> %s\n', ...
            mExp, nExp, mPeak, nPeak, ternary(posOK, 'OK', 'MISMATCH'));
    fprintf('      |peak|=%.6f vs |g1|=%.6f  (rel err %.3e);  phase rel err %.3e\n', ...
            peakMag, abs(g1), magErr, phaseErr);
    fprintf('      energy elsewhere / total = %.3e;  theory vs. full-pipeline = %.3e\n', ...
            energyRatio, err_theory_vs_pipe);

    cases(c).name         = caseDefs(c).name;
    cases(c).l1           = l1;
    cases(c).k1           = k1;
    cases(c).expected_pos = [mExp, nExp];
    cases(c).peak_pos     = [mPeak, nPeak];
    cases(c).err          = err_c;
    cases(c).Y            = Y_theory;

    if doPlot
        figure('Name', sprintf('Sim-12 %s', caseDefs(c).name), 'Color', 'w');
        subplot(1,2,1); imagesc(0:N-1, 0:M-1, abs(X));       axis square; colorbar;
        xlabel('n'); ylabel('m'); title('|X|  (delta input)');
        subplot(1,2,2); imagesc(0:N-1, 0:M-1, abs(Y_theory)); axis square; colorbar;
        xlabel('n'); ylabel('m');
        title(sprintf('|Y|  %s (l1=%d,k1=%d)', caseDefs(c).name, l1, k1));
    end
end

results.name    = 'Sim-12 single-path delta shift test';
results.cases   = cases;
results.err_rel = maxErr;
results.tol     = 1e-10;
results.pass    = maxErr < results.tol;

fprintf('\n    max error over all cases = %.3e   -> %s\n', maxErr, ternary(results.pass, 'PASS', 'FAIL'));

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
