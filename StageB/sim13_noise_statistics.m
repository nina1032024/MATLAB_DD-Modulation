function results = sim13_noise_statistics(p, numTrials, doPlot)
% =========================================================================
% Program : sim13_noise_statistics.m
% Description :
%   Sim-13  "雜訊統計一致性"  (chap4_list.md Sec. D, Stage 2)
%
%   Purpose  : confirm that AWGN keeps the same statistics, CN(0, sigma_w^2 I),
%              in every domain — because the operators that move a signal
%              between domains (P, I_N⊗F_M, I_M⊗F_N, ...) are all UNITARY,
%              which preserves both the per-sample variance and the
%              (identity) spatial correlation structure of white noise.
%   Theory   : w (time), w-check = (I_N⊗F_M)w (time-freq), w-tilde = P^T w
%              (delay-time), z = (I_M⊗F_N)w-tilde (delay-Doppler) are all
%              CN(0, sigma_w^2 I_NM).
%   Method   : generate numTrials iid realizations of w, push each through
%              the three unitary transforms, and for every domain estimate
%              (i) the per-sample variance (diagonal of the sample
%              correlation matrix) and (ii) the size of the off-diagonal
%              (cross-sample) correlation relative to sigma_w^2.
%   Criterion: max per-sample variance relative error < 1% ; RMS
%              off-diagonal correlation / sigma_w^2 < 5% (finite-sample
%              statistical threshold, not machine precision — it shrinks
%              like 1/sqrt(numTrials)).
%
% Input  : p         - parameter struct from otfs_params()   (optional)
%          numTrials - number of noise realizations           (default 1e5)
%          doPlot    - true to draw a per-domain variance bar chart
% Output : results - struct with .name .err_rel .tol .pass, plus
%              .varErr (1x4), .corrErr (1x4), .names
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Stage1'));

if nargin < 1 || isempty(p),         p = otfs_params(); end
if nargin < 2 || isempty(numTrials), numTrials = 1e5;   end
if nargin < 3 || isempty(doPlot),    doPlot = true;    end

M = p.M;  N = p.N;  NM = p.NM;
Pm = gen_perm_matrix(M, N);
IN_FM = kron(eye(N), p.F_M);
IM_FN = kron(eye(M), p.F_N);

sigma2 = 1;                          % normalized noise variance
rng(p.seed + 7);
W = sqrt(sigma2/2) * (randn(NM, numTrials) + 1j*randn(NM, numTrials));

fprintf('\n=== Sim-13 : noise statistics across domains (M=%d, N=%d, trials=%d) ===\n', M, N, numTrials);

Wcheck = IN_FM * W;
Wtilde = Pm.' * W;
Z      = IM_FN * Wtilde;

domains = {W, Wcheck, Wtilde, Z};
names   = {'w (time)', 'w_check (time-freq)', 'w_tilde (delay-time)', 'z (delay-Doppler)'};

varErr  = zeros(1, 4);
corrErr = zeros(1, 4);
for i = 1:4
    Wd = domains{i};
    Rcorr = (Wd * Wd') / numTrials;
    varErr(i)  = max(abs(diag(Rcorr) - sigma2)) / sigma2;
    offRMS     = norm(Rcorr - diag(diag(Rcorr)), 'fro') / sqrt(NM*(NM-1));
    corrErr(i) = offRMS / sigma2;
    fprintf('    %-24s : var. rel. err = %.4f%%   off-diag corr / sigma2 = %.4f\n', ...
            names{i}, 100*varErr(i), corrErr(i));
end

results.name    = 'Sim-13 noise statistics across domains';
results.varErr  = varErr;
results.corrErr = corrErr;
results.names   = names;
results.err_rel = max(varErr);
results.tol     = 0.01;
results.pass    = all(varErr < 0.01) && all(corrErr < 0.05);

fprintf('    -> %s\n', ternary(results.pass, 'PASS', 'FAIL'));

if doPlot
    figure('Name', 'Sim-13 noise variance per domain', 'Color', 'w');
    bar(100*varErr); grid on;
    set(gca, 'XTick', 1:4, 'XTickLabel', names, 'XTickLabelRotation', 20);
    ylabel('variance relative error [%]');
    yline(1, 'r--', '1% threshold');
    title('Sim-13 : per-domain noise variance stays at \sigma_w^2 under unitary transforms');
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
