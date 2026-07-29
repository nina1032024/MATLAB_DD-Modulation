function results = sim1_perfect_reconstruction(p, doPlot)
% =========================================================================
% Program : sim1_perfect_reconstruction.m
% Description :
%   Sim-1  "IDZT/DZT perfect reconstruction"  (chap4_list.md Sec. D, Stage 1)
%
%   Purpose  : verify that the modulator and demodulator are exact inverses,
%              i.e. that the sqrt(M)/sqrt(N) normalizations are consistent
%              (Sec. H.1 - the #1 source of errors).
%   Theory   : with no channel and no noise (r = s),  Y = X exactly.
%   Method   : X -> s -> r = s -> Y, then measure norm(Y - X, 'fro').
%              ALL 4 modulator paths x ALL 4 demodulator paths are tested
%              (16 combinations), so a normalization slip in any single
%              implementation shows up immediately.
%   Criterion: norm(Y - X, 'fro') < 1e-12   (and relative error < 1e-12)
%
% Input  : p      - parameter struct from otfs_params()  (optional)
%          doPlot - true to draw the |X| vs |Y| comparison   (default false)
% Output : results - struct with fields
%              .name       'Sim-1'
%              .err_abs    max absolute Frobenius error over all 16 combos
%              .err_rel    max relative error over all 16 combos
%              .tol        threshold used
%              .pass       logical
%              .table      4x4 matrix of absolute errors (mod x demod)
%              .modList / .demodList   method name lists
% =========================================================================

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;  N = p.N;

% ---- transmit symbols -------------------------------------------------
[X, ~, info] = gen_dd_symbols(p);
P = gen_perm_matrix(M, N);

modList   = {'isfft', 'idzt', 'elem', 'vec'};
demodList = {'wigner', 'dzt',  'elem', 'vec'};

errAbs = zeros(numel(modList), numel(demodList));
errRel = zeros(numel(modList), numel(demodList));

fprintf('\n=== Sim-1 : IDZT/DZT perfect reconstruction (M=%d, N=%d, %d-QAM/PSK) ===\n', ...
        M, N, p.modOrder);
fprintf('    no channel, no noise  ->  r = s  ->  Y must equal X\n\n');

for i = 1:numel(modList)
    for j = 1:numel(demodList)

        s = otfs_modulate(X, modList{i}, P);        % X -> s   (NM x 1)
        r = s;                                      % ideal channel: r = s
        Y = otfs_demodulate(r, M, N, demodList{j}, P);   % r -> Y  (M x N)

        errAbs(i,j) = norm(Y - X, 'fro');
        errRel(i,j) = relerr(Y, X);
    end
end

% ---- report -----------------------------------------------------------
fprintf('    absolute Frobenius error  norm(Y - X, ''fro'') :\n');
fprintf('    %-8s', 'mod\demod');
fprintf('%-12s', demodList{:});
fprintf('\n');
for i = 1:numel(modList)
    fprintf('    %-8s', modList{i});
    fprintf('%-12.2e', errAbs(i,:));
    fprintf('\n');
end

results.name     = 'Sim-1 IDZT/DZT perfect reconstruction';
results.err_abs  = max(errAbs(:));
results.err_rel  = max(errRel(:));
results.tol      = 1e-12;
results.pass     = results.err_abs < results.tol;
results.table    = errAbs;
results.modList  = modList;
results.demodList= demodList;
results.Es       = info.Es;

fprintf('\n    max |Y-X|_F = %.3e   (tol %.0e)   -> %s\n', ...
        results.err_abs, results.tol, ternary(results.pass, 'PASS', 'FAIL'));

% ---- optional figure --------------------------------------------------
if doPlot
    s = otfs_modulate(X, 'idzt', P);
    Y = otfs_demodulate(s, M, N, 'dzt', P);

    figure('Name', 'Sim-1 perfect reconstruction', 'Color', 'w');
    subplot(1,3,1); imagesc(0:N-1, 0:M-1, abs(X)); axis square; colorbar;
    xlabel('Doppler n'); ylabel('delay m'); title('|X| (transmitted DD)');
    subplot(1,3,2); imagesc(0:N-1, 0:M-1, abs(Y)); axis square; colorbar;
    xlabel('Doppler n'); ylabel('delay m'); title('|Y| (received DD)');
    subplot(1,3,3); imagesc(0:N-1, 0:M-1, abs(Y-X)); axis square; colorbar;
    xlabel('Doppler n'); ylabel('delay m'); title('|Y - X|  (~1e-16)');
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
