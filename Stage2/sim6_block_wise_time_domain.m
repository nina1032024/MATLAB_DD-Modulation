function results = sim6_block_wise_time_domain(p, doPlot)
% =========================================================================
% Program : sim6_block_wise_time_domain.m
% Description :
%   Sim-6  "Block-wise 時域關係 (4.39)"  (chap4_list.md Sec. D, Stage 2)
%
%   Purpose  : verify the block-partitioned view of r = G*s used throughout
%              Sec. 4.4: because l_max < M, a delay tap can shift a sample
%              across AT MOST one M-block boundary, so
%                  r_0 = G_{0,0} s_0
%                  r_n = G_{n,0} s_n + G_{n,1} s_{n-1},   n = 1...N-1   (4.39)
%              where G partitioned into N blocks of size M (block index =
%              time slot), G_{n,0} = self block (row n, col n), G_{n,1} =
%              one-block-back block (row n, col n-1).
%   Method   : slice G (RZP) into its N x N grid of M x M blocks; verify
%              every block OTHER than the diagonal and first sub-diagonal
%              is (numerically) zero; reconstruct r block-by-block via
%              (4.39) and compare against the full G*s.
%   Criterion: other-block Frobenius norm < 1e-12; block-wise r == G*s,
%              relative error < 1e-12.
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw the block-structure figure (default false)
% Output : results - struct with .name .err_rel .tol .pass, plus .maxOtherBlockNorm
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Stage1'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;  N = p.N;

fprintf('\n=== Sim-6 : block-wise time-domain relation (4.39)  (M=%d, N=%d) ===\n', M, N);

[g, l, k] = gen_channel_taps(p);
gs = gen_gs(p, g, l, k);
G  = gen_G(p, gs, 'RZP');

X  = gen_dd_symbols(p);
Pm = gen_perm_matrix(M, N);
s  = otfs_modulate(X, 'idzt', Pm);
r_full = G * s;

% ---- (a) all blocks other than diagonal / first sub-diagonal are ~0 ----
maxOtherBlockNorm = 0;
for n = 0:N-1
    for n2 = 0:N-1
        if n2 == n || n2 == n-1
            continue;
        end
        Bn = get_block(G, M, n, n2);
        maxOtherBlockNorm = max(maxOtherBlockNorm, norm(Bn, 'fro'));
    end
end

% ---- (b) reconstruct r block-by-block via (4.39) -----------------------
s_blocks = reshape(s, M, N);          % column n+1 = s_n  (q = m + n*M)
r_block  = zeros(M, N);
for n = 0:N-1
    Gn0 = get_block(G, M, n, n);      % self block
    r_n = Gn0 * s_blocks(:, n+1);
    if n >= 1
        Gn1 = get_block(G, M, n, n-1);      % one-block-back
        r_n = r_n + Gn1 * s_blocks(:, n);   % s_blocks(:,n) = s_{n-1}
    end
    r_block(:, n+1) = r_n;
end
r_block_vec = r_block(:);

err_blockwise = relerr(r_block_vec, r_full);

fprintf('    max Frobenius norm of non-(diag / 1st-subdiag) blocks : %.3e\n', maxOtherBlockNorm);
fprintf('    ||r_blockwise - G*s|| / ||G*s||                       : %.3e\n', err_blockwise);

results.name              = 'Sim-6 block-wise time-domain relation (4.39)';
results.maxOtherBlockNorm = maxOtherBlockNorm;
results.err_rel           = max(maxOtherBlockNorm, err_blockwise);
results.tol               = 1e-12;
results.pass              = results.err_rel < results.tol;

fprintf('    -> %s\n', ternary(results.pass, 'PASS', 'FAIL'));

if doPlot
    figure('Name', 'Sim-6 block structure of G', 'Color', 'w');
    imagesc(abs(G)); axis square; colorbar; hold on;
    for bIdx = 0:N
        plot([0.5 N*M+0.5], [bIdx*M+0.5 bIdx*M+0.5], 'w-', 'LineWidth', 0.5);
        plot([bIdx*M+0.5 bIdx*M+0.5], [0.5 N*M+0.5], 'w-', 'LineWidth', 0.5);
    end
    title('Sim-6 : only diagonal (G_{n,0}) and first sub-diagonal (G_{n,1}) blocks are nonzero');
    xlabel('column q'''); ylabel('row q');
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
