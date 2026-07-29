function results = sim8_delay_time_domain(p, doPlot)
% =========================================================================
% Program : sim8_delay_time_domain.m
% Description :
%   Sim-8  "延遲-時間域 I/O 關係 (4.4.3)"  (chap4_list.md Sec. D, Stage 2)
%
%   Purpose  : verify the delay-time domain channel relation, and the key
%              structural fact that makes low-complexity OTFS detection
%              possible: each sub-block K̃_{m,l} of H̃ is DIAGONAL — for the
%              delay bins where that literally holds.
%   Theory   : ỹ = H̃ x̃,  H̃ = P^T G P  (corrected, see build_channel_matrices.m),
%              x̃ = P^T s,  ỹ = P^T r.
%              H̃ is partitioned into M blocks of size N (block index = delay
%              m). K̃_{m,l} = block (m, [m-l]_M) of H̃, with entries
%                  nu_tilde_{m,l}[n] = g^s[l, m + n*M],   n = 0...N-1   (4.73)
%
%   *** Boundary refinement found while verifying this Sim (not in the
%   original chap4_list.md text, but forced by the numbers — RZP has NO
%   sample from "before the frame", so it has to come from somewhere) ***
%   Tracing q = m+nM through H̃ = P^T G P shows nu_tilde_{m,l}[n] sits at
%   intra-block position (n,n) ONLY when m >= l (no frame-boundary
%   crossing). When m < l, delay tap l reaches back far enough that the
%   needed input sample belongs to the PREVIOUS time-block (n-1), so the
%   nonzero entry sits at (n, n-1) instead — K̃_{m,l} is diagonal for m>=l,
%   but a *shifted* diagonal for m<l (and, under RZP, entirely zero at
%   n=0, since there is no "block -1"). This is the exact same physical
%   effect as the m<l_i branch in theory_dd_relation.m (4.121) — already
%   validated to 1e-15 in Sim-10, so it is trusted here too.
%
%   Method   : build H̃, verify ỹ == H̃ x̃; for every (m,l) with l=0...l_max,
%              extract K̃_{m,l} and compare it against the FULL analytic
%              structure above (not just its diagonal) — this is a complete
%              check, not an approximation.
%   Criterion: ||H̃ x̃ - ỹ||/||ỹ|| < 1e-10;  ||K̃_{m,l} - analytic|| < 1e-10
%              for every (m,l), including the shifted-diagonal blocks.
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw imagesc(abs(H̃)) (Fig. 4.11/4.14 style)
% Output : results - struct with .name .err_rel .tol .pass, plus
%              .err_io, .block_err_max, .n_diagonal_blocks, .n_shifted_blocks
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Stage1'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;  N = p.N;

fprintf('\n=== Sim-8 : delay-time domain I/O relation (4.4.3)  (M=%d, N=%d) ===\n', M, N);

[g, l, k] = gen_channel_taps(p);
gs = gen_gs(p, g, l, k);
G  = gen_G(p, gs, 'RZP');

X  = gen_dd_symbols(p);
Pm = gen_perm_matrix(M, N);
s  = otfs_modulate(X, 'idzt', Pm);
r  = G * s;

C = build_channel_matrices(p, G, Pm);
x_tilde      = Pm.' * s;
y_tilde      = Pm.' * r;
y_tilde_pred = C.Htilde * x_tilde;

err_io = relerr(y_tilde_pred, y_tilde);
fprintf('    ||Htilde*x_tilde - y_tilde|| / ||y_tilde|| = %.3e\n', err_io);

% ---- K_tilde_{m,l}: full structural comparison (diagonal or shifted) ----
l_max = max(l);
block_err_max  = 0;
n_diag_blocks  = 0;
n_shift_blocks = 0;
for m = 0:M-1
    for ll = 0:l_max
        mp = mod(m - ll, M);
        Kt = get_block(C.Htilde, N, m, mp);

        Expected = zeros(N, N);
        isDiagonal = (m >= ll);
        for n = 0:N-1
            q   = m + n*M;
            val = gs(ll+1, q+1);
            if isDiagonal
                nprime = n;
            else
                nprime = n - 1;
                if nprime < 0
                    continue;   % RZP: no predecessor block -> genuinely 0
                end
            end
            Expected(n+1, nprime+1) = val;
        end

        e = relerr(Kt, Expected);
        block_err_max = max(block_err_max, e);
        if isDiagonal, n_diag_blocks = n_diag_blocks + 1; else, n_shift_blocks = n_shift_blocks + 1; end
    end
end

fprintf('    blocks with m>=l (plain diagonal)      : %d\n', n_diag_blocks);
fprintf('    blocks with m<l  (shifted to n-1, RZP)  : %d\n', n_shift_blocks);
fprintf('    max relative error, K~_{m,l} vs full analytic structure : %.3e\n', block_err_max);

results.name           = 'Sim-8 delay-time domain I/O relation (4.4.3)';
results.err_io         = err_io;
results.block_err_max  = block_err_max;
results.n_diag_blocks  = n_diag_blocks;
results.n_shift_blocks = n_shift_blocks;
results.err_rel        = max(err_io, block_err_max);
results.tol            = 1e-10;
results.pass           = (err_io < 1e-10) && (block_err_max < 1e-10);

fprintf('    -> %s\n', ternary(results.pass, 'PASS', 'FAIL'));

if doPlot
    figure('Name', 'Sim-8 |Htilde|', 'Color', 'w');
    imagesc(abs(C.Htilde)); axis square; colorbar; hold on;
    for bIdx = 0:M
        plot([0.5 M*N+0.5], [bIdx*N+0.5 bIdx*N+0.5], 'w-', 'LineWidth', 0.5);
        plot([bIdx*N+0.5 bIdx*N+0.5], [0.5 M*N+0.5], 'w-', 'LineWidth', 0.5);
    end
    title('|$\tilde{H}$|  (delay-time domain, Fig. 4.11/4.14 style, RZP corner note)', 'Interpreter', 'latex');
    xlabel('column index'); ylabel('row index');
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
