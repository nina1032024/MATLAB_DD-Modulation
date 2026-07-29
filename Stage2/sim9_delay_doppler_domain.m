function results = sim9_delay_doppler_domain(p, doPlot)
% =========================================================================
% Program : sim9_delay_doppler_domain.m
% Description :
%   Sim-9  "延遲-都卜勒域 I/O 關係 (4.4.4)"  (chap4_list.md Sec. D, Stage 2)
%
%   Purpose  : verify the final delay-Doppler domain channel relation, and
%              that each sub-block K_{m,l} of H is CIRCULANT (the DFT
%              similarity transform F_N Λ F_N^† of a DIAGONAL matrix Λ is
%              always circulant — the fact behind OTFS's sparse, structured
%              DD-domain channel).
%   Theory   : y = H x,  H = (I_M kron F_N) H̃ (I_M kron F_N)^†,
%              K_{m,l} = F_N K̃_{m,l} F_N^†                              (4.60)
%              nu_{m,l}[k] = (1/N) sum_n nu_tilde_{m,l}[n] e^{-j2*pi*k*n/N} (4.77)
%
%   *** Boundary refinement, following directly from Sim-8 ***
%   K̃_{m,l} is a PLAIN diagonal only for m>=l (Sim-8); for m<l it is a
%   *shifted* diagonal (= a cyclic-shift permutation times a diagonal).
%   F_N * (shift * diag) * F_N^† = Λ_shift * (F_N * diag * F_N^†): a
%   diagonal-times-circulant product, which is generally NOT itself
%   circulant (unless Λ_shift were constant, which it isn't for a
%   nontrivial shift). So the "K_{m,l} circulant" claim is only
%   unconditionally true for the m>=l blocks. This mirrors the same m<l
%   boundary effect already found in Sim-8 and Sim-10/11 (loss of the
%   textbook's simplified per-block symmetry right at the frame edge).
%
%   Method   : y=(I_M⊗F_N)*(P^T r), compare against H*x AND against
%              otfs_demodulate(r,...,'dzt') (three independent routes to Y);
%              for every (m,l): if m>=l, verify K_{m,l} is circulant
%              (gated, hard criterion); if m<l, just report how far it sits
%              from circulant (informational, not gated — expected > 0).
%   Criterion: I/O relative errors < 1e-10; circulant relative error < 1e-10
%              for all "safe" (m>=l) blocks.
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw imagesc(abs(H)) (Fig. 4.12/4.15 style)
% Output : results - struct with .name .err_rel .tol .pass, plus
%              .err_io, .err_demod, .err_analysis,
%              .circ_err_max_safe, .circ_err_max_boundary
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Stage1'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;  N = p.N;

fprintf('\n=== Sim-9 : delay-Doppler domain I/O relation (4.4.4)  (M=%d, N=%d) ===\n', M, N);

[g, l, k] = gen_channel_taps(p);
gs = gen_gs(p, g, l, k);
G  = gen_G(p, gs, 'RZP');

[X, x] = gen_dd_symbols(p);
Pm = gen_perm_matrix(M, N);
s  = otfs_modulate(X, 'idzt', Pm);
r  = G * s;

C = build_channel_matrices(p, G, Pm);

% analysis operator A = (I_M kron F_N) P^T must itself recover x from s
% (channel-free identity, cf. Sim-1/Sim-4) -- a useful internal sanity check
x_from_s = C.IM_FN * (Pm.' * s);
err_analysis = relerr(x_from_s, x);

y      = C.IM_FN * (Pm.' * r);
y_pred = C.H * x;
err_io = relerr(y_pred, y);

Y_direct     = otfs_demodulate(r, M, N, 'dzt', Pm);
y_direct_vec = reshape(Y_direct.', N*M, 1);
err_demod    = relerr(y, y_direct_vec);

fprintf('    analysis operator recovers x from s (no channel)  : %.3e\n', err_analysis);
fprintf('    ||H*x - y|| / ||y||                                : %.3e\n', err_io);
fprintf('    ||y - vec(Y_demod.'')|| / ||y||                     : %.3e\n', err_demod);

% ---- K_{m,l} circulant check, split safe (m>=l) vs boundary (m<l) -------
l_max = max(l);
circ_err_max_safe     = 0;
circ_err_max_boundary = 0;
n_safe = 0; n_boundary = 0;
for m = 0:M-1
    for ll = 0:l_max
        mp = mod(m - ll, M);
        K = get_block(C.H, N, m, mp);
        c = K(:, 1);
        Kcirc = zeros(N, N);
        for col = 0:N-1
            Kcirc(:, col+1) = circshift(c, col);
        end
        e = norm(K - Kcirc, 'fro') / max(norm(K, 'fro'), eps);
        if m >= ll
            circ_err_max_safe = max(circ_err_max_safe, e);
            n_safe = n_safe + 1;
        else
            circ_err_max_boundary = max(circ_err_max_boundary, e);
            n_boundary = n_boundary + 1;
        end
    end
end
fprintf('    "safe" blocks (m>=l, %d blocks)      : max circulant rel. error = %.3e  (gated)\n', n_safe, circ_err_max_safe);
fprintf('    boundary blocks (m<l, %d blocks)      : max circulant rel. error = %.3e  (informational, expected > 0)\n', n_boundary, circ_err_max_boundary);

results.name               = 'Sim-9 delay-Doppler domain I/O relation (4.4.4)';
results.err_analysis       = err_analysis;
results.err_io             = err_io;
results.err_demod          = err_demod;
results.circ_err_max_safe     = circ_err_max_safe;
results.circ_err_max_boundary = circ_err_max_boundary;
results.err_rel      = max([err_analysis, err_io, err_demod, circ_err_max_safe]);
results.tol          = 1e-10;
results.pass         = results.err_rel < results.tol;

fprintf('    -> %s\n', ternary(results.pass, 'PASS', 'FAIL'));

if doPlot
    figure('Name', 'Sim-9 |H|', 'Color', 'w');
    imagesc(abs(C.H)); axis square; colorbar; hold on;
    for bIdx = 0:M
        plot([0.5 M*N+0.5], [bIdx*N+0.5 bIdx*N+0.5], 'w-', 'LineWidth', 0.5);
        plot([bIdx*N+0.5 bIdx*N+0.5], [0.5 M*N+0.5], 'w-', 'LineWidth', 0.5);
    end
    title('|H|  (delay-Doppler domain, Fig. 4.12/4.15 style)');
    xlabel('column index'); ylabel('row index');
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
