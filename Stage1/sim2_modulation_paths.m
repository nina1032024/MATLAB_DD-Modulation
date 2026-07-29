function results = sim2_modulation_paths(p, doPlot)
% =========================================================================
% Program : sim2_modulation_paths.m
% Description :
%   Sim-2  "Two modulation paths are equivalent" (chap4_list.md Sec. D, Stage 1)
%
%   Purpose  : verify that the ISFFT+Heisenberg description (Fig. 4.1) and
%              the IDZT description (Fig. 4.3) generate the SAME transmit
%              signal s, and likewise for the three receiver descriptions.
%   Theory   : because F_M' * F_M = I_M,
%                  Path A: X_tf = F_M X F_N' ; X_tilde = F_M' X_tf   (4.16)(4.17)
%                  Path B: X_tilde = X F_N'                          (4.17)
%                  Path C: element-wise (4.22)
%                  Path D: s = P (I_M kron F_N') vec(X.')            (4.32)
%              are identical, and on the receiver side (4.24)/(4.26)/(4.27)/(4.30)
%              likewise.
%   Method   : run all paths on the same X (and the same r) and compare
%              pairwise against Path B / 'dzt' as the reference.
%   Criterion: every pairwise relative error < 1e-12.
%
%   Also verified here (structural identities behind the equivalence):
%       ||F_M' * F_M - I_M||   = 0     unitarity of F_M
%       X_tf  = F_M * X_tilde  = F_M * X * F_N'          (4.16)
%       Y_tf  = F_M * Y_tilde                            (4.24)
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw the |X_tf| / |X_tilde| maps (default false)
% Output : results - struct with .name .err_tx .err_rx .err_rel .tol .pass
%                    .tx (per-path errors) .rx (per-path errors)
% =========================================================================

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;  N = p.N;
F_M = p.F_M;  F_N = p.F_N;

[X, x] = gen_dd_symbols(p);
P = gen_perm_matrix(M, N);

fprintf('\n=== Sim-2 : equivalence of the modulation / demodulation paths (M=%d, N=%d) ===\n', M, N);

% ======================================================================
% (a) transmitter : s_A == s_B == s_C == s_D
% ======================================================================
[s_A, auxA] = otfs_modulate(X, 'isfft', P);   % ISFFT + Heisenberg (Fig. 4.1)
[s_B, auxB] = otfs_modulate(X, 'idzt',  P);   % IDZT              (Fig. 4.3)
[s_C, ~   ] = otfs_modulate(X, 'elem',  P);   % element-wise (4.22)
[s_D, auxD] = otfs_modulate(X, 'vec',   P);   % vector form  (4.32)

e_tx = struct();
e_tx.A_vs_B = relerr(s_A, s_B);
e_tx.C_vs_B = relerr(s_C, s_B);
e_tx.D_vs_B = relerr(s_D, s_B);

% the (4.31) vectorization used inside Path D must match gen_dd_symbols
e_tx.x_vec  = relerr(auxD.x, x);

% cross-check of the intermediate quantities
e_tx.Xtf    = relerr(auxA.X_tf, F_M * X * F_N');       % (4.16)
e_tx.Xtilde = relerr(auxA.X_tilde, auxB.X_tilde);      % (4.17)
e_tx.unitary_FM = norm(F_M'*F_M - eye(M), 'fro');      % why A == B
e_tx.unitary_FN = norm(F_N'*F_N - eye(N), 'fro');

fprintf('\n  (a) transmitter, reference = Path B (IDZT)\n');
fprintf('      Path A (ISFFT+Heisenberg, 4.16-4.18) vs B : %.3e\n', e_tx.A_vs_B);
fprintf('      Path C (element-wise, 4.22)          vs B : %.3e\n', e_tx.C_vs_B);
fprintf('      Path D (vector form, 4.32)           vs B : %.3e\n', e_tx.D_vs_B);
fprintf('      x = vec(X.'') consistency  (4.31)         : %.3e\n', e_tx.x_vec);
fprintf('      X_tf = F_M*X*F_N''         (4.16)         : %.3e\n', e_tx.Xtf);
fprintf('      ||F_M''F_M - I|| , ||F_N''F_N - I||        : %.3e , %.3e\n', ...
        e_tx.unitary_FM, e_tx.unitary_FN);

% ======================================================================
% (b) receiver : Y_A == Y_B == Y_C == Y_D    (using r = s, cf. (4.27)(4.30))
% ======================================================================
r = s_B;                                       % no channel for this test

[Y_A, aY_A] = otfs_demodulate(r, M, N, 'wigner', P);   % (4.24)(4.25)
[Y_B, aY_B] = otfs_demodulate(r, M, N, 'dzt',    P);   % (4.26)
[Y_C, ~   ] = otfs_demodulate(r, M, N, 'elem',   P);   % (4.27)
[Y_D, aY_D] = otfs_demodulate(r, M, N, 'vec',    P);   % (4.30)

e_rx = struct();
e_rx.A_vs_B = relerr(Y_A, Y_B);
e_rx.C_vs_B = relerr(Y_C, Y_B);
e_rx.D_vs_B = relerr(Y_D, Y_B);
e_rx.Ytf    = relerr(aY_A.Y_tf, F_M * aY_B.Y_tilde);   % (4.24)
% y = vec(Y.') must hold for the vector path (dual of 4.31)
e_rx.y_vec  = relerr(aY_D.y, reshape(Y_D.', N*M, 1));

fprintf('\n  (b) receiver, reference = Path B (DZT)\n');
fprintf('      Path A (Wigner+SFFT, 4.24-4.25)      vs B : %.3e\n', e_rx.A_vs_B);
fprintf('      Path C (element-wise, 4.27)          vs B : %.3e\n', e_rx.C_vs_B);
fprintf('      Path D (vector form, 4.30)           vs B : %.3e\n', e_rx.D_vs_B);
fprintf('      Y_tf = F_M*Y_tilde        (4.24)         : %.3e\n', e_rx.Ytf);
fprintf('      y = vec(Y.'') consistency                 : %.3e\n', e_rx.y_vec);

% ======================================================================
results.name    = 'Sim-2 modulation / demodulation path equivalence';
results.tx      = e_tx;
results.rx      = e_rx;
results.err_tx  = max([e_tx.A_vs_B, e_tx.C_vs_B, e_tx.D_vs_B, e_tx.x_vec, e_tx.Xtf, e_tx.Xtilde]);
results.err_rx  = max([e_rx.A_vs_B, e_rx.C_vs_B, e_rx.D_vs_B, e_rx.Ytf,  e_rx.y_vec]);
results.err_rel = max(results.err_tx, results.err_rx);
results.tol     = 1e-12;
results.pass    = results.err_rel < results.tol;

fprintf('\n    max relative error = %.3e   (tol %.0e)   -> %s\n', ...
        results.err_rel, results.tol, ternary(results.pass, 'PASS', 'FAIL'));

% ---- optional figure --------------------------------------------------
if doPlot
    figure('Name', 'Sim-2 domains along the transmit chain', 'Color', 'w');
    subplot(1,3,1); imagesc(0:N-1, 0:M-1, abs(X));            axis square; colorbar;
    xlabel('Doppler n'); ylabel('delay m');  title('|X| : DD domain');
    subplot(1,3,2); imagesc(0:N-1, 0:M-1, abs(auxA.X_tf));    axis square; colorbar;
    xlabel('time k');    ylabel('freq l');   title('|X_{tf}| : ISFFT (4.16)');
    subplot(1,3,3); imagesc(0:N-1, 0:M-1, abs(auxB.X_tilde)); axis square; colorbar;
    xlabel('time n');    ylabel('delay m');  title('|$\tilde{X}$| : IDZT (4.17)', 'Interpreter','latex');
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
