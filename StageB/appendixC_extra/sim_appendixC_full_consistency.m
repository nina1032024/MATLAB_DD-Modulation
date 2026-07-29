function results = sim_appendixC_full_consistency(p, doPlot)
% =========================================================================
% Program : sim_appendixC_full_consistency.m
% Description :
%   Supplementary Sim (not numbered Sim-1~16): the CONSOLIDATED "r via 4
%   methods x Y via 3 methods" cross-check that chap4_list.md Sec. E.4
%   explicitly asks for ("你手上就有 r 的 4 種算法 × Y 的 3 種算法 + 1 個
%   理論值，任兩兩相減都應該 < 1e-10。這張誤差表放進進度報告，老師第 5
%   點就完整交代了") — this exact consolidated table was never actually
%   assembled before (Sim-5/8/9 validate each DOMAIN relation separately,
%   but never all 4 r-methods and all 3 Y-methods side by side in one
%   table, using the book's own Appendix-C code structure).
%
%   Reproduces Appendix C, MATLAB code 3 / 8 / 9 / 10 / 12, using our
%   ALREADY-CORRECTED conventions (chap4_list.md Sec. E.3):
%     - code 10 Method 1's off-by-one (`ell=0:(delay_spread-1)`) is NOT
%       reproduced — apply_channel_conv.m (Stage 2) already sums
%       ell=0:delay_spread inclusive.
%     - code 9's `H_tilda=P*G*P.'` typo is NOT reproduced —
%       build_channel_matrices.m (Stage 2) already uses H̃=P^T G P.
%     - all vectors are kept as NM x 1 COLUMN vectors throughout (code 3
%       Method 1's `reshape(X_tilda,1,N*M)` row-vector is not reproduced).
%
%   Two paths are also NEW here (never checked in Sim-2/Sim-7 before):
%     - modulation "Method 3" (code 3, Eq. 4.35): s = kron(Fn',Im)*P*x,
%       a DIFFERENT operator ordering than Method 2's P*kron(Im,Fn')*x.
%     - demodulation "Method 3" (code 12): Y from
%       y=(P.')*kron(Fn,Im)*r, again a different ordering than Method 2's
%       kron(Im,Fn)*(P.')*r.
%   Both rely on the same permutation<->Kronecker commutation identity
%   that makes the row-column interleaver P work in the first place
%   (a "stride permutation" identity familiar from mixed-radix FFT
%   factorizations) — worth checking numerically even though it follows
%   algebraically, since it exercises P in a structurally different way
%   than every other Sim so far.
%
%   NOTE ON INDEPENDENCE: of the 4 r-methods, Method 1 (TDL convolution)
%   and Method 2 (G*s) are built from independent first-principles code
%   (this is what Sim-5 already checks). Methods 3 and 4 (via H̃ and H)
%   are algebraically FORCED to agree with each other and with Method 2,
%   given how H̃ and H were constructed (a short proof is in
%   build_channel_matrices.m's derivation) — so they are regression/
%   construction checks, not independent physical validations. All are
%   still worth consolidating into one table, which is the actual
%   deliverable chap4_list.md Sec. E.4 asks for.
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw the 4x4 / 3x3 error-matrix heatmaps
% Output : results - struct with .name .err_rel .tol .pass, plus
%              .err_mod_new, .err_demod_new, .errR (4x4), .errY (3x3)
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..', 'Stage1'));
addpath(fullfile(thisDir, '..'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;      end

M = p.M;  N = p.N;  F_N = p.F_N;  Im = eye(M);

fprintf('\n=== Appendix C supplementary Sim : consolidated 4-route(r) x 3-route(Y) table ===\n');
fprintf('    (code 3 / 8 / 9 / 10 / 12, chap4_list.md Sec. E.4 "最小可跑驗證腳本")\n');

[g, l, k] = gen_channel_taps(p);
gs = gen_gs(p, g, l, k);
G  = gen_G(p, gs, 'RZP');

[X, x] = gen_dd_symbols(p);
Pm = gen_perm_matrix(M, N);

% ---- modulation: NEW Method-3 ordering (code 3, Eq. 4.35) --------------
s_ref = otfs_modulate(X, 'idzt', Pm);                  % trusted (Sim-1/Sim-2)
s_new = kron(F_N', Im) * (Pm * x);                     % code3 Method3: kron(Fn',Im)*P*x
err_mod_new = relerr(s_new, s_ref);
fprintf('  modulation Method 3 (kron(Fn'',Im)*P*x, NEW) vs Method 1 (idzt) : %.3e\n', err_mod_new);

s = s_ref;

% ---- channel: r via 4 methods (code 8+9+10) ----------------------------
r1 = apply_channel_conv(p, gs, s, 'RZP');               % Method 1: TDL conv (off-by-one fixed)
r2 = G * s;                                             % Method 2: time-domain matrix

Cm = build_channel_matrices(p, G, Pm);                  % H̃ = P^T G P (fixed)
X_tilde = X * F_N';
x_tilde = reshape(X_tilde.', N*M, 1);
y_tilde = Cm.Htilde * x_tilde;
r3 = Pm * y_tilde;                                      % Method 3: via H̃

y  = Cm.H * x;
r4 = Pm * kron(Im, F_N') * y;                           % Method 4: via H

rAll   = {r1, r2, r3, r4};
rNames = {'M1:TDL-conv', 'M2:G*s', 'M3:via-Htilde', 'M4:via-H'};
errR = zeros(4, 4);
for a = 1:4
    for b = 1:4
        errR(a, b) = relerr(rAll{a}, rAll{b});
    end
end

fprintf('\n  r (code 8-10) pairwise relative-error table:\n');
fprintf('  %-16s', ''); fprintf('%-14s', rNames{:}); fprintf('\n');
for a = 1:4
    fprintf('  %-16s', rNames{a});
    fprintf('%-14.2e', errR(a, :));
    fprintf('\n');
end

r = r2;

% ---- demod: NEW Method-3 ordering (code 12) ----------------------------
Y_ref = otfs_demodulate(r, M, N, 'dzt', Pm);            % trusted (Sim-1/Sim-2)
y3    = Pm.' * kron(F_N, Im) * r;                       % code12 Method3: (P.')*kron(Fn,Im)*r
Y_new = reshape(y3, N, M).';
err_demod_new = relerr(Y_new, Y_ref);
fprintf('\n  demod Method 3 ((P.'')*kron(Fn,Im)*r, NEW) vs Method 1 (dzt) : %.3e\n', err_demod_new);

% also assemble the 3-method Y table (Method1=dzt, Method2=vec [Sim-2], Method3=NEW)
y2 = kron(eye(M), F_N) * (Pm.' * r);
Y2 = reshape(y2, N, M).';
YAll   = {Y_ref, Y2, Y_new};
YNames = {'M1:dzt', 'M2:vec(P after)', 'M3:vec(P before)[NEW]'};
errY = zeros(3, 3);
for a = 1:3
    for b = 1:3
        errY(a, b) = relerr(YAll{a}(:), YAll{b}(:));
    end
end

fprintf('\n  Y (code 12) pairwise relative-error table:\n');
fprintf('  %-20s', ''); fprintf('%-16s', YNames{:}); fprintf('\n');
for a = 1:3
    fprintf('  %-20s', YNames{a});
    fprintf('%-16.2e', errY(a, :));
    fprintf('\n');
end

results.name         = 'Appendix C consolidated 4-route(r) x 3-route(Y) consistency table';
results.err_mod_new  = err_mod_new;
results.err_demod_new = err_demod_new;
results.errR         = errR;
results.errY         = errY;
results.rNames       = rNames;
results.YNames       = YNames;
results.err_rel      = max([err_mod_new, err_demod_new, errR(:).', errY(:).']);
results.tol          = 1e-9;
results.pass         = results.err_rel < results.tol;

fprintf('\n  -> max error over ALL modulation/r/Y cross-checks = %.3e   -> %s\n', ...
        results.err_rel, ternary(results.pass, 'PASS', 'FAIL'));

if doPlot
    figure('Name', 'Appendix C consolidated consistency', 'Color', 'w');
    subplot(1,2,1);
    imagesc(log10(errR + eps)); axis square; colorbar;
    set(gca, 'XTick', 1:4, 'XTickLabel', rNames, 'YTick', 1:4, 'YTickLabel', rNames, ...
             'XTickLabelRotation', 20);
    title('log_{10}(rel. error), r methods (code 8-10)');

    subplot(1,2,2);
    imagesc(log10(errY + eps)); axis square; colorbar;
    set(gca, 'XTick', 1:3, 'XTickLabel', YNames, 'YTick', 1:3, 'YTickLabel', YNames, ...
             'XTickLabelRotation', 20);
    title('log_{10}(rel. error), Y methods (code 12)');
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
