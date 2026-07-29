function results = sim4_parseval(p, doPlot)
% =========================================================================
% Program : sim4_parseval.m
% Description :
%   Sim-4  "Parseval / power check" (chap4_list.md Sec. D, Stage 1)
%
%   Purpose  : every transform in the OTFS chain (ISFFT, IDZT, Heisenberg
%              with a rectangular pulse, the permutation P) is UNITARY, so
%              the signal energy must be identical in all four domains.
%              This is an independent way to catch a missing sqrt(M)/sqrt(N)
%              that could hide inside a self-consistent mod/demod pair.
%
%   Theory   : ||x||^2 = ||s||^2 = ||x_check||^2 = ||x_tilde||^2
%              where (Sec. A.4)
%                  x        = vec(X.')                  DD domain     (4.31)
%                  s        = vec(X * F_N')             time domain   (4.18)
%                  x_check  = (I_N kron F_M) * s        time-frequency(4.41)
%                  x_tilde  = P.' * s                   delay-time    (4.53)
%              and additionally, by unitarity of F_M and F_N,
%                  ||X||_F^2 = ||X_tf||_F^2 = ||X_tilde||_F^2.
%              With E_s = 1 normalization all of these equal NM.
%
%   Method   : build every domain representation of the SAME frame and
%              compare energies against ||x||^2.
%   Criterion: relative error < 1e-12.
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw the per-domain energy bar chart (default false)
% Output : results - struct with .name .err_rel .tol .pass .energy (table)
% =========================================================================

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;  N = p.N;  NM = p.NM;
F_M = p.F_M;  F_N = p.F_N;

[X, x] = gen_dd_symbols(p);
P = gen_perm_matrix(M, N);

fprintf('\n=== Sim-4 : Parseval / energy conservation (M=%d, N=%d, NM=%d) ===\n', M, N, NM);

% ---------------------------------------------------------------------
% build every domain representation of the same frame
% ---------------------------------------------------------------------
X_tf    = F_M * X * F_N';                 % (4.16) time-frequency, M x N
X_tilde = F_M' * X_tf;                    % (4.17) delay-time,     M x N
s       = X_tilde(:);                     % (4.18) time domain,    NM x 1
x_check = kron(eye(N), F_M) * s;          % (4.41) time-frequency, NM x 1
x_tilde = P.' * s;                        % (4.53) delay-time,     NM x 1

% receiver-side counterparts (with r = s), to confirm the demodulator too
Y       = otfs_demodulate(s, M, N, 'dzt', P);
y       = reshape(Y.', NM, 1);

% ---------------------------------------------------------------------
% energies
% ---------------------------------------------------------------------
names = { 'x        (DD, 4.31)', ...
          's        (time, 4.18)', ...
          'x_check  (TF, 4.41)', ...
          'x_tilde  (delay-time, 4.53)', ...
          'y        (DD, receiver)', ...
          'X        (DD matrix)', ...
          'X_tf     (TF matrix, 4.16)', ...
          'X_tilde  (delay-time matrix)' };

E = [ norm(x)^2 ;
      norm(s)^2 ;
      norm(x_check)^2 ;
      norm(x_tilde)^2 ;
      norm(y)^2 ;
      norm(X,       'fro')^2 ;
      norm(X_tf,    'fro')^2 ;
      norm(X_tilde, 'fro')^2 ];

E_ref = E(1);                             % reference: ||x||^2
err   = abs(E - E_ref) / E_ref;

fprintf('\n    %-32s %14s %12s\n', 'domain', 'energy', 'rel. err');
for i = 1:numel(names)
    fprintf('    %-32s %14.8f %12.2e\n', names{i}, E(i), err(i));
end
fprintf('\n    theoretical value with E_s = 1 :  NM = %d\n', NM);

% ---------------------------------------------------------------------
% Sanity: the same statement expressed as operator unitarity
% ---------------------------------------------------------------------
e_op = struct();
e_op.FM  = norm(F_M'*F_M   - eye(M),  'fro');
e_op.FN  = norm(F_N'*F_N   - eye(N),  'fro');
e_op.P   = norm(P.'*P      - eye(NM), 'fro');
e_op.INF = norm(kron(eye(N),F_M)' * kron(eye(N),F_M) - eye(NM), 'fro');
e_op.A   = norm( (kron(eye(M),F_N)*P.')' * (kron(eye(M),F_N)*P.') - eye(NM), 'fro');
%           A = (I_M kron F_N) P.'   is the DD-domain analysis operator (4.59)

fprintf('\n    operator unitarity checks\n');
fprintf('      ||F_M''F_M - I||                 : %.3e\n', e_op.FM);
fprintf('      ||F_N''F_N - I||                 : %.3e\n', e_op.FN);
fprintf('      ||P''P - I||                     : %.3e\n', e_op.P);
fprintf('      ||(I_N kron F_M)^H(.) - I||     : %.3e\n', e_op.INF);
fprintf('      ||A^H A - I||,  A = (I_M kron F_N)P''  : %.3e\n', e_op.A);

results.name    = 'Sim-4 Parseval / energy conservation';
results.energy  = E;
results.names   = names;
results.err_rel = max([err(:); e_op.FM; e_op.FN; e_op.P; e_op.INF; e_op.A]);
results.tol     = 1e-12;
results.pass    = results.err_rel < results.tol;
results.E_ref   = E_ref;

fprintf('\n    max relative error = %.3e   (tol %.0e)   -> %s\n', ...
        results.err_rel, results.tol, ternary(results.pass, 'PASS', 'FAIL'));

% ---------------------------------------------------------------------
if doPlot
    figure('Name', 'Sim-4 Parseval', 'Color', 'w');
    bar(E); grid on;
    set(gca, 'XTick', 1:numel(names), 'XTickLabel', names, 'XTickLabelRotation', 30);
    ylabel('energy  ||.||^2');
    ylim([0, 1.25*E_ref]);
    hold on;
    plot([0.5, numel(names)+0.5], [NM NM], 'r--', 'LineWidth', 1.5);
    legend('measured', sprintf('NM = %d (E_s = 1)', NM), 'Location', 'south');
    title('Sim-4 : energy is identical in every domain (unitary transforms)');
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
