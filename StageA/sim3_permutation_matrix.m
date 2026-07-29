function results = sim3_permutation_matrix(p, doPlot)
% =========================================================================
% Program : sim3_permutation_matrix.m
% Description :
%   Sim-3  "Correctness of the permutation matrix P" (chap4_list.md, Stage 1)
%
%   Purpose  : P is the ONLY thing that reconciles the two vectorization
%              directions (Sec. H.2):
%                  s[m + n*M] : delay   index runs fastest
%                  x[n + m*N] : Doppler index runs fastest
%              A wrong P (or a transposed P, see Sec. E.3.2) silently breaks
%              every later domain identity, so it is checked on its own here.
%
%   Theory   : (4.33)(4.34)   P[m+nM , n+mN] = 1
%              (orthogonality)  P.' * P = I_NM,  P.' = inv(P)
%              (4.32)  s        = P * (I_M kron F_N') * x
%              (4.53)  x_tilde  = P.' * s     -> x_tilde[n+mN] = X_tilde[m,n]
%
%   Method   : 1. build P from the book blocks (4.33)(4.34);
%              2. build P again from the intuitive row-column interleaver
%                 ("write column-wise into an N x M array, read out row-wise");
%              3. compare the two, check orthogonality, check (4.32)/(4.53);
%              4. also verify P is a genuine permutation (one 1 per row/col).
%   Criterion: the two constructions are IDENTICAL (error exactly 0) and
%              norm(P.'*P - I_NM) < 1e-12.
%
%   Figures  : spy(P) and the Fig. 4.10 interleaving illustration (M=3, N=4).
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw spy(P) + Fig. 4.10       (default false)
% Output : results - struct with .name .err_rel .tol .pass and sub-errors
% =========================================================================

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;  N = p.N;  NM = N*M;
F_N = p.F_N;

fprintf('\n=== Sim-3 : permutation matrix P (M=%d, N=%d, NM=%d) ===\n', M, N, NM);

% ---------------------------------------------------------------------
% 1) two independent constructions
% ---------------------------------------------------------------------
[P, P_intuitive] = gen_perm_matrix(M, N);

e = struct();
e.two_builds = norm(P - P_intuitive, 'fro');          % must be exactly 0

% ---------------------------------------------------------------------
% 2) orthogonality / permutation-ness
% ---------------------------------------------------------------------
e.orthogonal  = norm(P.'*P - eye(NM), 'fro');         % (Sim-3 criterion)
e.orthogonal2 = norm(P*P.'  - eye(NM), 'fro');
e.is_binary   = norm(P(:) - double(P(:) > 0.5));      % entries are 0/1
e.row_sum     = norm(sum(P, 2) - 1);                  % exactly one 1 per row
e.col_sum     = norm(sum(P, 1).' - 1);                % exactly one 1 per col

% ---------------------------------------------------------------------
% 3) explicit index map  P[m+nM , n+mN] = 1        (4.33)(4.34)
% ---------------------------------------------------------------------
P_idx = zeros(NM, NM);
for m = 0:M-1
    for n = 0:N-1
        P_idx(m + n*M + 1, n + m*N + 1) = 1;
    end
end
e.index_map = norm(P - P_idx, 'fro');

% ---------------------------------------------------------------------
% 4) the identities P is actually used for
% ---------------------------------------------------------------------
[X, x] = gen_dd_symbols(p);

s_ref  = otfs_modulate(X, 'idzt', P);                 % s = vec(X * F_N')
s_vec  = P * kron(eye(M), F_N') * x;                  % (4.32)
e.eq432 = relerr(s_vec, s_ref);

x_tilde = P.' * s_ref;                                % (4.53)
X_tilde = X * F_N';                                   % (4.17)
e.eq453 = relerr(x_tilde, reshape(X_tilde.', NM, 1)); % x_tilde[n+mN] = X_tilde[m,n]

% s and x really are re-orderings of the same numbers
e.same_multiset = norm(sort(abs(s_ref)) - sort(abs(P.'*s_ref)));

fprintf('    build 1 (4.33)(4.34) vs build 2 (interleaver) : %.3e   (must be 0)\n', e.two_builds);
fprintf('    ||P''*P - I_NM||                               : %.3e\n', e.orthogonal);
fprintf('    ||P*P'' - I_NM||                               : %.3e\n', e.orthogonal2);
fprintf('    entries in {0,1} / row sums / col sums        : %.3e / %.3e / %.3e\n', ...
        e.is_binary, e.row_sum, e.col_sum);
fprintf('    index map P[m+nM, n+mN] = 1                   : %.3e\n', e.index_map);
fprintf('    (4.32)  s = P (I_M kron F_N'') x               : %.3e\n', e.eq432);
fprintf('    (4.53)  x_tilde = P'' s                        : %.3e\n', e.eq453);

results.name    = 'Sim-3 permutation matrix P';
results.err     = e;
results.err_rel = max(struct2array_local(e));
results.tol     = 1e-12;
results.pass    = results.err_rel < results.tol;
results.P       = P;

fprintf('\n    max error = %.3e   (tol %.0e)   -> %s\n', ...
        results.err_rel, results.tol, ternary(results.pass, 'PASS', 'FAIL'));

% ---------------------------------------------------------------------
% 5) figures
% ---------------------------------------------------------------------
if doPlot
    % ---- spy(P) for the simulation size ----------------------------
    figure('Name', 'Sim-3 spy(P)', 'Color', 'w');
    spy(P, 12);
    title(sprintf('spy(P),  M = %d, N = %d,  NM = %d', M, N, NM));
    xlabel(sprintf('column index  n + mN   (Doppler fast)\nnz = %d', nnz(P)));
    ylabel('row index  m + nM   (delay fast)');
    axis square;

    % ---- Fig. 4.10 style interleaving illustration, M = 3, N = 4 ----
    plot_interleaver_demo(3, 4);
end

end


% =========================================================================
function plot_interleaver_demo(M, N)
% Reproduce the row-column interleaving picture of Fig. 4.10.
%   x is read out of the DD grid with the DOPPLER index running fastest;
%   the interleaver writes those NM numbers column-wise into an N x M array
%   and reads them out row-wise, giving s with the DELAY index fastest.
NM = N*M;

% label the DD symbols 1..NM in the order they appear in x = vec(X.')
lbl_x = reshape(1:NM, N, M);         % lbl_x(n+1, m+1) = position of that
                                     % symbol inside the vector x
lbl_s = lbl_x.';                     % read row-wise -> layout of s (M x N)

figure('Name', 'Sim-3 Fig.4.10 row-column interleaver', 'Color', 'w');

subplot(1,2,1);
draw_grid(lbl_x, 'write column-wise:  x[n+mN] = X[m,n]', ...
          'block m  (delay)', 'index n  (Doppler)');

subplot(1,2,2);
draw_grid(lbl_s, 'read row-wise:  s[m+nM]', ...
          'block n  (time)', 'index m  (delay)');

% sgtitle needs R2018b; fall back silently on older releases
if exist('sgtitle', 'file')
    sgtitle(sprintf('Fig. 4.10 row-column interleaver  (M = %d, N = %d)', M, N));
end
end


% =========================================================================
function draw_grid(A, ttl, xlab, ylab)
[nr, nc] = size(A);
imagesc(A); colormap(gca, parula); axis equal tight;
set(gca, 'XTick', 1:nc, 'XTickLabel', 0:nc-1, ...
         'YTick', 1:nr, 'YTickLabel', 0:nr-1);
for i = 1:nr
    for j = 1:nc
        text(j, i, sprintf('%d', A(i,j)), 'HorizontalAlignment', 'center', ...
             'Color', 'w', 'FontWeight', 'bold');
    end
end
title(ttl); xlabel(xlab); ylabel(ylab);
end


% =========================================================================
function v = struct2array_local(s)
f = fieldnames(s);
v = zeros(numel(f), 1);
for i = 1:numel(f)
    v(i) = s.(f{i});
end
end


% =========================================================================
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
