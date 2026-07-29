function T = run_stage1(p, doPlot)
% =========================================================================
% Program : run_stage1.m
% Description :
%   Stage-1 driver: run Sim-1 ... Sim-4 of chap4_list.md
%   ("Stage 1 - self-consistency checks with no channel"), then print the
%   summary error table requested in Sec. F / Sec. G.9 of chap4_list.md.
%
%   Sim-1  IDZT/DZT perfect reconstruction        Y == X
%   Sim-2  equivalence of the modulation paths    s_A == s_B == s_C == s_D
%   Sim-3  correctness of the permutation P       two builds identical, P'P = I
%   Sim-4  Parseval / energy conservation         ||x|| = ||s|| = ||x_check|| = ||x_tilde||
%
%   Everything here is deterministic (rng seeded in gen_dd_symbols), so the
%   table can be pasted straight into the progress report.
%
% Usage :
%   run_stage1                       % M = 8, N = 6, no figures
%   run_stage1(otfs_params(), true)  % with all figures
%   run_stage1(otfs_params('M',64,'N',16))     % larger frame
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw all figures              (default false)
% Output : T      - summary table (or struct array if `table` is unavailable)
% =========================================================================

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

fprintf('\n');
fprintf('#########################################################################\n');
fprintf('#  Chapter 4  Stage 1 : self-consistency checks (no channel, no noise)  #\n');
fprintf('#########################################################################\n');
fprintf('  M = %d (delay bins)    N = %d (Doppler bins)    NM = %d\n', p.M, p.N, p.NM);
fprintf('  Delta_f = %.3g kHz   T = %.3g us   Ts = %.4g us   fc = %.3g GHz\n', ...
        p.df/1e3, p.T*1e6, p.Ts*1e6, p.fc/1e9);
fprintf('  delay resolution = %.4g us   Doppler resolution = %.4g Hz\n', ...
        p.d_res*1e6, p.v_res);
fprintf('  frame duration Tf = %.4g ms   modOrder = %d   seed = %d\n', ...
        p.Tf*1e3, p.modOrder, p.seed);

R = cell(1, 4);
R{1} = sim1_perfect_reconstruction(p, doPlot);
R{2} = sim2_modulation_paths(p, doPlot);
R{3} = sim3_permutation_matrix(p, doPlot);
R{4} = sim4_parseval(p, doPlot);

% ---------------------------------------------------------------------
% summary table
% ---------------------------------------------------------------------
n     = numel(R);
Name  = cell(n,1);
Err   = zeros(n,1);
Tol   = zeros(n,1);
Pass  = cell(n,1);

for i = 1:n
    Name{i} = R{i}.name;
    Err(i)  = R{i}.err_rel;
    Tol(i)  = R{i}.tol;
    if R{i}.pass, Pass{i} = 'PASS'; else, Pass{i} = 'FAIL'; end
end

fprintf('\n');
fprintf('=========================================================================\n');
fprintf('  Stage-1 summary\n');
fprintf('=========================================================================\n');
fprintf('  %-46s %11s  %6s\n', 'test', 'rel. error', 'result');
fprintf('  %s\n', repmat('-', 1, 67));
for i = 1:n
    fprintf('  %-46s %11.3e  %6s\n', Name{i}, Err(i), Pass{i});
end
fprintf('  %s\n', repmat('-', 1, 67));

allPass = all(strcmp(Pass, 'PASS'));
if allPass
    fprintf('  ALL STAGE-1 TESTS PASSED  -> modulator / demodulator are consistent.\n');
    fprintf('  Safe to continue with Stage 2 (Sim-5 ... Sim-13, on-grid channel).\n');
else
    fprintf('  *** SOME TESTS FAILED - fix Stage 1 before touching the channel. ***\n');
end
fprintf('=========================================================================\n\n');

% ---- return as a table when available (R2013b+), else a struct array ---
if exist('table', 'file') == 2 || exist('table', 'builtin') == 5
    T = table(Name, Err, Tol, Pass, 'VariableNames', {'Test','RelError','Tol','Result'});
else
    T = struct('Test', Name, 'RelError', num2cell(Err), ...
               'Tol', num2cell(Tol), 'Result', Pass);
end

end
