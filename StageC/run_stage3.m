function T = run_stage3(p, doPlot)
% =========================================================================
% Program : run_stage3.m
% Description :
%   Stage-3 driver: run Sim-14 ... Sim-16 of chap4_list.md
%   ("階段三：off-grid（分數延遲／分數都卜勒）"), then print the summary
%   error table.
%
%   Sim-14  分數都卜勒的洩漏驗證        numeric DFT of nu_tilde == closed-form D_N(.) (4.79/4.80)
%   Sim-15  分數延遲的洩漏驗證          sinc-interpolated G self-consistent + reduces to Stage-2 on-grid
%   Sim-16  通道稀疏度 S 分析           S vs. N, S vs. fractional degree (informational)
%
%   Sim-14 has an explicit "<1e-10" criterion in chap4_list.md; Sim-15/16
%   do not, so this driver gates Sim-15 on its two internal-consistency
%   checks (which ARE well-defined machine-precision criteria even though
%   chap4_list.md doesn't give explicit numbers) and reports Sim-16 purely
%   as an observation (see sim16_sparsity_analysis.m header for the scope
%   decisions this required).
%
% Usage :
%   T = run_stage3();                    % M=8, N=6, no figures
%   T = run_stage3(otfs_params(), true);  % with all figures
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw all figures              (default false)
% Output : T      - summary table (or struct array if `table` unavailable)
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Stage1'));
addpath(fullfile(thisDir, '..', 'Stage2'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

fprintf('\n');
fprintf('#########################################################################\n');
fprintf('#  Chapter 4  Stage 3 : off-grid (fractional delay / Doppler) taps    #\n');
fprintf('#########################################################################\n');
fprintf('  M = %d (delay bins)    N = %d (Doppler bins)    NM = %d\n', p.M, p.N, p.NM);

R = cell(1, 3);
R{1} = sim14_fractional_doppler_leakage(p, doPlot);
R{2} = sim15_fractional_delay_leakage(p, doPlot);
R{3} = sim16_sparsity_analysis(p, doPlot);

n    = numel(R);
Name = cell(n,1);  Err = zeros(n,1);  Tol = zeros(n,1);  Pass = cell(n,1);
for i = 1:n
    Name{i} = R{i}.name;
    Err(i)  = R{i}.err_rel;
    Tol(i)  = R{i}.tol;
    if R{i}.pass, Pass{i} = 'PASS'; else, Pass{i} = 'FAIL'; end
end

fprintf('\n');
fprintf('=========================================================================\n');
fprintf('  Stage-3 summary\n');
fprintf('=========================================================================\n');
fprintf('  %-70s %11s  %6s\n', 'test', 'rel. error', 'result');
fprintf('  %s\n', repmat('-', 1, 91));
for i = 1:n
    fprintf('  %-70s %11.3e  %6s\n', Name{i}, Err(i), Pass{i});
end
fprintf('  %s\n', repmat('-', 1, 91));

allPass = all(strcmp(Pass, 'PASS'));
if allPass
    fprintf('  ALL GATED STAGE-3 TESTS PASSED (Sim-16 is informational, no gate).\n');
    fprintf('  Off-grid fractional-tap machinery matches the closed-form leakage\n');
    fprintf('  theory and is consistent with Stage 2''s on-grid results.\n');
else
    fprintf('  *** SOME GATED TESTS FAILED - inspect before using off-grid results. ***\n');
end
fprintf('=========================================================================\n\n');

if exist('table', 'file') == 2 || exist('table', 'builtin') == 5
    T = table(Name, Err, Tol, Pass, 'VariableNames', {'Test','RelError','Tol','Result'});
else
    T = struct('Test', Name, 'RelError', num2cell(Err), 'Tol', num2cell(Tol), 'Result', Pass);
end

end
