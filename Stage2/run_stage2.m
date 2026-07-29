function T = run_stage2(p, doPlot)
% =========================================================================
% Program : run_stage2.m
% Description :
%   Stage-2 driver: run Sim-5 ... Sim-13 of chap4_list.md
%   ("階段二：on-grid 通道，四個 domain 的 I/O 關係驗證"), then print the
%   summary error table.
%
%   Sim-5   時域通道兩種實作一致              G*s        == direct convolution
%   Sim-6   Block-wise 時域關係 (4.39)        r_n        == G_{n,0}s_n+G_{n,1}s_{n-1}
%   Sim-7   時頻域 I/O 關係 (4.4.2)           Hcheck*x_check == y_check ; ICI/ISI
%   Sim-8   延遲-時間域 I/O 關係 (4.4.3)      Htilde*x_tilde == y_tilde ; K~ diagonal
%   Sim-9   延遲-都卜勒域 I/O 關係 (4.4.4)    H*x        == y ; K circulant
%   Sim-10  符號級 DD 域理論式驗證 (最終驗收) Y_theory   == Y_sim (RCP)
%   Sim-11  理想脈波 vs. 矩形脈波比較         (informational + 1 gated case)
%   Sim-12  單路徑位移測試                    delta shift/phase check
%   Sim-13  雜訊統計一致性                    variance / correlation preserved
%
%   All taps are ON-GRID (integer l_i, k_i) per Stage 2's scope; off-grid
%   (fractional) taps are Stage 3 (Sim-14 ... Sim-16), not covered here.
%
% Usage :
%   T = run_stage2();                        % M=8, N=6, no figures
%   T = run_stage2(otfs_params(), true);      % with all figures
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw all figures              (default false)
% Output : T      - summary table (or struct array if `table` unavailable)
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Stage1'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

fprintf('\n');
fprintf('#########################################################################\n');
fprintf('#  Chapter 4  Stage 2 : on-grid channel, four-domain I/O verification  #\n');
fprintf('#########################################################################\n');
fprintf('  M = %d (delay bins)    N = %d (Doppler bins)    NM = %d\n', p.M, p.N, p.NM);
fprintf('  Delta_f = %.3g kHz   T = %.3g us   fc = %.3g GHz\n', p.df/1e3, p.T*1e6, p.fc/1e9);

R = cell(1, 9);
R{1} = sim5_time_domain_convolution(p, doPlot);
R{2} = sim6_block_wise_time_domain(p, doPlot);
R{3} = sim7_time_frequency_domain(p, doPlot);
R{4} = sim8_delay_time_domain(p, doPlot);
R{5} = sim9_delay_doppler_domain(p, doPlot);
R{6} = sim10_symbol_level_theory(p, doPlot);
R{7} = sim11_ideal_vs_rectangular_pulse(p, doPlot);
R{8} = sim12_single_path_shift_test(p, doPlot);
R{9} = sim13_noise_statistics(p, [], doPlot);

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
fprintf('  Stage-2 summary\n');
fprintf('=========================================================================\n');
fprintf('  %-58s %11s  %6s\n', 'test', 'rel. error', 'result');
fprintf('  %s\n', repmat('-', 1, 79));
for i = 1:n
    fprintf('  %-58s %11.3e  %6s\n', Name{i}, Err(i), Pass{i});
end
fprintf('  %s\n', repmat('-', 1, 79));

allPass = all(strcmp(Pass, 'PASS'));
if allPass
    fprintf('  ALL STAGE-2 TESTS PASSED  -> on-grid channel matrices (G, Ȟ, H̃, H) are\n');
    fprintf('  consistent across all four domains and match the symbol-level theory.\n');
    fprintf('  Safe to continue with Stage 3 (Sim-14 ... Sim-16, off-grid / fractional taps).\n');
else
    fprintf('  *** SOME TESTS FAILED - inspect the failing Sim before moving to Stage 3. ***\n');
end
fprintf('=========================================================================\n\n');

if exist('table', 'file') == 2 || exist('table', 'builtin') == 5
    T = table(Name, Err, Tol, Pass, 'VariableNames', {'Test','RelError','Tol','Result'});
else
    T = struct('Test', Name, 'RelError', num2cell(Err), 'Tol', num2cell(Tol), 'Result', Pass);
end

end
