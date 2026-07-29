function results = sim16_sparsity_analysis(p, doPlot)
% =========================================================================
% Program : sim16_sparsity_analysis.m
% Description :
%   Sim-16  "通道稀疏度 S 分析"  (chap4_list.md Sec. D, Stage 3)
%
%   Purpose  : quantify how sparse the delay-Doppler channel matrix H is —
%              the property that underlies the low-complexity detectors of
%              Chapter 6 (chap4_list.md: "這一項可直接支撐後面第6章偵測
%              複雜度的討論").
%   Definition: S(row) = effective_path_count(|H(row,:)|^2, 0.99), i.e. how
%              many entries of that row are needed to capture 99% of its
%              energy; S = mean(S(row)) over all rows.
%
%   SCOPE NOTE (chap4_list.md gives no formula here, only two "觀察"
%   bullets — this script makes explicit, defensible choices for both):
%     - "S vs. N": M is held fixed and N is swept; the SAME fractional
%       taps (l_i, k_i values, not re-derived from a physical velocity) are
%       reused at every N, so this isolates how S scales with frame
%       dimension NM alone. The reported quantity that matters for
%       detection complexity is S relative to NM (S/NM), not S in
%       isolation — see chap4_list.md Sec. E.1 code 13/14 (LMMSE
%       complexity scales with matrix size).
%     - "S vs. 分數程度": at fixed N, sweep the fractional offset of the
%       taps and show S growing — this reproduces, via the FULL H matrix
%       (both delay AND Doppler leakage), the same qualitative trend Sim-15
%       already showed from the delay-only sinc coefficients alone (a
%       useful cross-check between the two).
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw S-vs-N and S-vs-fractional-offset curves
% Output : results - struct with .name .err_rel .tol .pass (no strict
%              gate; reported informationally per chap4_list.md), plus
%              .Nlist, .NMlist, .avgS_vsN, .fracList, .avgS_vsFrac
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Stage1'));
addpath(fullfile(thisDir, '..', 'Stage2'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;

fprintf('\n=== Sim-16 : channel sparsity S analysis  (M=%d) ===\n', M);

l_taps = [0.3, 1.7, 2.5];
k_taps = [0.2, -1.5, 2.5];

% ---- S vs. N (fixed fractional taps, sweep frame length) ---------------
Nlist   = [6 12 24];
avgS_vsN = zeros(size(Nlist));
NMlist   = zeros(size(Nlist));
for idx = 1:numel(Nlist)
    Ni = Nlist(idx);
    pi_ = otfs_params('M', M, 'N', Ni);
    [g, l, k] = gen_channel_taps(pi_, 'l', l_taps, 'k', k_taps, 'randomPhase', false, 'grid', 'off');
    Lrange = ceil(max(l)) + 9;
    gs = gen_gs_frac(pi_, g, l, k, Lrange);
    G  = gen_G(pi_, gs, 'RZP');
    Pm = gen_perm_matrix(M, Ni);
    Cm = build_channel_matrices(pi_, G, Pm);

    NMi = pi_.NM;
    rowS = zeros(1, NMi);
    for row = 1:NMi
        rowS(row) = effective_path_count(abs(Cm.H(row, :)).^2, 0.99);
    end
    avgS_vsN(idx) = mean(rowS);
    NMlist(idx)   = NMi;
end
fprintf('    N sweep        = %s\n', mat2str(Nlist));
fprintf('    NM             = %s\n', mat2str(NMlist));
fprintf('    mean S         = %s\n', mat2str(avgS_vsN, 4));
fprintf('    mean S / NM    = %s   (relative sparsity, -> 0 supports low-complexity detection)\n', ...
        mat2str(avgS_vsN ./ NMlist, 4));

% ---- S vs. fractional offset (fixed N, full H matrix) -------------------
Nfix = p.N;
p_fix = otfs_params('M', M, 'N', Nfix);
fracList = 0:0.1:0.5;
avgS_vsFrac = zeros(size(fracList));
for idx = 1:numel(fracList)
    off = fracList(idx);
    l_test = [0, 1, 2] + off;
    k_test = [0, 1, -2] + off;
    [g, l, k] = gen_channel_taps(p_fix, 'l', l_test, 'k', k_test, 'randomPhase', false, 'grid', 'off');
    Lrange = ceil(max(l)) + 9;
    gs = gen_gs_frac(p_fix, g, l, k, Lrange);
    G  = gen_G(p_fix, gs, 'RZP');
    Pm = gen_perm_matrix(M, Nfix);
    Cm = build_channel_matrices(p_fix, G, Pm);

    rowS = zeros(1, p_fix.NM);
    for row = 1:p_fix.NM
        rowS(row) = effective_path_count(abs(Cm.H(row, :)).^2, 0.99);
    end
    avgS_vsFrac(idx) = mean(rowS);
end
fprintf('    fractional offset sweep = %s\n', mat2str(fracList));
fprintf('    mean S (fixed N=%d)     = %s\n', Nfix, mat2str(avgS_vsFrac, 4));

results.name        = 'Sim-16 channel sparsity S analysis (informational, no numeric threshold in chap4_list.md)';
results.Nlist       = Nlist;
results.NMlist      = NMlist;
results.avgS_vsN    = avgS_vsN;
results.fracList    = fracList;
results.avgS_vsFrac = avgS_vsFrac;
results.err_rel     = 0;
results.tol         = Inf;
results.pass        = true;     % purely observational, chap4_list.md gives no criterion

fprintf('    -> observational (no pass/fail criterion given for this Sim)\n');

if doPlot
    figure('Name', 'Sim-16 sparsity analysis', 'Color', 'w');
    subplot(1,2,1);
    yyaxis left;  plot(Nlist, avgS_vsN, '-o', 'LineWidth', 1.5); ylabel('mean S');
    yyaxis right; plot(Nlist, avgS_vsN ./ NMlist, '-s', 'LineWidth', 1.5); ylabel('mean S / NM');
    xlabel('N'); grid on; title('Sim-16 : S vs. frame length N');

    subplot(1,2,2);
    plot(fracList, avgS_vsFrac, '-s', 'LineWidth', 1.5); grid on;
    xlabel('fractional offset of taps'); ylabel('mean S');
    title(sprintf('Sim-16 : S vs. fractional degree  (N=%d fixed)', Nfix));
end

end
