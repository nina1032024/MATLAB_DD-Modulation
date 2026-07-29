function results = sim14_fractional_doppler_leakage(p, doPlot)
% =========================================================================
% Program : sim14_fractional_doppler_leakage.m
% Description :
%   Sim-14  "分數都卜勒的洩漏驗證"  (chap4_list.md Sec. D, Stage 3)
%
%   Purpose  : verify the closed-form leakage pattern that appears when a
%              path's Doppler tap k_i does NOT land exactly on an integer
%              Doppler bin.
%   Theory   : for integer delay l, single path (g_i, l, k_i),
%                nu_tilde_{m,l}[n] = g_i * z^{k_i(m-l)} * e^{j*2*pi*k_i*n/N}   (4.73)
%              Its forward DFT over n (4.77) is a geometric sum whose closed
%              form is the periodic sinc / Dirichlet kernel (4.79)-(4.80):
%                nu_{m,l}[k] = (1/N) * g_i * z^{k_i(m-l)} * D_N(k_i - k)
%              (full re-derivation in periodic_sinc.m's header). When k_i is
%              an integer this collapses EXACTLY to a single nonzero bin at
%              k=k_i (D_N(0)=N cancels the 1/N); when k_i is fractional, the
%              energy leaks into every Doppler bin, enveloped by D_N(.).
%   Method   : build g^s[l,q] via gen_gs.m (Stage 2, works for any k_i,
%              integer or fractional, as long as l is integer), extract
%              nu_tilde_{m,l}[n] directly from it, take its DFT numerically,
%              and compare against the closed-form (4.80) — for BOTH an
%              integer tap (kappa=2) and a fractional tap (kappa=2.5), to
%              reproduce Fig. 4.16.
%   Criterion: relative error, numeric DFT vs. closed form, < 1e-10 for
%              both cases.
%
% Input  : p      - parameter struct from otfs_params()   (optional)
%          doPlot - true to draw the Fig. 4.16 style reproduction (default false)
% Output : results - struct with .name .err_rel .tol .pass
% =========================================================================

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'Stage1'));
addpath(fullfile(thisDir, '..', 'Stage2'));

if nargin < 1 || isempty(p),      p = otfs_params();  end
if nargin < 2 || isempty(doPlot), doPlot = true;     end

M = p.M;  N = p.N;
m0 = 0;  l0 = 0;                     % arbitrary (magnitude doesn't depend on m,l)

fprintf('\n=== Sim-14 : fractional Doppler leakage (4.79)/(4.80)  (M=%d, N=%d) ===\n', M, N);

cases = struct('name', {'integer (kappa=2)', 'fractional (kappa=2.5)'}, 'kappa', {2, 2.5});
errAll = zeros(1, numel(cases));
nuNumAll = cell(1, numel(cases));
nuAnaAll = cell(1, numel(cases));

for c = 1:numel(cases)
    kappa = cases(c).kappa;
    [g, l, k] = gen_channel_taps(p, 'l', l0, 'k', kappa, 'randomPhase', false, 'grid', 'off');

    gs = gen_gs(p, g, l, k);                     % Stage 2's function; l integer, k may be fractional
    nu_tilde = zeros(N, 1);
    for n = 0:N-1
        q = m0 + n*M;
        nu_tilde(n+1) = gs(l+1, q+1);
    end

    nu_num = zeros(N, 1);
    for kk = 0:N-1
        nu_num(kk+1) = (1/N) * sum(nu_tilde.' .* exp(-1j*2*pi*kk*(0:N-1)/N));
    end

    kOut = (0:N-1).';
    nu_ana = (1/N) * g * p.z^(k*(m0-l)) * periodic_sinc(k - kOut, N);

    e = relerr(nu_num, nu_ana);
    errAll(c) = e;
    nuNumAll{c} = nu_num;
    nuAnaAll{c} = nu_ana;

    fprintf('    [%s]  ||nu_num - nu_analytic|| / ||nu_analytic|| = %.3e\n', cases(c).name, e);
end

results.name    = 'Sim-14 fractional Doppler leakage (4.79/4.80)';
results.err_rel = max(errAll);
results.tol     = 1e-10;
results.pass    = results.err_rel < results.tol;

fprintf('    -> %s\n', ternary(results.pass, 'PASS', 'FAIL'));

if doPlot
    figure('Name', 'Sim-14 Fig. 4.16 reproduction', 'Color', 'w');
    kc = linspace(-0.5, N-0.5, 400);
    for c = 1:numel(cases)
        subplot(numel(cases), 1, c);
        stem(0:N-1, abs(nuNumAll{c}), 'filled'); hold on;
        env = (1/N) * abs(periodic_sinc(cases(c).kappa - kc, N));
        plot(kc, env, 'r--', 'LineWidth', 1.2);
        grid on; xlim([-0.5, N-0.5]);
        xlabel('Doppler bin k'); ylabel('|\nu_{m,l}[k]|');
        title(sprintf('Sim-14 : %s  (envelope = D_N(\\kappa-k)/N)', cases(c).name));
        legend('numeric (from G)', 'analytic envelope', 'Location', 'northeast');
    end
end

end

% -------------------------------------------------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
