%% OTFS: convolution path vs. closed-form channel relation consistency check
%
%   Path A (ground truth):  X --IDZT--> s --G(off-grid)--> r --DZT--> Y_DD1
%   Path B (model):         X --closed-form H_DD / Table 4.3--> Y_DD2
%
%   Test 1     integer taps, four variants, matrix form + symbol-wise form -> should be machine precision
%   Test 2a-i  Doppler-only sweep (ell integer, KERNEL='none')            -> Dirichlet leakage, closed-form overlay
%   Test 2a-ii Delay-only sweep, 3x3 kernel identifiability heatmap       -> pulse shape must match model kernel
%   Test 2a-ii(b) guard-length sweep for the sinc kernel                  -> non-causal error floor, not fixed by L_g
%   Test 2a-iii phase convention: z^{kappa(m-ell_i)} vs z^{kappa(m-l)}    -> ~30 dB gap on fractional delay
%   Test 2d    LS channel estimation, bias-variance crossover (integer vs fractional model, MSE vs SNR)
%   Test 3     TF domain: exact Hcheck identity vs. ideal single-tap product
%   Test 4     Doppler spread vector visualization (Fig. 4.16)
%
%   Requires otfs_build_G.m, otfs_dd_closed_form.m, dirichlet_N.m, delay_kernel.m
clear; close all; clc;

% mirror everything printed to the terminal into output.txt as well, via
% pf() below rather than diary() -- diary's on-disk buffering does not
% reliably survive a mid-script crash, so we own the file handle directly
% and guarantee it gets closed (and thus flushed) with a try/catch.
out_file = fullfile(fileparts(mfilename('fullpath')), 'output.txt');
fid = fopen(out_file, 'w');
if fid == -1
    warning('otfs_verify:log','could not open %s for writing -- console-only output', out_file);
else
    fprintf('logging console output to: %s\n', out_file);
end

try

%% ---------------- system parameters ----------------
M   = 16;          % delay bins  (number of subcarriers)
N   = 8;           % Doppler bins (number of time slots)
Lg  = 6;           % CP/ZP length (per block)
df  = 15e3;        % subcarrier spacing Delta f  (Hz)
T   = 1/df;  Ts = T/M;  Tf = N*T;
NM  = M*N;

pf(fid, 'M=%d N=%d  Ts=%.2f us  T=%.1f us  Tf=%.2f ms\n', M,N,Ts*1e6,T*1e6,Tf*1e3);
pf(fid, 'delay resolution %.2f us,  Doppler resolution %.1f Hz\n\n', Ts*1e6, 1/Tf);

% normalized unitary DFT matrices (no Signal Processing Toolbox required)
FN = exp(-1j*2*pi*(0:N-1).'*(0:N-1)/N)/sqrt(N);
FM = exp(-1j*2*pi*(0:M-1).'*(0:M-1)/M)/sqrt(M);

% row-column interleaver P:  s = P*(I_M kron F_N')*x
pidx = reshape(reshape(1:NM, N, M).', [], 1);
P    = sparse(1:NM, pidx, 1, NM, NM);

% transmit symbols (QAM would also work)
rng(1);
X = (randn(M,N) + 1j*randn(M,N))/sqrt(2);

%% ================= Test 1: integer taps, four variants =================
pf(fid, '--- Test 1: integer delay/Doppler taps ---\n');
li = [0; 2; 5];            % integer delay taps
ki = [1; -2; 3];           % integer Doppler taps (for CP/ZP these are "post-scaling")
gains = [0.8; 0.5+0.2j; 0.3];

for v = {'CP','ZP','RCP','RZP'}
    p = base_params(M,N,Lg,v{1});
    gam      = strcmpi(v{1},'CP')*Lg/M + strcmpi(v{1},'ZP')*Lg/M + 1;
    p.gains  = gains;
    p.ell    = li;
    p.kappa  = ki/gam;     % makes kappa*gamma_g = ki an integer
    p.lmax   = max(li);

    % --- Path A ---
    G  = otfs_build_G(p);
    s  = reshape(X*FN', [], 1);            % IDZT + vec  (4.21)
    r  = G*s;                              % channel convolution
    Y1 = reshape(r, M, N) * FN;            % DZT         (4.29)

    % --- Path B-1: matrix identity (4.60) ---
    x = reshape(X.', [], 1);               % vec(X^T)
    H = kron(speye(M),FN) * P.' * G * P * kron(speye(M),FN');
    Y2mat = reshape(H*x, N, M).';

    % --- Path B-2: symbol-wise Table 4.3 ---
    Y2sym = otfs_dd_closed_form(X, p, 'integer');

    pf(fid, '%4s : matrix %8.2e   symbol-wise %8.2e (%.1f dB)\n', ...
        v{1}, nmse(Y1,Y2mat), nmse(Y1,Y2sym), 10*log10(nmse(Y1,Y2sym)));
end
pf(fid, '  (RZP''s symbol-wise result is the (4.122) large-N approximation, so it will never be exactly 0)\n\n');

%% ================= Test 2a-i: Doppler-only sweep =================
% ell is locked to integers -> KERNEL='none' means the delay kernel is a
% pure Kronecker delta, i.e. no interpolation assumption is exercised at
% all here. Only the Doppler-domain Dirichlet kernel D_N(.) is on trial.
pf(fid, '--- Test 2a-i: Doppler-only fractional sweep (ell integer, KERNEL=none) ---\n');
p = base_params(M,N,Lg,'CP');
p.gains  = gains;  p.ell = li;  p.lmax = max(li);
p.kernel = 'none';
gam = 1 + Lg/M;

fr     = sort([0:0.05:0.5, 0.005]);
e_i    = zeros(size(fr));  e_f = zeros(size(fr));  e_pred = zeros(size(fr));
for a = 1:numel(fr)
    p.kappa = (ki + fr(a))/gam;            % fractional Doppler
    G  = otfs_build_G(p);
    s  = reshape(X*FN', [], 1);
    r  = G*s;
    Y1 = reshape(r, M, N) * FN;
    e_i(a)    = nmse(Y1, otfs_dd_closed_form(X,p,'integer'));
    e_f(a)    = nmse(Y1, otfs_dd_closed_form(X,p,'fractional'));
    e_pred(a) = 2 - 2*real(dirichlet_N(fr(a),N))/N;    % closed-form leakage prediction (S3/S4 in two_path_modify.md)
    pf(fid, '  frac=%.3f : integer %7.2f dB   fractional %8.2e   closed-form %7.2f dB\n', ...
        fr(a), 10*log10(e_i(a)), e_f(a), 10*log10(e_pred(a)));
end

% --- single-path diagnostic: with only one unit-gain path, the integer
% curve should match the closed-form prediction to ~0.1 dB; the residual
% gap seen in the 3-path curve above is then attributable to cross-path
% terms (the closed form 2-2Re{D_N}/N is derived for a single path).
pf(fid, '  -- single-path check (g=1, ell=2, kappa_gam=-2+eps) --\n');
p1 = base_params(M,N,Lg,'CP');
p1.gains = 1;  p1.ell = 2;  p1.lmax = 2;  p1.kernel = 'none';
e_i_1p = zeros(size(fr));
for a = 1:numel(fr)
    p1.kappa = (-2 + fr(a))/gam;
    G  = otfs_build_G(p1);
    s  = reshape(X*FN', [], 1);
    r  = G*s;
    Y1 = reshape(r, M, N) * FN;
    e_i_1p(a) = nmse(Y1, otfs_dd_closed_form(X,p1,'integer'));
    pf(fid, '    frac=%.3f : integer %7.2f dB   closed-form %7.2f dB   diff %5.2f dB\n', ...
        fr(a), 10*log10(e_i_1p(a)), 10*log10(e_pred(a)), 10*log10(e_i_1p(a))-10*log10(e_pred(a)));
end

e_i_db    = 10*log10(e_i);
e_f_db    = 10*log10(max(e_f,1e-32));
e_pred_db = 10*log10(e_pred);  e_pred_db(~isfinite(e_pred_db)) = NaN;   % eps=0 -> log10(0) = -Inf
e_i1p_db  = 10*log10(e_i_1p);

figure;
subplot(1,2,1);
plot(fr, e_i_db, '-o', 'LineWidth',1.2); grid on; hold on;
plot(fr, e_f_db, '-s', 'LineWidth',1.2);
plot(fr, e_pred_db, 'k--', 'LineWidth',1.2);
plot(fr, e_i1p_db, ':d', 'LineWidth',1.0);
xlabel('fractional Doppler offset \epsilon'); ylabel('NMSE (dB)');
ylim([-320 20]); title('full range');

subplot(1,2,2);
plot(fr, e_i_db, '-o', 'LineWidth',1.2); grid on; hold on;
plot(fr, e_f_db, '-s', 'LineWidth',1.2);
plot(fr, e_pred_db, 'k--', 'LineWidth',1.2);
plot(fr, e_i1p_db, ':d', 'LineWidth',1.0);
xlabel('fractional Doppler offset \epsilon'); ylabel('NMSE (dB)');
ylim([-40 5]); title('zoom');
legend('Integer-tap (4.118), 3 paths','Fractional (4.105)', ...
       'Closed-form 2-2Re(D_N(\epsilon))/N','Integer-tap, 1 path','Location','southeast');
sgtitle('Result 2a-i - Doppler-only sweep (ell integer, sinc = delta, not exercised here)');

%% ================= Test 2a-ii: Delay-only, kernel identifiability (3x3) =================
% Doppler is locked to integers; only the delay kernel is on trial. Each
% row is a ground-truth Tx/Rx pulse, each column is the kernel the model
% assumes -- only the diagonal (matched kernel) should be exact.
pf(fid, '\n--- Test 2a-ii: Delay-only sweep, ground-truth pulse vs model kernel (3x3) ---\n');
kernels = {'sinc','tri','zoh'};
p2 = base_params(M,N,Lg,'CP');
p2.gains = gains;  p2.ell = li + 0.4;  p2.kappa = ki/gam;
p2.lmax  = max(li) + 1 + p2.Lsinc;     % must cover the widest (sinc) kernel's support

NMSE_33 = zeros(3,3);
for aTruth = 1:3
    p2.kernel = kernels{aTruth};       % ground-truth pulse
    G  = otfs_build_G(p2);
    s  = reshape(X*FN', [], 1);
    r  = G*s;
    Y1 = reshape(r, M, N) * FN;
    for bModel = 1:3
        pm = p2; pm.kernel = kernels{bModel};    % model kernel
        Y2 = otfs_dd_closed_form(X, pm, 'fractional');
        NMSE_33(aTruth,bModel) = 10*log10(nmse(Y1,Y2));
    end
    pf(fid, '  truth=%-5s : model sinc %7.1f dB   tri %7.1f dB   zoh %7.1f dB\n', ...
        kernels{aTruth}, NMSE_33(aTruth,1), NMSE_33(aTruth,2), NMSE_33(aTruth,3));
end

figure; imagesc(NMSE_33); colorbar; clim([-40 0]); colormap(parula);
set(gca,'XTick',1:3,'XTickLabel',kernels,'YTick',1:3,'YTickLabel',kernels);
xlabel('model kernel'); ylabel('ground-truth pulse');
title('Result F3 - delay-kernel identifiability, NMSE (dB), darker = better');
for aTruth = 1:3
    for bModel = 1:3
        if NMSE_33(aTruth,bModel) > -20, txtColor = 'k'; else, txtColor = 'w'; end
        text(bModel,aTruth,sprintf('%.1f',NMSE_33(aTruth,bModel)), ...
             'HorizontalAlignment','center','Color',txtColor,'FontWeight','bold');
    end
end

%% ========= Test 2a-ii(b): guard-length sweep for the sinc kernel (F4) =========
% sinc is non-causal (negative taps); CP's per-block circular wrap cannot
% represent that correctly, leaving an error floor that more guard length
% alone does not fix. Contrast with the 'tri' diagonal above (~ -290 dB).
pf(fid, '\n--- Guard-length sweep: sinc kernel error floor (F4, non-causality) ---\n');
Lg_list = [6 10 14 20 30];
NMSE_Lg = zeros(size(Lg_list));
for a = 1:numel(Lg_list)
    gamA = 1 + Lg_list(a)/M;
    pg = base_params(M,N,Lg_list(a),'CP');
    pg.gains  = gains;  pg.ell = li + 0.4;  pg.kappa = ki/gamA;
    pg.kernel = 'sinc'; pg.lmax = max(li) + 1 + pg.Lsinc;
    G  = otfs_build_G(pg);
    s  = reshape(X*FN', [], 1);
    r  = G*s;
    Y1 = reshape(r, M, N) * FN;
    Y2 = otfs_dd_closed_form(X, pg, 'fractional');
    NMSE_Lg(a) = 10*log10(nmse(Y1,Y2));
    pf(fid, '  Lg=%2d : %6.1f dB\n', Lg_list(a), NMSE_Lg(a));
end

figure; plot(Lg_list, NMSE_Lg, '-o', 'LineWidth',1.2); grid on; hold on;
ylim([-30 0]);
yline(mean(NMSE_Lg), 'k--', 'LineWidth',1);
xlabel('guard length L_g'); ylabel('NMSE (dB)');
title('Result F4 - sinc kernel error floor vs. guard length');
text(Lg_list(1), -20, sprintf('%.2f dB total variation over L_g = %d \\rightarrow %d', ...
     max(NMSE_Lg)-min(NMSE_Lg), Lg_list(1), Lg_list(end)), ...
     'FontSize',12,'EdgeColor','k','BackgroundColor','w','Margin',4);

%% ================= Test 2a-iii: phase convention =================
% z^{kappa(m-ell_i)} (true delay) vs the book's z^{kappa(m-l)} (tap index),
% kernel fixed to 'tri' so any gap here is purely the phase bug (F5).
pf(fid, '\n--- Test 2a-iii: phase convention z^(kappa*(m-ell_i)) vs z^(kappa*(m-l)) ---\n');
p3 = base_params(M,N,Lg,'CP');
p3.gains  = gains;
p3.ell    = li + 0.4;
p3.kappa  = (ki + 0.3)/gam;
p3.kernel = 'tri';
p3.lmax   = max(li) + 1;

G  = otfs_build_G(p3);
s  = reshape(X*FN', [], 1);
r  = G*s;
Y1 = reshape(r, M, N) * FN;

NMSE_true = 10*log10(nmse(Y1, otfs_dd_closed_form(X, p3, 'fractional')));  % corrected: uses ell_i
NMSE_book = 10*log10(nmse(Y1, otfs_dd_bookphase(X, p3)));                  % book convention: uses l

pf(fid, '  correct phase z^(kappa*(m-ell_i)) : %6.1f dB\n', NMSE_true);
pf(fid, '  book   phase z^(kappa*(m-l))      : %6.1f dB\n', NMSE_book);

figure; barh([NMSE_true, NMSE_book]);
xlim([-320 0]); xlabel('NMSE (dB)'); grid on;
set(gca,'YTickLabel',{'$z^{\kappa(m-\ell_i)}$ (correct)','$z^{\kappa(m-l)}$ (book)'}, ...
    'TickLabelInterpreter','latex');
title('Result F5: phase convention -- true delay $\ell_i$ vs. tap index $l$', ...
    'Interpreter','latex');

%% ================= Test 2d: bias-variance crossover (LS channel estimation) =================
% Not a ground-truth-consistency test like 2a-i/ii/iii -- a channel
% ESTIMATION question. Per pilot-delay bin m, the received Y[m,:] depends
% only on nu_{m,ell_i}[k] for that same m (paths' delays ell_i assumed
% known/support-known), so LS decouples completely across m: solve M
% independent (T*N) x (P*N) systems instead of one big one.
%   Model A (integer)    : estimate only k = round(kappa_i*gamma_g)  -> P unknowns per m
%   Model B (fractional) : estimate all N bins                       -> P*N unknowns per m
% See two_path_modify.md S5 "Test 2d" for the full spec.
pf(fid, '\n--- Test 2d: bias-variance crossover, LS channel estimation ---\n');

T      = 8;                    % pilot frames, channel block-static
TRIALS = 200;
Ppath  = numel(li);            % P = 3 known path delays li = [0;2;5]
z2d    = exp(1j*2*pi/NM);

rng(2);
Xs = exp(1j*pi/2*randi([0 3], M, N, T));     % QPSK pilots, T frames

idx = zeros(Ppath*N,2);                      % column layout: (i=1,k=0..N-1),(i=2,...),...
ci = 0;
for i = 1:Ppath
    for k = 0:N-1
        ci = ci+1; idx(ci,:) = [i k];
    end
end

% --- pilot dictionary + pinv for Model B: independent of eps/SNR/trial, hoist out ---
AB    = cell(M,1);
pinvB = cell(M,1);
for m = 0:M-1
    Am = zeros(T*N, Ppath*N);
    for cc = 1:size(idx,1)
        i = idx(cc,1); k = idx(cc,2);
        r = 0;
        col = zeros(T*N,1);
        for t = 1:T
            for n = 0:N-1
                r = r+1;
                col(r) = Xs(mod(m-li(i),M)+1, mod(n-k,N)+1, t);
            end
        end
        Am(:,cc) = col;
    end
    AB{m+1}    = Am;
    pinvB{m+1} = pinv(Am);
end

eps_list = [0.005 0.05 0.3];
snrs     = -10:5:40;

figure;
tl = tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');
for eIdx = 1:numel(eps_list)
    epsv = eps_list(eIdx);
    kg   = ki + epsv;                 % kappa_i * gamma_g per path
    kap  = kg/gam;                    % kappa_i (unscaled)
    kA   = mod(round(kg), N);         % Model A's single dominant bin per path

    % --- true Doppler spread vector nu_{m,i,k}, (4.105)/S3 with corrected phase ---
    nu_t = zeros(M, Ppath, N);
    for i = 1:Ppath
        Dk = dirichlet_N(kg(i) - (0:N-1), N);          % 1 x N
        for m = 0:M-1
            nu_t(m+1,i,:) = gains(i) * z2d^(kap(i)*(m-li(i))) * Dk / N;
        end
    end

    % true-signal energy for NMSE normalization (independent of SNR/trial)
    th_all = cell(M,1);  sigK = 0;
    for m = 0:M-1
        th_all{m+1} = reshape(squeeze(nu_t(m+1,:,:)).', [], 1);   % (P*N) x 1, order matches idx
        sigK = sigK + norm(th_all{m+1})^2;
    end

    % Model A's selected columns/pinv depend on eps (via kA) -> hoist out of SNR/trial loops
    sel = false(size(idx,1),1);
    for i = 1:Ppath
        sel = sel | (idx(:,1)==i & idx(:,2)==kA(i));
    end
    pinvA = cell(M,1);
    for m = 0:M-1
        pinvA{m+1} = pinv(AB{m+1}(:,sel));
    end

    % LS bias floor for Model A -- NOT the 2a-i formula (that one includes
    % the retained bin's amplitude error too; LS re-fits the amplitude, so
    % only the un-recoverable leaked energy remains). See S5 note.
    bias_ls = 1 - abs(dirichlet_N(epsv, N))^2 / N^2;

    MSE_A = zeros(size(snrs));  MSE_B = zeros(size(snrs));
    for sIdx = 1:numel(snrs)
        sig2 = 10^(-snrs(sIdx)/10);        % unit pilot power per symbol
        eA = 0; eB = 0;
        for tr = 1:TRIALS
            for m = 0:M-1
                th    = th_all{m+1};
                noise = sqrt(sig2/2)*(randn(T*N,1)+1j*randn(T*N,1));
                y     = AB{m+1}*th + noise;
                eB    = eB + norm(pinvB{m+1}*y - th)^2;
                eA    = eA + norm(pinvA{m+1}*y - th(sel))^2 + (norm(th)^2 - norm(th(sel))^2);
            end
        end
        MSE_A(sIdx) = 10*log10(eA/(TRIALS*sigK));
        MSE_B(sIdx) = 10*log10(eB/(TRIALS*sigK));
    end

    % crossover SNR: linearly interpolate on the (MSE_A-MSE_B) sign change
    % rather than snapping to the nearest sampled SNR grid point (5 dB steps
    % would otherwise silently round every crossover to a multiple of 5).
    d  = MSE_A - MSE_B;
    jc = find(d(1:end-1) < 0 & d(2:end) > 0, 1, 'first');
    if isempty(jc)
        cross_snr = NaN;
    else
        cross_snr = snrs(jc) + (snrs(jc+1)-snrs(jc)) * (-d(jc))/(d(jc+1)-d(jc));
    end
    pf(fid, '  eps=%.3f : bias floor %6.2f dB, crossover SNR ~ %.1f dB\n', epsv, 10*log10(bias_ls), cross_snr);

    nexttile;
    if ~isnan(cross_snr)
        patch([-10 cross_snr cross_snr -10], [-50 -50 10 10], [0.9 0.95 0.9], 'EdgeColor','none');
        uistack(findobj(gca,'Type','patch'), 'bottom');
    end
    hold on;
    hA = plot(snrs, MSE_A, '-o');
    hB = plot(snrs, MSE_B, '-s');
    yline(10*log10(bias_ls), 'k--');
    if ~isnan(cross_snr)
        xline(cross_snr, ':', sprintf('%.1f dB', cross_snr));
    end
    ylim([-50 10]); xlim([-10 50]); grid on;
    xlabel('pilot SNR (dB)'); ylabel('channel MSE (dB)');
    title(sprintf('$\\varepsilon = %.3f$', epsv), 'Interpreter','latex');
    if eIdx == 1
        % low-SNR gap between the two curves is set by the parameter-count
        % ratio alone (variance ~ #params), independent of eps: ~10log10(N)
        iArrow = 2;
        xa = snrs(iArrow); yA = MSE_A(iArrow); yB = MSE_B(iArrow);
        plot([xa xa], [yA yB], 'k-', 'LineWidth',1);
        plot(xa, yA, 'k^', 'MarkerFaceColor','k', 'MarkerSize',5);
        plot(xa, yB, 'kv', 'MarkerFaceColor','k', 'MarkerSize',5);
        text(xa+2, (yA+yB)/2, sprintf('%.1f dB \\approx 10log_{10}(N)', yB-yA), 'FontSize',9);
    end
end

% one shared legend in its own tile below the row, so it can never
% overlap plot content (unlike an in-axes 'northeast'/'Position' legend)
lgd = legend([hA hB], {'Integer model (P params)','Fractional model (P{\cdot}N params)'}, ...
    'Orientation','horizontal');
lgd.Layout.Tile = 'south';
title(tl, 'Result F6 - bias-variance crossover: channel estimation MSE vs. pilot SNR');

%% ================= Test 3: time-frequency domain =================
pf(fid, '\n--- Test 3: TF domain ---\n');
p = base_params(M,N,Lg,'CP');
p.gains  = gains;  p.ell = li;  p.lmax = max(li);
p.kernel = 'tri';
p.kappa = (ki + 0.3)/gam;
G  = otfs_build_G(p);
s  = reshape(X*FN', [], 1);
r  = G*s;

% (a) exact identity (4.42): always holds
Hc  = kron(speye(N),FM) * G * kron(speye(N),FM');
xc  = kron(speye(N),FM) * s;
yc  = kron(speye(N),FM) * r;
pf(fid, '  exact Hcheck identity NMSE = %.2e\n', nmse(yc, Hc*xc));

% (b) ideal-pulse single-tap product (4.13)-(4.14): only an approximation
Xtf  = FM * X * FN';
Ytf1 = FM * reshape(r, M, N);
l = (0:M-1).';  k = 0:N-1;
Htf = zeros(M,N);
for i = 1:numel(p.gains)
    Htf = Htf + p.gains(i) * exp(1j*2*pi*p.kappa(i)*gam*k/N) .* exp(-1j*2*pi*l*p.ell(i)/M);
end
pf(fid, '  ideal single-tap product  NMSE = %.2f dB  <-- biorthogonality loss\n', ...
    10*log10(nmse(Ytf1, Htf.*Xtf)));

% (c) quantify off-diagonal energy
Hn0 = Hc(1:M,1:M);
pf(fid, '  Hcheck_{0,0} off-diagonal energy fraction = %.1f %%\n', ...
    100*(norm(Hn0,'fro')^2 - norm(diag(Hn0))^2)/norm(Hn0,'fro')^2);

%% ================= Test 4: Doppler spread vector =================
m0 = 6; l0 = 2;
figure;
for a = 1:2
    fracv = (a-1)*0.5;
    p.kappa = (ki + fracv)/gam;
    nu = zeros(1,N);
    for n = 0:N-1
        q = m0 + n*(M+Lg);
        nu(n+1) = sum(p.gains .* exp(1j*2*pi/NM).^(p.kappa*(q-l0)) .* ...
                      sinc_local(l0-p.ell, p.Lsinc));
    end
    nu_dd = fft(nu)/N;                      % (4.102)
    subplot(2,1,a); stem(0:N-1, abs(nu_dd), 'filled'); grid on;
    xlabel('Doppler bin k'); ylabel('|\nu_{m,l}[k]|');
    title(sprintf('m=%d, l=%d, fractional offset = %.1f', m0, l0, fracv));
end

if fid > 0, fclose(fid); end

catch ME
    if fid > 0, fclose(fid); end
    rethrow(ME);
end

%% ---------------- local functions ----------------
function pf(fid, varargin)
%PF  print to the console and, if open, mirror the same text into fid.
fprintf(varargin{:});
if fid > 0
    fprintf(fid, varargin{:});
end
end

function p = base_params(M,N,Lg,variant)
p = struct('M',M,'N',N,'Lg',Lg,'lmax',Lg,'Lsinc',8,'variant',variant,'kernel','tri', ...
           'gains',[],'ell',[],'kappa',[]);
end

function e = nmse(a,b)
e = sum(abs(a(:)-b(:)).^2) / sum(abs(a(:)).^2);
end

function y = sinc_local(x,L)
x = x(:);
y = sin(pi*x)./(pi*x);
y(x==0) = 1;
y(abs(x-round(x))<1e-12 & x~=0) = 0;
y(abs(x)>L) = 0;
end

function Y = otfs_dd_bookphase(X, p)
%OTFS_DD_BOOKPHASE  Same as otfs_dd_closed_form(X,p,'fractional') but with
%   the book's (4.105)/(4.118) phase convention z^{kappa_i(m-l)} (tap
%   index l) instead of the corrected z^{kappa_i(m-ell_i)} (true delay).
%   Kept only in the verification script, for the Test 2a-iii comparison
%   -- see F4/F5 in two_path_modify.md.
if ~any(strcmpi(p.variant,{'CP','ZP'}))
    error('otfs_dd_bookphase currently supports CP / ZP only');
end
M = p.M;  N = p.N;  NM = M*N;
z = exp(1j*2*pi/NM);
switch upper(p.variant)
    case {'CP','ZP'}, gam = 1 + p.Lg/M;
    otherwise,        gam = 1;
end
g   = p.gains(:);
ell = p.ell(:);
kap = p.kappa(:);
if isfield(p,'kernel') && ~isempty(p.kernel), kernel = p.kernel; else, kernel = 'tri'; end
Y = zeros(M,N);
k = 0:N-1;
for m = 0:M-1
    for l = 0:p.lmax
        if strcmpi(p.variant,'ZP') && m < l, continue; end
        sw = delay_kernel(l-ell, kernel);
        if strcmpi(kernel,'sinc'), sw(abs(l-ell) > p.Lsinc) = 0; end
        if all(sw==0), continue; end
        nu = (1/N) * ( ( g .* z.^(kap*(m-l)) .* sw ).' * dirichlet_N(kap*gam - k, N) );
        mm = mod(m-l,M);
        for n = 0:N-1
            Y(m+1,n+1) = Y(m+1,n+1) + nu * X(mm+1, mod(n-k,N)+1).';
        end
    end
end
end
