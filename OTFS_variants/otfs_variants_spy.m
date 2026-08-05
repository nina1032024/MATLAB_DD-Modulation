%% otfs_variants_spy.m  (rev. 2)
%  Matrix structure of the four OTFS variants (CP / ZP / RCP / RZP) in three domains
%    rows = domain: time G, delay-time H~, delay-Doppler H
%    columns = variant
%
%  Corresponds to book Fig. 4.8/4.11/4.12/4.14/4.15/4.17/4.18/4.19 and Table 4.2/4.4
%
%  ---- rev.2 fixes relative to rev.1 ----
%  A1 (important): no longer forces gamma_g = 1. CP/ZP sample g^s at
%    q = m + n(M+Lg) per (4.93)/(4.109); RCP/RZP sample at q = m + nM per
%    (4.38)/(4.83). So the "matrix row index" and the "physical sampling
%    instant" are now explicitly separated.
%    Consequence: (4.107) requires kappa*gamma_g to be an integer, not kappa itself.
%    With Lg=2, M=8 => gamma_g=1.25, so a physically-integer Doppler is still
%    off-grid for CP/ZP.
%    Set USE_GAMMA=false to fall back to rev.1's gamma_g=1 (for comparison only,
%    does not represent a real system).
%
%  ---- modeling assumptions that still hold (should be noted in any slides) ----
%  A2: ell and kap are paired one-to-one, so each delay tap has exactly one
%      Doppler value (|K_l|=1). In general
%      K_{m,l} = sum_{kap in K_l} nu_l(kap) z^{kap(m-l)} * Pi^{[kap]_N},
%      i.e. |K_l| circular diagonals; the 1-1/N collapse shown here is the
%      best case. Set DEMO_A2=true to run a |K_1|=2 counter-example.
%  A3: sinc is truncated at |l-ell|<=Ls. sinc(l-ell) is nonzero for every
%      integer l, so l_max is a modeling choice, not a channel property;
%      report the span together with Ls.
%  A4: spy only shows the support set, not magnitude; a -40 dB sinc tail is
%      drawn just as solid black as a strong tap. Set SHOW_DB=true to switch
%      to a dB plot, which reveals the off-grid "all black" region as a thick
%      diagonal band.
%
%  ---- verified nnz (M=8, N=6, g=[1 .75 .55]) ----
%    on-grid, Lg=0 (gamma_g=1)   CP 144/144/144   ZP 126/126/126
%                                RCP 144/144/144  RZP 141/141/234
%    on-grid, Lg=2 (gamma_g=1.25) CP 144/144/624  ZP 126/126/516
%                                RCP 144/144/144  RZP 141/141/234   <-- ranking flips
%    off-grid, Ls=3               CP 288/288/1728 ZP 198/198/1188
%                                RCP 288/288/1728 RZP 273/273/1698

clear; close all; clc;

%% ---------------- switches ----------------
USE_GAMMA = true;    % true: CP/ZP use the gamma_g from (4.93)/(4.109); false: gamma_g=1 everywhere
SHOW_DB   = false;   % true: plot dB magnitude instead of the binary support set
DEMO_A2   = false;   % true: additionally run the |K_1|=2 counter-example
DB_FLOOR  = -40;     % dynamic range floor for the dB plot

%% ---------------- parameters ----------------
M = 8;  N = 6;  L = M*N;
LM = 7;                      % max l in the delay tap table
g  = [1.00, 0.75, 0.55];     % gains of the three paths

% cfg.Lg: guard length per block, must be >= the effective l_max
% on-grid (M1)
cfg(1).ell = [0 1 2];        cfg(1).kap = [0 1 2];
cfg(1).Ls  = 0.5;            cfg(1).Lg  = 2;
cfg(1).tag = 'integer delay and Doppler  (on-grid, M1)';
cfg(1).fn  = 'otfs_variants_spy_ongrid.png';

% off-grid (M3): with Ls=3 the effective span of ell=[0.3 1.7 2.4] is l = 0..5, hence Lg=5
cfg(2).ell = [0.3 1.7 2.4];  cfg(2).kap = [0.3 1.7 -1.2];
cfg(2).Ls  = 3;              cfg(2).Lg  = 5;
cfg(2).tag = 'fractional delay and Doppler  (off-grid, M3)';
cfg(2).fn  = 'otfs_variants_spy_offgrid.png';

if DEMO_A2                   % A2 counter-example: delay tap l=1 carries two Doppler values
  cfg(3).ell = [1 1 2];      cfg(3).kap = [0 2 1];
  cfg(3).Ls  = 0.5;          cfg(3).Lg  = 2;
  cfg(3).tag = 'on-grid with |K_1| = 2  (A2 counter-example)';
  cfg(3).fn  = 'otfs_variants_spy_A2.png';
end

variants = {'CP','ZP','RCP','RZP'};
rowName  = {'Time domain  G','Delay-time  $\tilde{H}$','Delay-Doppler  H'};
blkSize  = [M, N, N];

%% ---------------- shared matrices ----------------
% permutation matrix P: s = P * x_tilde, i.e. x_tilde(m*N+n) = s(n*M+m), (4.33)
P = zeros(L,L);
for m = 0:M-1
  for n = 0:N-1
    P(n*M+m+1, m*N+n+1) = 1;
  end
end
FN = fft(eye(N))/sqrt(N);        % toolbox-free replacement for dftmtx
A  = kron(eye(M), FN);           % the I_M kron F_N from (4.60)

%% ---------------- main loop ----------------
for c = 1:numel(cfg)
  Lg  = cfg(c).Lg * USE_GAMMA;             % falls back to Lg=0 when USE_GAMMA=false
  Q   = (M + Lg)*N;                        % need values up to the longest guard-inclusive time
  gs  = gs_table(g, cfg(c).ell, cfg(c).kap, cfg(c).Ls, LM, Q, L);

  figure('Position',[80 80 1120 890],'Color','w');
  fprintf('\n=== %s ===\n', cfg(c).tag);
  fprintf('    Lg = %d,  gamma_g(CP/ZP) = %.4g,  gamma_g(RCP/RZP) = 1\n', Lg, 1+Lg/M);
  fprintf('    kappa * gamma_g (CP/ZP) = [%s]\n', num2str(cfg(c).kap*(1+Lg/M),'%.3g  '));
  fprintf('%-6s %10s %12s %10s\n','variant','nnz(G)','nnz(Ht)','nnz(H)');

  for j = 1:4
    v  = variants{j};
    G  = build_G(v, gs, M, N, L, Lg);
    Ht = P.' * G * P;              % (4.55)
    H  = A * Ht * A';              % (4.60)
    mats = {G, Ht, H};  nz = zeros(1,3);

    for i = 1:3
      X = mats{i};  blk = blkSize(i);
      mask  = abs(X) > 1e-9*max(abs(X(:)));
      nz(i) = nnz(mask);

      ax = subplot(3,4,(i-1)*4 + j);
      if SHOW_DB
        Xd = 20*log10(abs(X)/max(abs(X(:))) + eps);
        imagesc(Xd); caxis([DB_FLOOR 0]); colormap(ax, flipud(gray));
      else
        imagesc(~mask); caxis([0 1]); colormap(ax, gray);
      end
      axis image; hold on;
      for t = 0:blk:L
        plot([0.5 L+0.5], [t+0.5 t+0.5], 'Color',[.22 .54 .87], 'LineWidth',0.6);
        plot([t+0.5 t+0.5], [0.5 L+0.5], 'Color',[.22 .54 .87], 'LineWidth',0.6);
      end
      set(ax,'XTick',[],'YTick',[],'Box','on');
      xlabel(sprintf('nnz = %d', nz(i)), 'FontSize',9);
      if i == 1, title([v '-OTFS'], 'FontSize',13, 'FontWeight','normal'); end
      if j == 1, ylabel(rowName{i}, 'Interpreter','latex', 'FontSize',12); end
    end
    fprintf('%-6s %10d %12d %10d\n', v, nz(1), nz(2), nz(3));
  end

  sgtitle(sprintf('M=%d, N=%d, L_g=%d (\\gamma_g=%.4g for CP/ZP), %s', ...
                  M, N, Lg, 1+Lg/M, cfg(c).tag), 'FontSize',12);
  exportgraphics(gcf, cfg(c).fn, 'Resolution',200);   % on releases before R2020a, use print instead
end

%% ---------------- P.10: which cells the on-grid assumption discards ----------------
% Baseline = true channel (fractional kappa); comparison = kappa rounded to an integer.
% ell is fixed at [0 1 2], only Doppler differs, so G / Ht have identical support sets.
plot_lost(g, M, N, L, LM, P, A, USE_GAMMA*2);

%% ================= subfunctions =================
function gs = gs_table(g, ell, kap, Ls, LM, Q, NM)
% discrete time-varying tap g^s[l,q], book Eq. (4.6)
%   gs(l+1,q+1) = sum_i g_i * sinc(l-ell_i) * z^{kap_i (q-l)},  z = exp(j2pi/NM)
% NM is always M*N (the Doppler normalization base; it does not change with the guard)
% Q is the number of samples to evaluate; variants with a per-block guard need Q = (M+Lg)*N > NM
% Ls is the sinc truncation half-width; on-grid, taking 0.5 degenerates to (4.8)
  gs = zeros(LM+1, Q);
  q  = 0:Q-1;
  for i = 1:numel(g)
    for l = 0:LM
      if abs(l - ell(i)) > Ls, continue; end
      w = sinc_n(l - ell(i));
      if abs(w) < 1e-13, continue; end
      gs(l+1,:) = gs(l+1,:) + g(i)*w*exp(1j*2*pi*kap(i)*(q-l)/NM);
    end
  end
end

function G = build_G(v, gs, M, N, L, Lg)
% assemble the time-domain channel matrix G per variant, (4.38)/(4.83)/(4.93)/(4.109)
%   q  : matrix row index with the guard removed, always m + n*M
%   qc : physical sample time. CP/ZP insert the guard between blocks, so qc = m + n*(M+Lg)
%        RCP/RZP have a single guard at the very front of the frame, which
%        does not affect inter-block spacing, so qc = q
  G = zeros(L,L);
  active = find(any(abs(gs) > 1e-13, 2)).' - 1;   % delay taps l (0-based) that are nonzero
  perBlockGuard = any(strcmp(v, {'CP','ZP'}));
  for n = 0:N-1
    for m = 0:M-1
      q  = m + n*M;
      if perBlockGuard, qc = m + n*(M+Lg); else, qc = q; end
      for l = active
        val = gs(l+1, qc+1);
        switch v
          case 'CP'                               % wraps within the block
            cc = n*M + mod(m-l,M) + 1;
            G(q+1, cc) = G(q+1, cc) + val;
          case 'ZP'                               % zeroed within the block
            if m >= l
              cc = n*M + (m-l) + 1;
              G(q+1, cc) = G(q+1, cc) + val;
            end
          case 'RCP'                              % wraps at the frame level
            cc = mod(q-l, L) + 1;
            G(q+1, cc) = G(q+1, cc) + val;
          case 'RZP'                              % zeroed at the frame level
            if q - l >= 0
              G(q+1, q-l+1) = G(q+1, q-l+1) + val;
            end
        end
      end
    end
  end
end

function plot_lost(g, M, N, L, LM, P, A, Lg)
% three-color plot: dark blue = nonzero in both; orange = nonzero in the true
% channel but zero in the on-grid approximation; white = zero in both.
% also light gray = nonzero on-grid but zero in the true channel (can happen
% for RZP due to cancellation of the all-ones term)
  variants = {'CP','ZP','RCP','RZP'};
  rowName  = {'Time domain  G','Delay-time  $\tilde{H}$','Delay-Doppler  H'};
  blkSize  = [M, N, N];
  Q  = (M+Lg)*N;
  gsA = gs_table(g, [0 1 2], [0.3 1.7 -1.2], 0.5, LM, Q, L);   % true channel (fractional kappa)
  gsB = gs_table(g, [0 1 2], [0   1   2   ], 0.5, LM, Q, L);   % on-grid approximation

  cmap = [1 1 1; 0.85 0.85 0.85; 0.90 0.55 0.20; 0.09 0.20 0.33];  % 0/1/2/3
  figure('Position',[80 80 1120 890],'Color','w');
  fprintf('\n=== which entries are destroyed by rounding kappa to integers (Lg=%d) ===\n', Lg);
  for j = 1:4
    v = variants{j};
    MA = cell(1,3); MB = cell(1,3);
    gsPair = {gsA, gsB};                 % MATLAB doesn't allow literal indexing like {a,b}{t}
    for t = 1:2
      gs = gsPair{t};
      G  = build_G(v, gs, M, N, L, Lg);
      Ht = P.' * G * P;  H = A * Ht * A';
      mats = {G, Ht, H};
      for i = 1:3
        X = mats{i};  msk = abs(X) > 1e-9*max(abs(X(:)));
        if t==1, MA{i} = msk; else, MB{i} = msk; end
      end
    end
    for i = 1:3
      base = MA{i};  on = MB{i};
      code = 3*(base & on) + 2*(base & ~on) + 1*(~base & on);
      surv = nnz(base & on);  lost = nnz(base & ~on);  tot = nnz(base);
      ax = subplot(3,4,(i-1)*4 + j);
      imagesc(code); caxis([0 3]); colormap(ax, cmap); axis image; hold on;
      for t2 = 0:blkSize(i):L
        plot([0.5 L+0.5],[t2+0.5 t2+0.5],'Color',[.7 .8 .9],'LineWidth',0.5);
        plot([t2+0.5 t2+0.5],[0.5 L+0.5],'Color',[.7 .8 .9],'LineWidth',0.5);
      end
      set(ax,'XTick',[],'YTick',[],'Box','on');
      xlabel(sprintf('base %d | surv %d | lost %d (%.0f%%)', ...
                     tot, surv, lost, 100*lost/max(tot,1)), 'FontSize',8);
      if i==1, title([v '-OTFS'],'FontSize',13,'FontWeight','normal'); end
      if j==1, ylabel(rowName{i},'Interpreter','latex','FontSize',12); end
      if i==3
        fprintf('  %-4s base=%4d  surv=%4d  lost=%4d (%.1f%%)\n', v, tot, surv, lost, 100*lost/tot);
      end
    end
  end
  sgtitle(sprintf(['Baseline = true channel (\\kappa=[0.3,1.7,-1.2]);  ' ...
                   'comparison = \\kappa rounded to [0,1,2].  M=%d, N=%d, ' ...
                   '\\ell=[0,1,2] in both'], M, N), 'FontSize',11);
  exportgraphics(gcf, 'otfs_variants_lost.png', 'Resolution',200);
end

function y = sinc_n(x)
% normalized sinc: sin(pi x)/(pi x), equal to 1 at x=0 (avoids needing the Signal Processing Toolbox)
  y = ones(size(x));
  k = (x ~= 0);
  y(k) = sin(pi*x(k)) ./ (pi*x(k));
end
