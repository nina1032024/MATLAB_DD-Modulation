%% OTFS: convolution path vs. closed-form channel relation consistency check
%
%   Path A (ground truth):  X --IDZT--> s --G(off-grid)--> r --DZT--> Y_DD1
%   Path B (model):         X --closed-form H_DD / Table 4.3--> Y_DD2
%
%   Test 1  integer taps, four variants, matrix form + symbol-wise form -> should be machine precision
%   Test 2  fractional Doppler sweep                                    -> measures leakage
%   Test 3  TF domain: exact Hcheck identity vs. ideal single-tap product
%   Test 4  Doppler spread vector visualization (Fig. 4.16)
%
%   Requires otfs_build_G.m and otfs_dd_closed_form.m
clear; close all; clc;

%% ---------------- system parameters ----------------
M   = 16;          % delay bins  (number of subcarriers)
N   = 8;           % Doppler bins (number of time slots)
Lg  = 6;           % CP/ZP length (per block)
df  = 15e3;        % subcarrier spacing Delta f  (Hz)
T   = 1/df;  Ts = T/M;  Tf = N*T;
NM  = M*N;

fprintf('M=%d N=%d  Ts=%.2f us  T=%.1f us  Tf=%.2f ms\n', M,N,Ts*1e6,T*1e6,Tf*1e3);
fprintf('delay resolution %.2f us,  Doppler resolution %.1f Hz\n\n', Ts*1e6, 1/Tf);

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
fprintf('--- Test 1: integer delay/Doppler taps ---\n');
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

    fprintf('%4s : matrix %8.2e   symbol-wise %8.2e (%.1f dB)\n', ...
        v{1}, nmse(Y1,Y2mat), nmse(Y1,Y2sym), 10*log10(nmse(Y1,Y2sym)));
end
fprintf('  (RZP''s symbol-wise result is the (4.122) large-N approximation, so it will never be exactly 0)\n\n');

%% ================= Test 2: fractional Doppler sweep =================
fprintf('--- Test 2: fractional Doppler (CP-OTFS) ---\n');
p = base_params(M,N,Lg,'CP');
p.gains = gains;  p.ell = li;  p.lmax = max(li);
gam = 1 + Lg/M;

fr  = sort([0:0.05:0.5, 0.005]);
e_i = zeros(size(fr));  e_f = zeros(size(fr));
for a = 1:numel(fr)
    p.kappa = (ki + fr(a))/gam;            % fractional Doppler
    G  = otfs_build_G(p);
    s  = reshape(X*FN', [], 1);
    r  = G*s;
    Y1 = reshape(r, M, N) * FN;
    e_i(a) = nmse(Y1, otfs_dd_closed_form(X,p,'integer'));
    e_f(a) = nmse(Y1, otfs_dd_closed_form(X,p,'fractional'));
    fprintf('  frac=%.2f : integer %7.2f dB   fractional %8.2e\n', ...
        fr(a), 10*log10(e_i(a)), e_f(a));
end

figure; plot(fr, 10*log10(e_i), '-o', 'LineWidth',1.2); grid on; hold on;
plot(fr, 10*log10(max(e_f,1e-32)), '-s', 'LineWidth',1.2);
xlabel('fractional Doppler offset'); ylabel('NMSE (dB)');
legend('Integer-tap approximation (4.118)','Fractional formula (4.105)','Location','east');
title('Closed-Form vs. Convolution-Based Ground Truth');

%% ================= Test 3: time-frequency domain =================
fprintf('\n--- Test 3: TF domain ---\n');
p.kappa = (ki + 0.3)/gam;
G  = otfs_build_G(p);
s  = reshape(X*FN', [], 1);
r  = G*s;

% (a) exact identity (4.42): always holds
Hc  = kron(speye(N),FM) * G * kron(speye(N),FM');
xc  = kron(speye(N),FM) * s;
yc  = kron(speye(N),FM) * r;
fprintf('  exact Hcheck identity NMSE = %.2e\n', nmse(yc, Hc*xc));

% (b) ideal-pulse single-tap product (4.13)-(4.14): only an approximation
Xtf  = FM * X * FN';
Ytf1 = FM * reshape(r, M, N);
l = (0:M-1).';  k = 0:N-1;
Htf = zeros(M,N);
for i = 1:numel(p.gains)
    Htf = Htf + p.gains(i) * exp(1j*2*pi*p.kappa(i)*gam*k/N) .* exp(-1j*2*pi*l*p.ell(i)/M);
end
fprintf('  ideal single-tap product  NMSE = %.2f dB  <-- biorthogonality loss\n', ...
    10*log10(nmse(Ytf1, Htf.*Xtf)));

% (c) quantify off-diagonal energy
Hn0 = Hc(1:M,1:M);
fprintf('  Hcheck_{0,0} off-diagonal energy fraction = %.1f %%\n', ...
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

%% ---------------- local functions ----------------
function p = base_params(M,N,Lg,variant)
p = struct('M',M,'N',N,'Lg',Lg,'lmax',Lg,'Lsinc',8,'variant',variant, ...
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
