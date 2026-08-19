function G = otfs_build_G(p)
%OTFS_BUILD_G  Time-domain channel matrix G such that r = G*s + w.
%
%   This is the "ground truth": it directly discretizes the off-grid
%   (fractional delay / Doppler) channel into a convolution matrix. Every
%   other domain's relation must be derivable from this G.
%
%   p.M, p.N      : delay bins / Doppler bins
%   p.Lg          : guard length (CP/ZP length per block; unused by RCP/RZP)
%   p.lmax        : maximum delay tap (must be <= p.Lg for CP/ZP)
%   p.Lsinc       : sinc truncation half-width (only matters for KERNEL='sinc')
%   p.kernel      : 'none' | 'zoh' | 'tri' | 'sinc'  (default 'tri', matching
%                   the rectangular Tx/Rx pulse shaping the system declares;
%                   'sinc' is the old ideal band-limited default; 'none' is
%                   for integer-delay-only tests with zero interpolation)
%   p.variant     : 'CP' | 'ZP' | 'RCP' | 'RZP'
%   p.gains       : P x 1  complex gains g_i
%   p.ell         : P x 1  normalized delay  ell_i = tau_i/Ts        (may be fractional)
%   p.kappa       : P x 1  normalized Doppler kappa_i = nu_i*N*T     (may be fractional)
%
%   Note: the CP/ZP variants return the NM x NM equivalent matrix "with the
%   guard removed", consistent with the book's (4.93)/(4.109).
%
%   NOTE - phase convention: the Doppler tilt phase is z^{kappa_i (q - ell_i)},
%   using the *true* (possibly fractional) delay ell_i, not the discrete tap
%   index l. Physically the phase comes straight from e^{j2*pi*nu_i*(qTs -
%   tau_i)} (see two_path_modify.md S1) -- it is a per-path, q-dependent
%   rotation with a constant offset set by tau_i, and has nothing to do with
%   which tap l the delay kernel happens to be evaluated at. Substituting l
%   for ell_i here (as the book's (4.6) does) is only exact when ell_i is an
%   integer (l == ell_i); for fractional delay it is a ~30 dB error and, if
%   left in the "ground truth" G, silently makes G itself wrong.

M = p.M;  N = p.N;  NM = M*N;  Lg = p.Lg;  lmax = p.lmax;
z  = exp(1j*2*pi/NM);
g   = p.gains(:);
ell = p.ell(:);
kap = p.kappa(:);

if isfield(p,'kernel') && ~isempty(p.kernel)
    kernel = p.kernel;
else
    kernel = 'tri';
end

% candidate delay taps that can carry nonzero energy for this KERNEL
taps = kernel_taps(ell, kernel, p.Lsinc);
taps = taps(taps <= lmax);

% g^s[l,q] = sum_i g_i * z^{kappa_i (q-ell_i)} * delay_kernel(l - ell_i, KERNEL)   (4.6, corrected)
gs = @(l,q) sum( g .* z.^(kap.*(q-ell)) .* eval_kernel(l-ell, kernel, p.Lsinc) );

G = zeros(NM,NM);

switch upper(p.variant)

    case 'RZP'   % one ZP per frame -> pure linear convolution, lower triangular
        for q = 0:NM-1
            for l = taps
                if l < 0 || l > q, continue; end
                G(q+1, q-l+1) = gs(l,q);
            end
        end

    case 'RCP'   % one CP per frame -> circular over the whole frame
        for q = 0:NM-1
            for l = taps
                c = mod(q-l,NM)+1;
                G(q+1,c) = G(q+1,c) + gs(l,q);
            end
        end

    case 'CP'    % one CP per block -> each M x M diagonal block is circular on its own
        for n = 0:N-1
            for m = 0:M-1
                q = m + n*(M+Lg);            % true time index including guard
                for l = taps
                    c = n*M + mod(m-l,M) + 1;
                    G(n*M+m+1,c) = G(n*M+m+1,c) + gs(l,q);
                end
            end
        end

    case 'ZP'    % one ZP per block -> each M x M diagonal block is lower triangular
        for n = 0:N-1
            for m = 0:M-1
                q = m + n*(M+Lg);
                for l = taps
                    if l < 0 || l > m, continue; end
                    G(n*M+m+1, n*M+(m-l)+1) = gs(l,q);
                end
            end
        end

    otherwise
        error('otfs_build_G:variant','variant must be CP / ZP / RCP / RZP');
end
end

% ---------------------------------------------------------------
function y = eval_kernel(x, kernel, Lsinc)
%EVAL_KERNEL  delay_kernel(x, kernel), with the 'sinc' case truncated to
%             |x| <= Lsinc (the other kernels already have compact support).
y = delay_kernel(x, kernel);
if strcmpi(kernel,'sinc')
    y(abs(x) > Lsinc) = 0;
end
end

% ---------------------------------------------------------------
function ls = kernel_taps(ell, kernel, Lsinc)
%KERNEL_TAPS  Integer delay-tap indices l that can carry nonzero
%             delay_kernel(l-ell_i, KERNEL) energy for at least one path.
%             Deliberately over-covers (a couple of taps wider than the
%             kernel's exact support) and lets the caller's kernel
%             evaluation zero out anything outside the true support --
%             that is more robust than hand-deriving floor/ceil offsets
%             per kernel (e.g. 'zoh' needs ceil(ell_i), not floor(ell_i),
%             and it is easy to get that off-by-one wrong).
switch lower(kernel)
    case 'none',        Lspan = 0;
    case {'zoh','tri'},  Lspan = 2;
    case 'sinc',         Lspan = Lsinc;
    otherwise
        error('otfs_build_G:kernel','kernel must be none/zoh/tri/sinc');
end
fl = floor(ell(:).');
ls = [];
for x = fl
    ls = [ls, (x-Lspan):(x+Lspan+1)]; %#ok<AGROW>
end
ls = unique(ls);
end
