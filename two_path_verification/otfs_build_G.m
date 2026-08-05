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
%   p.Lsinc       : sinc truncation half-width (only matters for fractional delay)
%   p.variant     : 'CP' | 'ZP' | 'RCP' | 'RZP'
%   p.gains       : P x 1  complex gains g_i
%   p.ell         : P x 1  normalized delay  ell_i = tau_i/Ts        (may be fractional)
%   p.kappa       : P x 1  normalized Doppler kappa_i = nu_i*N*T     (may be fractional)
%
%   Note: the CP/ZP variants return the NM x NM equivalent matrix "with the
%   guard removed", consistent with the book's (4.93)/(4.109).

M = p.M;  N = p.N;  NM = M*N;  Lg = p.Lg;  lmax = p.lmax;
z  = exp(1j*2*pi/NM);
g   = p.gains(:);
ell = p.ell(:);
kap = p.kappa(:);

% g^s[l,q] = sum_i g_i * z^{kappa_i (q-l)} * sinc(l - ell_i)     (4.6)
gs = @(l,q) sum( g .* z.^(kap*(q-l)) .* sincw(l-ell, p.Lsinc) );

G = zeros(NM,NM);

switch upper(p.variant)

    case 'RZP'   % one ZP per frame -> pure linear convolution, lower triangular
        for q = 0:NM-1
            for l = 0:min(lmax,q)
                G(q+1, q-l+1) = gs(l,q);
            end
        end

    case 'RCP'   % one CP per frame -> circular over the whole frame
        for q = 0:NM-1
            for l = 0:lmax
                c = mod(q-l,NM)+1;
                G(q+1,c) = G(q+1,c) + gs(l,q);
            end
        end

    case 'CP'    % one CP per block -> each M x M diagonal block is circular on its own
        for n = 0:N-1
            for m = 0:M-1
                q = m + n*(M+Lg);            % true time index including guard
                for l = 0:lmax
                    c = n*M + mod(m-l,M) + 1;
                    G(n*M+m+1,c) = G(n*M+m+1,c) + gs(l,q);
                end
            end
        end

    case 'ZP'    % one ZP per block -> each M x M diagonal block is lower triangular
        for n = 0:N-1
            for m = 0:M-1
                q = m + n*(M+Lg);
                for l = 0:min(lmax,m)
                    G(n*M+m+1, n*M+(m-l)+1) = gs(l,q);
                end
            end
        end

    otherwise
        error('otfs_build_G:variant','variant must be CP / ZP / RCP / RZP');
end
end

% ---------------------------------------------------------------
function y = sincw(x,L)
%SINCW  truncated sinc(x) = sin(pi x)/(pi x), set to 0 outside |x|>L
x = x(:);
y = sin(pi*x)./(pi*x);
y(x==0) = 1;
y(abs(x-round(x))<1e-12 & x~=0) = 0;   % force exact zero at integer points, avoiding leftover 0/0 numerical noise
y(abs(x)>L) = 0;
end
