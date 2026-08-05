function Y = otfs_dd_closed_form(X, p, mode)
%OTFS_DD_CLOSED_FORM  Path B: compute Y_DD directly from the closed-form channel relation
%
%   Y = otfs_dd_closed_form(X, p, 'integer')     % Table 4.3 / (4.118)-(4.122)
%   Y = otfs_dd_closed_form(X, p, 'fractional')  % (4.105)/(4.78) with the periodic sinc
%
%   'integer'    : assumes delay/Doppler fall exactly on the grid, rounding
%                  ell/kappa to the nearest integer. This is the "approximate"
%                  version -- the gap to ground truth is exactly the leakage.
%   'fractional' : fully expands using the periodic sinc S_N(.); in theory
%                  should match ground truth exactly (up to sinc truncation
%                  error). Currently implemented for CP / ZP only.
%
%   NOTE - important correction: the book's (4.118) writes the phase as
%     z^{k_i(m-l_i)}, but (4.97) implies the correct phase should use the
%     *unscaled* kappa_i; only the Doppler shift index k_i needs to be
%     multiplied by gamma_g. This function uses the corrected version --
%     otherwise CP/ZP will not match.

M = p.M;  N = p.N;  NM = M*N;
z = exp(1j*2*pi/NM);
switch upper(p.variant)
    case {'CP','ZP'}, gam = 1 + p.Lg/M;   % (4.97) Doppler scaling
    otherwise,        gam = 1;
end
g   = p.gains(:);
ell = p.ell(:);
kap = p.kappa(:);
Y   = zeros(M,N);

switch lower(mode)

% ================= symbol-wise form for integer taps (Table 4.3) =================
case 'integer'
    li = round(ell);
    ki = round(kap*gam);
    for i = 1:numel(g)
        for m = 0:M-1
            mm = mod(m-li(i),M);
            for n = 0:N-1
                nn = mod(n-ki(i),N);
                switch upper(p.variant)
                    case 'CP'
                        a = z^(kap(i)*(m-li(i)));
                    case 'ZP'
                        if m >= li(i), a = z^(kap(i)*(m-li(i))); else, a = 0; end
                    case 'RCP'
                        if m >= li(i), a = z^(kap(i)*(m-li(i)));
                        else,          a = z^(kap(i)*mm) * exp(-1j*2*pi*n/N); end
                    case 'RZP'
                        if m >= li(i), a = z^(kap(i)*(m-li(i)));
                        else,          a = ((N-1)/N) * z^(kap(i)*(m-li(i))) * ...
                                           exp(-1j*2*pi*nn/N); end
                end
                Y(m+1,n+1) = Y(m+1,n+1) + g(i)*a*X(mm+1,nn+1);
            end
        end
    end

% ============ fractional: full Doppler spread vector (CP/ZP) ============
case 'fractional'
    if ~any(strcmpi(p.variant,{'CP','ZP'}))
        error('fractional mode currently supports CP / ZP only');
    end
    k = 0:N-1;
    for m = 0:M-1
        for l = 0:p.lmax
            if strcmpi(p.variant,'ZP') && m < l, continue; end
            sw = sincw(l-ell, p.Lsinc);
            if all(sw==0), continue; end
            % nu_{m,l}[k]  (4.105), 1 x N
            nu = (1/N) * ( ( g .* z.^(kap*(m-l)) .* sw ).' * SN(kap*gam - k, N) );
            mm = mod(m-l,M);
            for n = 0:N-1
                Y(m+1,n+1) = Y(m+1,n+1) + nu * X(mm+1, mod(n-k,N)+1).';
            end
        end
    end

otherwise
    error('mode must be integer or fractional');
end
end

% ---------------------------------------------------------------
function y = SN(x,N)
%SN  periodic sinc  S_N(x) = sum_{n=0}^{N-1} exp(j2*pi*x*n/N)      (4.79)
den   = sin(pi*x/N);
small = abs(den) < 1e-12;                    % x = 0 (mod N)
y = sin(pi*x)./(den + small) .* exp(1j*pi*x*(N-1)/N);
y(small) = N;                                % the limiting value is N, not sqrt(N)
end

% ---------------------------------------------------------------
function y = sincw(x,L)
x = x(:);
y = sin(pi*x)./(pi*x);
y(x==0) = 1;
y(abs(x-round(x))<1e-12 & x~=0) = 0;
y(abs(x)>L) = 0;
end
