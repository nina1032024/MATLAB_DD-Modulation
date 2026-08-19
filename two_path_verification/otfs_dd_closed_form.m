function Y = otfs_dd_closed_form(X, p, mode)
%OTFS_DD_CLOSED_FORM  Path B: compute Y_DD directly from the closed-form channel relation
%
%   Y = otfs_dd_closed_form(X, p, 'integer')     % Table 4.3 / (4.118)-(4.122)
%   Y = otfs_dd_closed_form(X, p, 'fractional')  % (4.105)/(4.78) with the periodic Dirichlet kernel
%
%   'integer'    : assumes delay/Doppler fall exactly on the grid, rounding
%                  ell/kappa to the nearest integer. This is the "approximate"
%                  version -- the gap to ground truth is exactly the leakage.
%   'fractional' : fully expands using the periodic Dirichlet kernel D_N(.);
%                  in theory should match ground truth exactly (up to the
%                  delay-kernel's own truncation/support error). Currently
%                  implemented for CP / ZP only.
%
%   p.kernel     : 'none' | 'zoh' | 'tri' | 'sinc'  (default 'tri'), same
%                  meaning as in otfs_build_G -- must match the ground
%                  truth pulse shape for the fractional mode to be exact.
%
%   NOTE - important correction: the book's (4.118)/(4.105) write the phase
%     as z^{kappa_i(m-l)} using the *tap index* l, but the correct phase
%     (from S2 of the delay-time channel derivation) uses the *true* delay
%     ell_i, i.e. z^{kappa_i(m-ell_i)}. The two coincide when ell_i is an
%     integer (l == ell_i), which is why this was invisible on
%     integer-delay tests, but they differ by ~30 dB once ell_i is
%     fractional. This function uses the corrected z^{kappa_i(m-ell_i)}.

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

if isfield(p,'kernel') && ~isempty(p.kernel)
    kernel = p.kernel;
else
    kernel = 'tri';
end

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
    taps = kernel_taps(ell, kernel, p.Lsinc);
    taps = taps(taps <= p.lmax);
    for m = 0:M-1
        for l = taps
            if l < 0, continue; end                        % CP/ZP delay-time channel is causal in l
            if strcmpi(p.variant,'ZP') && m < l, continue; end
            sw = eval_kernel(l-ell, kernel, p.Lsinc);
            if all(sw==0), continue; end
            % nu_{m,l}[k]  (4.105), 1 x N -- phase uses the TRUE delay
            % ell_i (not the tap index l); see the F4 note above.
            nu = (1/N) * ( ( g .* z.^(kap.*(m-ell)) .* sw ).' * dirichlet_N(kap*gam - k, N) );
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
%             Deliberately over-covers a couple of taps and lets the
%             kernel evaluation zero out anything outside the true
%             support, rather than hand-deriving floor/ceil offsets per
%             kernel (e.g. 'zoh' needs ceil(ell_i), not floor(ell_i)).
switch lower(kernel)
    case 'none',        Lspan = 0;
    case {'zoh','tri'},  Lspan = 2;
    case 'sinc',         Lspan = Lsinc;
    otherwise
        error('otfs_dd_closed_form:kernel','kernel must be none/zoh/tri/sinc');
end
fl = floor(ell(:).');
ls = [];
for x = fl
    ls = [ls, (x-Lspan):(x+Lspan+1)]; %#ok<AGROW>
end
ls = unique(ls);
end
