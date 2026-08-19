function p = delay_kernel(x, name)
%DELAY_KERNEL  Delay-domain composite pulse kernel p(x) = p(l - ell_i).
%
%   name = 'none' : integer delay, no interpolation at all (Kronecker delta)
%   name = 'zoh'  : rect Tx pulse + ideal (zero-hold) sampler, 1 tap
%   name = 'tri'  : rect Tx pulse + matched rect Rx filter, 2 taps
%   name = 'sinc' : ideal band-limited Tx/Rx pulse pair (current default
%                   ground truth), infinite support -- callers are
%                   responsible for truncating |x| beyond p.Lsinc.

switch lower(name)
    case 'none'
        p = double(abs(x) < 1e-12);
    case 'zoh'
        p = double(x >= 0 & x < 1);
    case 'tri'
        p = max(0, 1 - abs(x));
    case 'sinc'
        p = sin(pi*x)./(pi*x);
        p(x==0) = 1;
        p(abs(x-round(x))<1e-12 & x~=0) = 0;   % exact zero at nonzero integers, avoid 0/0 float noise
    otherwise
        error('delay_kernel:name','name must be one of none / zoh / tri / sinc');
end
end
