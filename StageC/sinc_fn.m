function y = sinc_fn(x)
% =========================================================================
% Program : sinc_fn.m
% Description :
%   Normalized sinc function, sinc(x) = sin(pi*x)/(pi*x), sinc(0) = 1.
%   Implemented locally (instead of calling MATLAB's `sinc`) so Stage 3 does
%   not depend on the Signal Processing Toolbox — same convention as
%   unitary_dft.m avoiding `dftmtx` in Stage 1.
%
% Input  : x - array of any size (real)
% Output : y - sinc(x), same size as x
% =========================================================================

y = ones(size(x));
nz = (x ~= 0);
y(nz) = sin(pi*x(nz)) ./ (pi*x(nz));

end
