function e = relerr(a, b)
% =========================================================================
% Program : relerr.m
% Description :
%   Relative error between two quantities, used as the common pass/fail
%   metric for every Sim in chap4_list.md Sec. B.3:
%
%       relerr(a,b) = ||a - b||_2 / ||b||_2        (vectorized, Frobenius)
%
%   If ||b|| is (numerically) zero the ABSOLUTE error ||a-b|| is returned
%   instead, so that comparisons against an all-zero reference still make
%   sense (e.g. checking that an off-diagonal block vanishes).
%
% Input  : a, b - arrays of identical size (real or complex)
% Output : e    - scalar relative (or absolute) error
% =========================================================================

if ~isequal(size(a), size(b))
    error('relerr:sizeMismatch', 'Inputs must have the same size: [%s] vs [%s].', ...
          num2str(size(a)), num2str(size(b)));
end

a = a(:);  b = b(:);
nb = norm(b);

if nb < eps
    e = norm(a - b);          % absolute error against a zero reference
else
    e = norm(a - b) / nb;     % relative error
end

end
