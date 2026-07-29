function F = unitary_dft(K)
% =========================================================================
% Program : unitary_dft.m
% Description :
%   Build the K-by-K UNITARY (normalized) DFT matrix used throughout
%   Chapter 4:
%
%       F_K[a,b] = (1/sqrt(K)) * exp(-j*2*pi*a*b/K),   a,b = 0...K-1
%
%   Properties: F_K' * F_K = I_K  (unitary), F_K.' = F_K (symmetric).
%
%   Why not dftmtx()?  dftmtx() needs the Signal Processing Toolbox and
%   returns the UN-normalized DFT matrix.  Book (Sec. 4.3) uses the unitary
%   version everywhere; mixing the two is the single most common source of
%   sqrt(M)/sqrt(N) errors (see chap4_list.md, Sec. H.1).
%
% Input  : K - transform size (positive integer)
% Output : F - K x K unitary DFT matrix (complex)
% =========================================================================

if ~isscalar(K) || K < 1 || K ~= floor(K)
    error('unitary_dft:badSize', 'K must be a positive integer scalar.');
end

a = (0:K-1).';           % row index  (K x 1)
b = (0:K-1);             % col index  (1 x K)
F = exp(-1j*2*pi*(a*b)/K) / sqrt(K);

end
