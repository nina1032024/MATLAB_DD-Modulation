function [X, x, info] = gen_dd_symbols(p)
% =========================================================================
% Program : gen_dd_symbols.m
% Description :
%   Generate the delay-Doppler (DD) domain information symbol matrix X and
%   its vectorized form x  (chap4_list.md Sec. A.2, Appendix-C code 2).
%
%       X[m,n] in A,     m = 0...M-1 (delay),  n = 0...N-1 (Doppler)
%       x = vec(X.'),    x[n + m*N] = X[m,n]                     (4.31)
%
%   *** The two vectorizations are DIFFERENT (Sec. H.2) ***
%       x : Doppler index runs fastest  -> reshape(X.', NM, 1)
%       s : delay   index runs fastest  -> X_tilde(:)   (see otfs_modulate)
%   The permutation matrix P exists precisely to convert between them.
%
%   Constellation is normalized to unit average symbol energy, E_s = 1.
%
% Input  : p    - parameter struct from otfs_params()
% Output : X    - M x N  DD-domain symbol matrix
%          x    - NM x 1 DD-domain symbol vector, x[n+mN] = X[m,n]
%          info - .alphabet  constellation points
%                 .sym_idx   M x N matrix of constellation indices
%                 .Es        measured average symbol energy
% =========================================================================

M = p.M;  N = p.N;
rng(p.seed);                                   % reproducibility (Sec. B.3)

switch p.modOrder
    case 4      % ---- QPSK, Gray-mapped, Es = 1 -------------------------
        alphabet = [ 1+1j, -1+1j, -1-1j, 1-1j ] / sqrt(2);
    case 16     % ---- 16-QAM, Es = 1 ------------------------------------
        lvl = [-3 -1 1 3];
        alphabet = reshape(lvl.' + 1j*lvl, 1, []) / sqrt(10);
    case 2      % ---- BPSK ----------------------------------------------
        alphabet = [1, -1];
    otherwise
        error('gen_dd_symbols:modOrder', 'Unsupported modOrder = %d.', p.modOrder);
end

sym_idx = randi(numel(alphabet), M, N);        % M x N indices
X       = alphabet(sym_idx);                   % M x N DD symbol matrix

% (4.31): x = vec(X.')  -> Doppler index runs fastest
x = reshape(X.', M*N, 1);

info.alphabet = alphabet;
info.sym_idx  = sym_idx;
info.Es       = mean(abs(X(:)).^2);

end
