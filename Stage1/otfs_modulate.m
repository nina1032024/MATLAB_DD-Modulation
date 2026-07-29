function [s, aux] = otfs_modulate(X, method, P)
% =========================================================================
% Program : otfs_modulate.m
% Description :
%   OTFS modulator: DD-domain symbol matrix X (M x N) -> time-domain
%   transmit vector s (NM x 1).  Four MATHEMATICALLY EQUIVALENT
%   implementations are provided; Sim-2 / Sim-3 check that they agree.
%   (chap4_list.md Sec. A.2, Appendix-C code 3.)
%
%   'isfft'  Path A - Fig. 4.1 (ISFFT + Heisenberg)
%              X_tf    = F_M * X * F_N'          (4.16)  ISFFT
%              X_tilde = F_M' * X_tf             (4.17)  Heisenberg (G_tx=I)
%              s       = X_tilde(:)              (4.18)  vec, delay fastest
%
%   'idzt'   Path B - Fig. 4.3 (inverse discrete Zak transform)
%              X_tilde = X * F_N'                (4.17)  since F_M'*F_M = I
%              s       = X_tilde(:)              (4.18)
%
%   'elem'   Path C - element-wise form (4.22), pure for-loops
%              s[m+nM] = (1/sqrt(N)) * sum_p X[m,p] * exp(j*2*pi*p*n/N)
%
%   'vec'    Path D - vector/matrix form (4.32)
%              x = vec(X.')                      (4.31)  Doppler fastest
%              s = P * (G_tx kron F_N') * x,   G_tx = I_M  (rectangular pulse)
%
%   *** Output is ALWAYS an NM x 1 COLUMN vector. ***
%   (Appendix-C code 3 Method 1 returns a 1 x NM row vector, which is why
%    the book later has to write `G*s.'`; see chap4_list.md Sec. E.3.3.
%    Keeping one convention avoids a whole class of bugs.)
%
% Input  : X      - M x N DD-domain symbol matrix
%          method - 'isfft' | 'idzt' | 'elem' | 'vec'   (default 'idzt')
%          P      - NM x NM permutation matrix, only needed for 'vec';
%                   built internally if omitted.
% Output : s      - NM x 1 time-domain transmit vector, s[q], q = m + n*M
%          aux    - .X_tilde  M x N delay-time domain matrix (4.17)
%                   .X_tf     M x N time-frequency matrix    (4.16), 'isfft' only
%                   .x        NM x 1 DD vector (4.31),               'vec'   only
% =========================================================================

if nargin < 2 || isempty(method), method = 'idzt'; end

[M, N] = size(X);
F_M = unitary_dft(M);
F_N = unitary_dft(N);

aux = struct();

switch lower(method)

    % ---------------------------------------------------------------
    case 'isfft'        % Path A : ISFFT (4.16) then Heisenberg (4.17)(4.18)
    % ---------------------------------------------------------------
        X_tf    = F_M * X * F_N';        % time-frequency domain, M x N
        X_tilde = F_M' * X_tf;           % delay-time  domain,    M x N
        s       = X_tilde(:);            % column-wise vec: q = m + n*M
        aux.X_tf = X_tf;

    % ---------------------------------------------------------------
    case 'idzt'         % Path B : IDZT directly (4.17)
    % ---------------------------------------------------------------
        X_tilde = X * F_N';              % F_M' * F_M = I_M cancels out
        s       = X_tilde(:);

    % ---------------------------------------------------------------
    case 'elem'         % Path C : element-wise (4.22), explicit loops
    % ---------------------------------------------------------------
        s = zeros(N*M, 1);
        for m = 0:M-1
            for n = 0:N-1
                acc = 0;
                for pp = 0:N-1
                    acc = acc + X(m+1, pp+1) * exp(1j*2*pi*pp*n/N);
                end
                s(m + n*M + 1) = acc / sqrt(N);
            end
        end
        X_tilde = reshape(s, M, N);

    % ---------------------------------------------------------------
    case 'vec'          % Path D : vector form (4.32)
    % ---------------------------------------------------------------
        if nargin < 3 || isempty(P)
            P = gen_perm_matrix(M, N);
        end
        x = reshape(X.', N*M, 1);        % (4.31) x[n+mN] = X[m,n]
        G_tx = eye(M);                   % rectangular pulse (Sec. A.2)
        s = P * kron(G_tx, F_N') * x;    % (4.32)
        X_tilde = reshape(s, M, N);
        aux.x = x;

    otherwise
        error('otfs_modulate:method', ...
              'Unknown method ''%s''. Use isfft | idzt | elem | vec.', method);
end

aux.X_tilde = X_tilde;

end
