function [P, P_intuitive] = gen_perm_matrix(M, N)
% =========================================================================
% Program : gen_perm_matrix.m
% Description :
%   Build the NM x NM row-column interleaver (permutation) matrix P of
%   (4.33)/(4.34), in TWO independent ways so that they can be cross-checked
%   against each other (chap4_list.md Sim-3, Appendix-C code 3).
%
%   What P does
%   -----------
%   Let u = (I_M kron F_N') * x  be indexed by  j = n + m*N  (Doppler fast)
%   and let s be indexed by      q = m + n*M  (delay fast).  Then
%
%       s = P * u        i.e.   s[m + n*M] = u[n + m*N]
%       P[m + n*M , n + m*N] = 1,   all other entries 0
%
%   so that (4.32) reads  s = P * (G_tx kron F_N') * x,  and the receiver
%   uses  x_tilde = P.' * s   (4.53).  P is orthogonal: P.'*P = I_NM.
%
%   Method 1 - book form (4.33)(4.34)
%       P = [ E_11 E_21 ... E_M1 ;
%             E_12 E_22 ... E_M2 ;
%              :              :
%             E_1N E_2N ... E_MN ]
%       where E_mn is the M x N matrix with a single 1 at position (m,n).
%       (N block-rows) x (M block-cols), each block M x N  ->  NM x NM.
%
%   Method 2 - "write column-wise into an N x M matrix, read out row-wise"
%       This is the intuitive row-column interleaver picture of Fig. 4.10.
%       Applying that operation to each column of I_NM gives P directly.
%
%   Both methods must produce the SAME matrix (Sim-3 pass criterion).
%
% Input  : M - delay grid size,  N - Doppler grid size
% Output : P           - NM x NM permutation matrix (book form 4.33/4.34)
%          P_intuitive - NM x NM permutation matrix (interleaver form)
% =========================================================================

NM = N*M;

% ---------------------------------------------------------------------
% Method 1 : literal block construction from (4.33)(4.34)
% ---------------------------------------------------------------------
P = zeros(NM, NM);
for m = 1:M
    for n = 1:N
        E = zeros(M, N);
        E(m, n) = 1;                                   % E_mn
        rows = (n-1)*M + (1:M);                        % block-row n
        cols = (m-1)*N + (1:N);                        % block-col m
        P(rows, cols) = E;
    end
end

% ---------------------------------------------------------------------
% Method 2 : row-column interleaver applied to the columns of I_NM
%   u -> reshape(u, N, M) writes u column-wise into an N x M array
%        (i.e. entry (n,m) holds u[n + m*N]);
%        transposing and re-vectorizing reads it out row-wise, giving
%        index m + n*M.  That is exactly the action of P.
% ---------------------------------------------------------------------
P_intuitive = zeros(NM, NM);
I = eye(NM);
for j = 1:NM
    U = reshape(I(:, j), N, M);                        % write column-wise
    P_intuitive(:, j) = reshape(U.', NM, 1);           % read   row-wise
end

end
