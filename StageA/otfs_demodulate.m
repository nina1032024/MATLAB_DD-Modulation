function [Y, aux] = otfs_demodulate(r, M, N, method, P)
% =========================================================================
% Program : otfs_demodulate.m
% Description :
%   OTFS demodulator: time-domain receive vector r (NM x 1) -> DD-domain
%   matrix Y (M x N).  Four MATHEMATICALLY EQUIVALENT implementations,
%   mirroring otfs_modulate.m.  (chap4_list.md Sec. A.4, Appendix-C code 12.)
%
%   All paths start by undoing the Heisenberg/vec step (G_rx = I_M):
%       Y_tilde = reshape(r, M, N)      % delay-time domain, q = m + n*M
%
%   'wigner' Path A - Wigner transform + SFFT (Fig. 4.2)
%              Y_tf = F_M * Y_tilde                 (4.24)  Wigner
%              Y    = F_M' * Y_tf * F_N             (4.25)  SFFT
%
%   'dzt'    Path B - forward discrete Zak transform (Fig. 4.4)
%              Y    = Y_tilde * F_N                 (4.26)
%
%   'elem'   Path C - element-wise form (4.27), pure for-loops
%              Y[m,n] = (1/sqrt(N)) sum_p Y_tilde[m,p]*exp(-j*2*pi*n*p/N)
%
%   'vec'    Path D - vector form (4.30)
%              y = (I_M kron F_N) * P.' * r
%              Y = reshape(y, N, M).'      % undo x = vec(X.')
%
% Input  : r      - NM x 1 time-domain receive vector
%          M, N   - delay / Doppler grid sizes
%          method - 'wigner' | 'dzt' | 'elem' | 'vec'   (default 'dzt')
%          P      - NM x NM permutation matrix, only needed for 'vec';
%                   built internally if omitted.
% Output : Y      - M x N DD-domain received symbol matrix
%          aux    - .Y_tilde  M x N delay-time domain matrix
%                   .Y_tf     M x N time-frequency matrix, 'wigner' only
%                   .y        NM x 1 DD vector, y[n+mN] = Y[m,n], 'vec' only
% =========================================================================

if nargin < 4 || isempty(method), method = 'dzt'; end

r = r(:);
if numel(r) ~= N*M
    error('otfs_demodulate:size', 'r must have N*M = %d elements (got %d).', N*M, numel(r));
end

F_M = unitary_dft(M);
F_N = unitary_dft(N);

Y_tilde = reshape(r, M, N);        % delay-time domain (G_rx = I_M)
aux = struct();

switch lower(method)

    % ---------------------------------------------------------------
    case 'wigner'       % Path A : Wigner (4.24) then SFFT (4.25)
    % ---------------------------------------------------------------
        Y_tf = F_M * Y_tilde;                  % time-frequency domain
        Y    = F_M' * Y_tf * F_N;              % SFFT back to DD
        aux.Y_tf = Y_tf;

    % ---------------------------------------------------------------
    case 'dzt'          % Path B : forward DZT (4.26)
    % ---------------------------------------------------------------
        Y = Y_tilde * F_N;

    % ---------------------------------------------------------------
    case 'elem'         % Path C : element-wise (4.27)
    % ---------------------------------------------------------------
        Y = zeros(M, N);
        for m = 0:M-1
            for n = 0:N-1
                acc = 0;
                for pp = 0:N-1
                    acc = acc + Y_tilde(m+1, pp+1) * exp(-1j*2*pi*n*pp/N);
                end
                Y(m+1, n+1) = acc / sqrt(N);
            end
        end

    % ---------------------------------------------------------------
    case 'vec'          % Path D : vector form (4.30)
    % ---------------------------------------------------------------
        if nargin < 5 || isempty(P)
            P = gen_perm_matrix(M, N);
        end
        G_rx = eye(M);                          % rectangular pulse
        y = kron(G_rx, F_N) * (P.' * r);        % y[n+mN] = Y[m,n]
        Y = reshape(y, N, M).';                 % undo (4.31)
        aux.y = y;

    otherwise
        error('otfs_demodulate:method', ...
              'Unknown method ''%s''. Use wigner | dzt | elem | vec.', method);
end

aux.Y_tilde = Y_tilde;

end
