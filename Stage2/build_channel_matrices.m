function C = build_channel_matrices(p, G, P)
% =========================================================================
% Program : build_channel_matrices.m
% Description :
%   From the time-domain channel matrix G, build the channel matrix in
%   every other domain defined in chap4_list.md Sec. A.4 / Table row 66-71:
%
%       time-frequency   : Ȟ = (I_N kron F_M) G (I_N kron F_M)^†      (4.43)
%       delay-time       : H̃ = P^T G P                                (4.55, corrected)
%       delay-Doppler    : H  = (I_M kron F_N) H̃ (I_M kron F_N)^†     (4.60)
%
%   *** Correction vs. Appendix-C code 9 (chap4_list.md Sec. E.3.2) ***
%   The book's reference code computes `H_tilde = P*G*P.'`, but the
%   receiver actually forms x_tilde = P^T * s (4.53) and r_tilde = P^T * r,
%   so consistency requires H̃ = P^T G P, NOT P G P^T. Sim-8 numerically
%   confirms which one is right by checking ỹ = H̃ x̃ directly; this file
%   implements the corrected version.
%
% Input  : p - parameter struct (needs p.M, p.N, p.F_M, p.F_N)
%          G - NM x NM time-domain channel matrix (from gen_G.m)
%          P - NM x NM permutation matrix (from gen_perm_matrix.m)
% Output : C - struct with fields
%              .Hcheck  NM x NM  Ȟ, time-frequency domain channel (4.43)
%              .Htilde  NM x NM  H̃, delay-time domain channel (4.55, corrected)
%              .H       NM x NM  H,  delay-Doppler domain channel (4.60)
%              .IN_FM   NM x NM  (I_N kron F_M), the time -> time-freq operator
%              .IM_FN   NM x NM  (I_M kron F_N), the delay-time -> DD operator
% =========================================================================

M = p.M;  N = p.N;
F_M = p.F_M;  F_N = p.F_N;

IN_FM = kron(eye(N), F_M);
IM_FN = kron(eye(M), F_N);

C.IN_FM  = IN_FM;
C.IM_FN  = IM_FN;
C.Hcheck = IN_FM * G * IN_FM';
C.Htilde = P.' * G * P;
C.H      = IM_FN * C.Htilde * IM_FN';

end
