function [r_noisy, w] = add_awgn(r, sigma2)
% =========================================================================
% Program : add_awgn.m
% Description :
%   Add circularly-symmetric complex AWGN to a time-domain vector
%   (chap4_list.md Sec. A.3 (4.36): r[q] = ... + w[q], w ~ CN(0, sigma_w^2)).
%
%   sigma2 is the noise variance PER COMPLEX SAMPLE, E|w|^2 = sigma2, split
%   evenly between the real and imaginary parts.  To relate sigma2 to a
%   target SNR = E_s/N_0 (Sec. B.1), use  sigma2 = Es / 10^(EsN0_dB/10).
%
% Input  : r      - NM x 1 (or any size) noiseless signal
%          sigma2 - noise variance per complex sample
% Output : r_noisy - r + w
%          w        - the noise realization actually added (same size as r)
% =========================================================================

sz = size(r);
w = sqrt(sigma2/2) * (randn(sz) + 1j*randn(sz));
r_noisy = r + w;

end
