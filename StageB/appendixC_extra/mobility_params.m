function info = mobility_params(p, max_speed_kmh)
% =========================================================================
% Program : mobility_params.m
% Description :
%   Compute the maximum Doppler spread from a maximum UE speed
%   (Appendix C, MATLAB code 4 — not yet used anywhere in Stage 1-3 before
%   this file; only referenced indirectly via B.2's "最大車速 500 km/h"
%   note). Matches the book's code exactly:
%
%       max_UE_speed = max_UE_speed_kmh * (1000/3600)
%       nu_max       = max_UE_speed * fc / c            (one-sided Doppler, Hz)
%       k_max        = nu_max / Doppler_resolution      (normalized Doppler tap)
%
% Input  : p             - parameter struct from otfs_params()
%          max_speed_kmh - maximum UE speed [km/h] (default 500, per
%                          chap4_list.md Sec. B.2; the book's own code 4
%                          example uses 100 km/h)
% Output : info - struct with
%              .max_speed_kmh, .max_speed_mps, .nu_max_hz, .k_max
%              .aliasing_ok = (k_max < N/2)   -- Sec. H point 8 check
% =========================================================================

if nargin < 2 || isempty(max_speed_kmh), max_speed_kmh = 500; end

max_speed_mps = max_speed_kmh * (1000/3600);
nu_max_hz     = max_speed_mps * p.fc / p.c;
k_max         = nu_max_hz / p.v_res;

info.max_speed_kmh = max_speed_kmh;
info.max_speed_mps  = max_speed_mps;
info.nu_max_hz      = nu_max_hz;
info.k_max          = k_max;
info.aliasing_ok    = k_max < p.N/2;

end
