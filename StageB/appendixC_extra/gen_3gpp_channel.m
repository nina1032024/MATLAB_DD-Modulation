function [g, l, k, info] = gen_3gpp_channel(p, model, k_max, seed)
% =========================================================================
% Program : gen_3gpp_channel.m
% Description :
%   Generate channel taps (g_i, l_i, k_i) from a 3GPP standard
%   power-delay-profile (Appendix C, MATLAB code 5 + code 6 — a capability
%   completely absent from Stage 1-3 before this file: every Sim so far
%   used only synthetic on-grid/off-grid taps, never a real standard model).
%
%   Delay profiles (code 5), delays in seconds, PDP in dB:
%       EPA : 7 taps, delays 0...410 ns   (pedestrian, small delay spread)
%       EVA : 9 taps, delays 0...2510 ns  (vehicular, medium delay spread)
%       ETU : 9 taps, delays 0...5000 ns  (urban, large delay spread)
%
%   Processing (code 6), applied exactly as in the book:
%       1. dB -> linear, normalize so sum(pdp_linear) = 1
%       2. complex Rayleigh gain per path: g_i = sqrt(pdp_i) * CN(0,1)
%       3. delay tap l_i = round(delay_i / delay_resolution)   -- ON-GRID
%          (quantizing to the OTFS delay grid; see Sec. H point 5 / E.3.5:
%          for a small M, several distinct physical delays can collapse
%          onto the SAME l_i — quantified in sim_appendixC_channel_models.m)
%       4. Doppler tap k_i = k_max * cos(2*pi*rand(1,taps))   -- Jakes
%          spectrum (classic Clarke/Jakes arcsine-density Doppler model).
%          k_i is FRACTIONAL by construction (unless it happens to land on
%          an integer) — gen_gs.m (Stage 2) handles this fine since only
%          l_i needs to be an integer for its row-indexing, not k_i.
%
% Input  : p     - parameter struct from otfs_params()
%          model - 'EPA' | 'EVA' | 'ETU'
%          k_max - one-sided max normalized Doppler (default: from
%                  mobility_params(p), i.e. 500 km/h at p.fc)
%          seed  - RNG seed (default p.seed + 200, separate stream from
%                  the DD symbol / synthetic-tap generators)
% Output : g    - 1 x taps complex Rayleigh path gains
%          l    - 1 x taps integer delay taps (on-grid, may repeat!)
%          k    - 1 x taps Doppler taps (fractional, Jakes spectrum)
%          info - .delays_s, .pdp_dB, .pdp_linear, .taps, .l_raw (before
%                 rounding), .distinct_l (how many UNIQUE l_i survive
%                 quantization — the Sec. H.5 concern, quantified here)
% =========================================================================

if nargin < 3 || isempty(k_max)
    mp = mobility_params(p);
    k_max = mp.k_max;
end
if nargin < 4 || isempty(seed), seed = p.seed + 200; end

switch upper(model)
    case 'EPA'
        delays = [0, 30, 70, 90, 110, 190, 410] * 1e-9;
        pdp_dB = [0.0, -1.0, -2.0, -3.0, -8.0, -17.2, -20.8];
    case 'EVA'
        delays = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510] * 1e-9;
        pdp_dB = [0.0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];
    case 'ETU'
        delays = [0, 50, 120, 200, 230, 500, 1600, 2300, 5000] * 1e-9;
        pdp_dB = [-1.0, -1.0, -1.0, 0.0, 0.0, 0.0, -3.0, -5.0, -7.0];
    otherwise
        error('gen_3gpp_channel:model', 'Unknown model ''%s''. Use EPA | EVA | ETU.', model);
end

pdp_linear = 10 .^ (pdp_dB / 10);
pdp_linear = pdp_linear / sum(pdp_linear);
taps = numel(pdp_dB);

rng(seed);
g = sqrt(pdp_linear) .* (sqrt(1/2) * (randn(1, taps) + 1j*randn(1, taps)));

l_raw = delays / p.d_res;
l = round(l_raw);

k = k_max * cos(2*pi*rand(1, taps));

info.delays_s    = delays;
info.pdp_dB      = pdp_dB;
info.pdp_linear  = pdp_linear;
info.taps        = taps;
info.l_raw       = l_raw;
info.distinct_l  = numel(unique(l));
info.k_max       = k_max;

end
