function [g, l, k, info] = gen_channel_taps(p, varargin)
% =========================================================================
% Program : gen_channel_taps.m
% Description :
%   Generate an ON-GRID multipath channel (g_i, l_i, k_i), i = 1...P
%   (chap4_list.md Sec. B.2, Appendix-C code 7 with k_i forced to integers
%   per Sec. E.3.4: "第一階段請加 k_i = round(k_i)").
%
%   Default taps match Sec. B.2's verification table:
%       l_i = [0, 1, 2]              (integer delay taps)
%       k_i = [0, 1, -2]             (integer Doppler taps, may be negative)
%       g_i : equal power, |g_i| = 1/sqrt(P), sum|g_i|^2 = 1, random phase
%
%   *** k_i is kept SIGNED (not wrapped mod N) ***
%   The time-domain phase z^{k_i(q-l)} in (4.8) depends on the RAW k_i
%   (z has period NM, not N), so wrapping k_i early would silently corrupt
%   every downstream matrix (gs, G, H-check, H-tilde, H). Only when a
%   symbol-domain index such as X[...,[n-k_i]_N] is taken (e.g. in
%   theory_dd_relation.m) does k_i need mod(.,N) — that wrapping happens
%   locally at the point of use, never here.
%
% Usage :
%   [g,l,k] = gen_channel_taps(p);                     % default 3-path channel
%   [g,l,k] = gen_channel_taps(p,'l',3,'k',0);          % single tap, on axis
%   [g,l,k] = gen_channel_taps(p,'l',[0 1 2],'k',[0 0 0]); % static (no Doppler)
%
% Name-value options :
%   'l'            delay taps l_i (integers)              (default [0 1 2])
%   'k'            Doppler taps k_i (integers, may be < 0) (default [0 1 -2])
%   'powerProfile' 'equal' | 'exp'                         (default 'equal')
%   'grid'         'on' | 'off' -> 'on' rounds k to integers (default 'on')
%   'randomPhase'  true/false -> random uniform path phase  (default true)
%
% Input  : p - parameter struct from otfs_params()
% Output : g    - 1xP complex path gains, sum|g|^2 = 1
%          l    - 1xP delay taps (integers, 0-indexed)
%          k    - 1xP Doppler taps (signed integers)
%          info - .l_max .k_max .numPaths  (+ warnings if constraints violated)
% =========================================================================

ip = inputParser;
ip.addParameter('l', [0 1 2]);
ip.addParameter('k', [0 1 -2]);
ip.addParameter('powerProfile', 'equal');
ip.addParameter('grid', 'on');
ip.addParameter('randomPhase', true);
ip.parse(varargin{:});
o = ip.Results;

l = o.l(:).';
k = o.k(:).';
if numel(l) ~= numel(k)
    error('gen_channel_taps:size', 'l and k must have the same number of paths.');
end
if strcmpi(o.grid, 'on')
    k = round(k);
end
Np = numel(l);

switch lower(o.powerProfile)
    case 'equal'
        mag = ones(1, Np) / sqrt(Np);
    case 'exp'
        decay = exp(-(0:Np-1));
        mag = sqrt(decay / sum(decay));
    otherwise
        error('gen_channel_taps:powerProfile', 'Unknown powerProfile ''%s''.', o.powerProfile);
end

if o.randomPhase
    rng(p.seed + 100);                 % separate stream from symbol RNG
    ph = 2*pi*rand(1, Np);
else
    ph = zeros(1, Np);
end
g = mag .* exp(1j*ph);

info.numPaths = Np;
info.l_max    = max(l);
info.k_max    = max(abs(k));

if info.l_max >= p.M
    warning('gen_channel_taps:lmax', 'l_max (%d) should be < M (%d).', info.l_max, p.M);
end
if info.k_max >= p.N/2
    warning('gen_channel_taps:kmax', ...
            'k_max (%d) should be < N/2 (%.1f) to avoid Doppler aliasing.', info.k_max, p.N/2);
end

end
