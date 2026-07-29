function p = otfs_params(varargin)
% =========================================================================
% Program : otfs_params.m
% Description :
%   Build the OTFS frame parameter struct for Chapter 4 simulations
%   (chap4_list.md Sec. B.1, Appendix-C code 1 & 4).
%
%   Default values are the VERIFICATION-stage values (M = 8, N = 6), which
%   match Fig. 4.6 / 4.7 of the book and are small enough that the matrix
%   structures (G, H_tilde, H) are still visible with imagesc().
%
% Usage :
%   p = otfs_params();                       % M = 8, N = 6 (verification)
%   p = otfs_params('M',64,'N',16);          % performance-stage size
%   p = otfs_params('M',3,'N',4,'seed',1);   % Fig. 4.10 interleaver demo
%
% Name-value options :
%   'M'       delay  grid size                            (default 8)
%   'N'       Doppler grid size                           (default 6)
%   'df'      subcarrier spacing Delta_f [Hz]             (default 15e3)
%   'fc'      carrier frequency [Hz]                      (default 4e9)
%   'modOrder' constellation order Q (4 = QPSK)           (default 4)
%   'seed'    RNG seed for reproducibility                (default 0)
%
% Output fields (in addition to the above) :
%   p.NM       = N*M            total frame length (samples)
%   p.T        = 1/df           block duration [s]
%   p.B        = M*df           bandwidth [Hz]
%   p.Ts       = 1/B = T/M      sample interval [s]
%   p.Tf       = N*M*Ts = N*T   frame duration [s]
%   p.d_res    = 1/(M*df)       delay resolution [s]
%   p.v_res    = 1/(N*T)        Doppler resolution [Hz]
%   p.F_M,p.F_N unitary DFT matrices  (see unitary_dft.m)
%   p.tol      = 1e-12          numerical pass threshold (Sec. B.3)
%   p.c        = 3e8            speed of light [m/s]
% =========================================================================

ip = inputParser;
ip.addParameter('M',        8,     @(v) isscalar(v) && v>0 && v==floor(v));
ip.addParameter('N',        6,     @(v) isscalar(v) && v>0 && v==floor(v));
ip.addParameter('df',       15e3,  @isscalar);
ip.addParameter('fc',       4e9,   @isscalar);
ip.addParameter('modOrder', 4,     @isscalar);
ip.addParameter('seed',     0,     @isscalar);
ip.addParameter('tol',      1e-12, @isscalar);
ip.parse(varargin{:});
o = ip.Results;

p.M        = o.M;
p.N        = o.N;
p.NM       = o.N * o.M;
p.df       = o.df;
p.fc       = o.fc;
p.modOrder = o.modOrder;
p.seed     = o.seed;
p.tol      = o.tol;
p.c        = 3e8;

% --- derived time/frequency quantities (Sec. B.1) ---------------------
p.T     = 1 / p.df;                 % block duration
p.B     = p.M * p.df;               % bandwidth
p.Ts    = 1 / p.B;                  % sample interval (= T/M)
p.Tf    = p.NM * p.Ts;              % frame duration  (= N*T)
p.d_res = 1 / (p.M * p.df);         % delay resolution
p.v_res = 1 / (p.N * p.T);          % Doppler resolution

% --- unitary DFT matrices (Sec. B.4) ----------------------------------
p.F_M = unitary_dft(p.M);
p.F_N = unitary_dft(p.N);

% --- z = exp(j*2*pi/(NM)); NOTE the denominator is NM, not N (Sec. H.3)
p.z = exp(1j*2*pi/p.NM);

end
