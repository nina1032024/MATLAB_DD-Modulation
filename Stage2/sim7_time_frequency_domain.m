function results = sim7_time_frequency_domain(p, doPlot)
    % =========================================================================
    % Program : sim7_time_frequency_domain.m
    % Description :
    %   Sim-7  "time-frequency domain I/O relation (4.4.2)"
    %          (chap4_list.md Sec. D, Stage 2)
    %
    %   Purpose  : verify the time-frequency domain channel relation and
    %              quantify ICI (inter-carrier interference) / ISI
    %              (inter-symbol/inter-block interference) energy.
    %   Theory   : y-check = H-check x-check,
    %              H-check = (I_N kron F_M) G (I_N kron F_M)^†,
    %              H-check_{n,n'} = F_M G_{n,n'} F_M^†          (4.40-4.43)(4.46)
    %              Using the Sim-6 block notation: H-check_{n,0} = F_M G_{n,0} F_M^†
    %              (self block -> ICI, off-diagonal energy within this M x M
    %              block), H-check_{n,1} = F_M G_{n,1} F_M^† (one-block-back -> ISI).
    %   Method   : x-check=(I_N⊗F_M)s, y-check=(I_N⊗F_M)r, compare against
    %              H-check*x-check. Then, for interior block n=2, compute the
    %              ICI ratio = ||off-diag(H_{n,0})||^2 / ||H_{n,0}||^2 and the
    %              ISI energy = ||H_{n,1}||^2, sweeping k_max and l_max.
    %   Criterion: ||H-check x-check - y-check|| / ||y-check|| < 1e-10, plus the
    %              two analytic identities below.
    %
    %   -----------------------------------------------------------------------
    %   ANALYTIC CHECKS ADDED (these turn Sim-7 from "plot and eyeball" into a
    %   real verification):
    %
    %   (a) ISI energy is INVARIANT to Doppler.  F_M is unitary, so
    %           ||H-check_{n,1}||_F^2 = ||F_M G_{n,1} F_M^†||_F^2 = ||G_{n,1}||_F^2
    %       and Doppler only multiplies each entry of G by a unit-modulus phase,
    %       which cannot change |entry|.  Hence
    %           ||H-check_{n,1}||_F^2 = sum_i |g_i|^2 * l_i                    (ISI_THEORY)
    %       -> sweeping k_max produces a FLAT ISI curve by construction.  To make
    %          ISI move you must sweep the delay spread l_max, which is why the
    %          l_max sweep below was added.
    %
    %   (b) ICI ratio at k_max = 0 has a closed-form upper bound.  With no
    %       Doppler, G_{n,0} is a truncated Toeplitz matrix C - E, where C is the
    %       circulant completion and E holds the l_i entries missing from each
    %       tap's wrapped corner.  All of E's energy is off-diagonal, so
    %           ICI_ratio|_{k=0}  <=  ||E||_F^2 / ||G_{n,0}||_F^2
    %                              =  sum_i |g_i|^2 l_i / sum_i |g_i|^2 (M - l_i)
    %       This bound applies ONLY at k_max = 0; once Doppler is present
    %       G_{n,0} is no longer Toeplitz and the bound does not hold.
    %   -----------------------------------------------------------------------
    %
    %   Note on the "static channel -> H-check_{n,0} diagonal" observation
    %   (chap4_list.md Sec. D, Sim-7 bullet): the book ties this EXACTLY to
    %   diagonal blocks becoming circulant, which under Sec. 4.4's isolated-
    %   frame RZP/RCP models only happens at the FRAME boundary, not per-block
    %   -> true per-block circulant structure requires CP-OTFS (Sec. 4.5.3,
    %   Sim-17, Stage 4, deliberately NOT built here). What IS verified here:
    %   with k_i = 0 for all paths, the ICI ratio drops sharply (structured
    %   Toeplitz block) even though it is not exactly zero without a real CP.
    %
    % Input  : p      - parameter struct from otfs_params()   (optional)
    %          doPlot - true to draw imagesc(abs(H-check)) and the ICI/ISI curves
    %                   (default true)
    % Output : results - struct with .name .err_rel .tol .pass, plus
    %              .ici_ratio (1xN), .isi_energy (1xN), .ici_static,
    %              .kmax_list, .ici_curve, .isi_curve,
    %              .lmax_list, .ici_curve_l, .isi_curve_l,
    %              .isi_theory, .isi_theory_err, .ici_bound_k0, .checks
    % =========================================================================
    
    thisDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(thisDir, '..', 'Stage1'));
    
    if nargin < 1 || isempty(p),      p = otfs_params();  end
    if nargin < 2 || isempty(doPlot), doPlot = true;      end
    
    M = p.M;  N = p.N;  F_M = p.F_M;
    
    fprintf('\n=== Sim-7 : time-frequency domain I/O relation (4.4.2)  (M=%d, N=%d) ===\n', M, N);
    
    [g, l, k] = gen_channel_taps(p);
    gs = gen_gs(p, g, l, k);
    G  = gen_G(p, gs, 'RZP');
    
    X  = gen_dd_symbols(p);
    Pm = gen_perm_matrix(M, N);
    s  = otfs_modulate(X, 'idzt', Pm);
    r  = G * s;
    
    C = build_channel_matrices(p, G, Pm);
    x_check      = C.IN_FM * s;
    y_check      = C.IN_FM * r;
    y_check_pred = C.Hcheck * x_check;
    
    err_io = relerr(y_check_pred, y_check);
    fprintf('    ||Hcheck*x_check - y_check|| / ||y_check|| = %.3e\n', err_io);
    
    % ---- ICI / ISI energy per time-block ------------------------------------
    ici_ratio  = zeros(1, N);
    isi_energy = zeros(1, N);
    for n = 0:N-1
        Gn0 = get_block(G, M, n, n);
        Hn0 = F_M * Gn0 * F_M';
        offE  = norm(Hn0 - diag(diag(Hn0)), 'fro')^2;
        totE  = norm(Hn0, 'fro')^2;
        ici_ratio(n+1) = offE / max(totE, eps);
        if n >= 1
            Gn1 = get_block(G, M, n, n-1);
            Hn1 = F_M * Gn1 * F_M';
            isi_energy(n+1) = norm(Hn1, 'fro')^2;
        end
    end
    fprintf('    ICI ratio per block  (n=0..%d): %s\n', N-1, mat2str(ici_ratio, 3));
    fprintf('    ISI energy per block (n=0..%d): %s\n', N-1, mat2str(isi_energy, 3));
    
    % ---- ANALYTIC CHECK (a): ISI energy = sum_i |g_i|^2 * l_i ----------------
    isi_theory = sum(abs(g(:)).^2 .* l(:));
    isi_meas   = isi_energy(2:end);                      % blocks n >= 1 only
    isi_theory_err = max(abs(isi_meas - isi_theory)) / max(isi_theory, eps);
    fprintf('\n    [check a] ISI energy: theory sum|g_i|^2 l_i = %.6f\n', isi_theory);
    fprintf('              measured (n>=1) = %s\n', mat2str(isi_meas, 6));
    fprintf('              max rel. deviation = %.3e  -> %s\n', ...
            isi_theory_err, ternary(isi_theory_err < 1e-10, 'PASS', 'FAIL'));
    
    % ---- static channel (k_i = 0) : ICI drops but is not exactly 0 (no CP) --
    nBlockCheck = min(2, N-1);
    [g0, l0, k0] = gen_channel_taps(p, 'l', l, 'k', zeros(size(k)), 'randomPhase', false);
    gs0 = gen_gs(p, g0, l0, k0);
    G0  = gen_G(p, gs0, 'RZP');
    Gn0_static = get_block(G0, M, nBlockCheck, nBlockCheck);
    Hn0_static = F_M * Gn0_static * F_M';
    ici_static = norm(Hn0_static - diag(diag(Hn0_static)), 'fro') / norm(Hn0_static, 'fro');
    fprintf('\n    static channel (all k_i=0), ICI ratio at block n=%d : %.3e  (>0: needs CP-OTFS to vanish exactly)\n', ...
            nBlockCheck, ici_static);
    
    % ---- ANALYTIC CHECK (b): closed-form ICI bound at k = 0 -----------------
    ici_bound_k0 = sum(abs(g0(:)).^2 .* l0(:)) / sum(abs(g0(:)).^2 .* (M - l0(:)));
    ici_ratio_k0 = norm(Hn0_static - diag(diag(Hn0_static)), 'fro')^2 / norm(Hn0_static, 'fro')^2;
    fprintf('    [check b] ICI ratio at k=0 = %.6f,  analytic bound = %.6f  -> %s\n', ...
            ici_ratio_k0, ici_bound_k0, ...
            ternary(ici_ratio_k0 <= ici_bound_k0 * (1 + 1e-9), 'PASS (within bound)', 'FAIL (exceeds bound)'));
    
    % ---- sweep 1 : ICI/ISI vs k_max (delay profile fixed) --------------------
    % ISI is EXPECTED to be flat here; see analytic check (a).
    lFixed    = [0 1 2];
    kmax_list = 0:(floor(N/2)-1);
    ici_curve = zeros(size(kmax_list));
    isi_curve = zeros(size(kmax_list));
    for idx = 1:numel(kmax_list)
        kt = kmax_list(idx);
        if kt == 0
            kTest = zeros(1, numel(lFixed));
        else
            kTest = [0, kt, -kt];
        end
        [gT, lT, kT] = gen_channel_taps(p, 'l', lFixed, 'k', kTest, 'randomPhase', false);
        gsT = gen_gs(p, gT, lT, kT);
        GT  = gen_G(p, gsT, 'RZP');
        Hn0T = F_M * get_block(GT, M, nBlockCheck, nBlockCheck)   * F_M';
        Hn1T = F_M * get_block(GT, M, nBlockCheck, nBlockCheck-1) * F_M';
        ici_curve(idx) = norm(Hn0T - diag(diag(Hn0T)), 'fro')^2 / norm(Hn0T, 'fro')^2;
        isi_curve(idx) = norm(Hn1T, 'fro')^2;
    end
    fprintf('\n    k_max sweep = %s   (delay profile fixed at l = %s)\n', ...
            mat2str(kmax_list), mat2str(lFixed));
    fprintf('    ICI ratio   = %s\n', mat2str(ici_curve, 3));
    fprintf('    ISI energy  = %s   <- flat by construction, see check (a)\n', mat2str(isi_curve, 3));
    
    % ---- sweep 2 : ICI/ISI vs l_max (Doppler fixed) -------------------------
    % This is the sweep that actually moves ISI.
    kFixed      = [0 1 -1];
    lmax_list   = 1:min(M-1, 6);
    ici_curve_l = zeros(size(lmax_list));
    isi_curve_l = zeros(size(lmax_list));
    for idx = 1:numel(lmax_list)
        lm = lmax_list(idx);
        lTest = unique([0, floor(lm/2), lm]);
        kTest = kFixed(1:numel(lTest));
        [gL, lL, kL] = gen_channel_taps(p, 'l', lTest, 'k', kTest, 'randomPhase', false);
        gsL = gen_gs(p, gL, lL, kL);
        GL  = gen_G(p, gsL, 'RZP');
        Hn0L = F_M * get_block(GL, M, nBlockCheck, nBlockCheck)   * F_M';
        Hn1L = F_M * get_block(GL, M, nBlockCheck, nBlockCheck-1) * F_M';
        ici_curve_l(idx) = norm(Hn0L - diag(diag(Hn0L)), 'fro')^2 / norm(Hn0L, 'fro')^2;
        isi_curve_l(idx) = norm(Hn1L, 'fro')^2;
    end
    fprintf('\n    l_max sweep = %s   (Doppler fixed at k = %s)\n', ...
            mat2str(lmax_list), mat2str(kFixed));
    fprintf('    ICI ratio   = %s\n', mat2str(ici_curve_l, 3));
    fprintf('    ISI energy  = %s\n', mat2str(isi_curve_l, 3));
    
    % ---- results ------------------------------------------------------------
    checks.io_relation = err_io < 1e-10;
    checks.isi_theory  = isi_theory_err < 1e-10;
    checks.ici_bound   = ici_ratio_k0 <= ici_bound_k0 * (1 + 1e-9);
    
    results.name           = 'Sim-7 time-frequency domain I/O relation (4.4.2)';
    results.err_rel        = err_io;
    results.tol            = 1e-10;
    results.ici_ratio      = ici_ratio;
    results.isi_energy     = isi_energy;
    results.ici_static     = ici_static;
    results.kmax_list      = kmax_list;
    results.ici_curve      = ici_curve;
    results.isi_curve      = isi_curve;
    results.lmax_list      = lmax_list;
    results.ici_curve_l    = ici_curve_l;
    results.isi_curve_l    = isi_curve_l;
    results.isi_theory     = isi_theory;
    results.isi_theory_err = isi_theory_err;
    results.ici_bound_k0   = ici_bound_k0;
    results.ici_ratio_k0   = ici_ratio_k0;
    results.checks         = checks;
    results.pass           = checks.io_relation && checks.isi_theory && checks.ici_bound;
    
    fprintf('\n    -> %s\n', ternary(results.pass, 'PASS', 'FAIL'));
    
    % ---- plots --------------------------------------------------------------
    if doPlot
        % Fig. 3(b) : |H-check|
        figure('Name', 'Sim-7 |Hcheck|', 'Color', 'w');
        imagesc(abs(C.Hcheck)); axis square; colorbar;
        % NOTE: the vertical bars MUST be inside math mode.  In text mode with the
        % OT1 encoding, '|' is typeset as an em-dash, which is why the old title
        % rendered as "—H—".
        title('$|\check{H}|$ (time-frequency domain)', 'Interpreter', 'latex');
        xlabel('column index'); ylabel('row index');
    
        % Fig. 6 : ICI / ISI vs Doppler spread
        figure('Name', 'Sim-7 ICI/ISI vs k_{max}', 'Color', 'w');
        yyaxis left;
        plot(kmax_list, ici_curve, '-o', 'LineWidth', 1.5);
        ylabel('ICI ratio (block n=2)');
        ylim([0, max(ici_curve) * 1.2 + eps]);
        yyaxis right;
        plot(kmax_list, isi_curve, '-s', 'LineWidth', 1.5);
        ylabel('ISI energy (block n=2)');
        % Without an explicit ylim, MATLAB autoscales the ~1e-13 floating-point
        % ripple of a constant curve into a full-height axis.
        ylim([0, max(isi_curve) * 1.5 + eps]);
        xlabel('k_{max}'); grid on;
        xticks(kmax_list);
        title('ICI grows with Doppler spread; ISI is invariant to it');
        legend({'ICI ratio', 'ISI energy'}, 'Location', 'northwest');
    
        % ICI / ISI vs delay spread -- the sweep that moves ISI
        figure('Name', 'Sim-7 ICI/ISI vs l_{max}', 'Color', 'w');
        yyaxis left;
        plot(lmax_list, ici_curve_l, '-o', 'LineWidth', 1.5);
        ylabel('ICI ratio (block n=2)');
        ylim([0, max(ici_curve_l) * 1.2 + eps]);
        yyaxis right;
        plot(lmax_list, isi_curve_l, '-s', 'LineWidth', 1.5);
        ylabel('ISI energy (block n=2)');
        ylim([0, max(isi_curve_l) * 1.2 + eps]);
        xlabel('l_{max}'); grid on;
        xticks(lmax_list);
        title('ISI grows with delay spread');
        legend({'ICI ratio', 'ISI energy'}, 'Location', 'northwest');
    end
    
    end
    
    % -------------------------------------------------------------------------
    function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
    end