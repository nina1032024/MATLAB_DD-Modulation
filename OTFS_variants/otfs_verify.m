function ok = otfs_verify(g, ell, kap, Ls, LM, M, N, Lg, P, A, tag)
    %OTFS_VERIFY  Six-part self-verification confirming the OTFS three-domain matrices are correct
    %
    %  Simplest call (no arguments): just type at the command window
    %     >> otfs_verify
    %  This automatically runs three default scenarios: on-grid Lg=0 / on-grid Lg=2 / off-grid Lg=5.
    %
    %  Or call it from inside the main loop of otfs_variants_spy.m:
    %     otfs_verify(g, cfg(c).ell, cfg(c).kap, cfg(c).Ls, LM, M, N, Lg, P, A, cfg(c).tag);
    %
    %  This file carries its own gs_table / build_G / sinc_n and does not
    %  depend on otfs_variants_spy.m.
    %
    %  T1 algebraic identities   P is a permutation, A is unitary, Frobenius norm and support set are preserved across the three domains
    %  T2 time-domain assembly   does G*s match the book's element-wise formulas (4.64)/(4.82)/(4.91)/(4.109)
    %  T3 block structure        do K~_{m,l}, K_{m,l} match Table 4.2 (holds both on-grid and off-grid)
    %  T4 closed-form I/O        does H*x match the Table 4.3 phase table (on-grid only, Lg=0 only)
    %  T5 closed-form nnz        does the nonzero count match the hand-derived value (on-grid only, Lg=0 only)
    %  T6 grid test              is kappa*gamma_g an integer? If not, predicts no 1-1/N collapse will occur

      %% ===== no arguments: run the default scenarios then return =====
      if nargin == 0
        M = 8;  N = 6;  L = M*N;  LM = 7;
        g = [1.00, 0.75, 0.55];

        P = zeros(L,L);                                  % (4.33)
        for m = 0:M-1
          for n = 0:N-1
            P(n*M+m+1, m*N+n+1) = 1;
          end
        end
        FN = fft(eye(N))/sqrt(N);
        A  = kron(eye(M), FN);                           % the I_M kron F_N from (4.60)

        sc(1) = struct('ell',[0 1 2],       'kap',[0 1 2],        'Ls',0.5, 'Lg',0, ...
                       'tag','on-grid,  Lg=0  (gamma_g = 1)');
        sc(2) = struct('ell',[0 1 2],       'kap',[0 1 2],        'Ls',0.5, 'Lg',2, ...
                       'tag','on-grid,  Lg=2  (gamma_g = 1.25)');
        sc(3) = struct('ell',[0.3 1.7 2.4], 'kap',[0.3 1.7 -1.2], 'Ls',3,   'Lg',5, ...
                       'tag','off-grid, Lg=5  (gamma_g = 1.625)');

        allok = true;
        for c = 1:numel(sc)
          allok = otfs_verify(g, sc(c).ell, sc(c).kap, sc(c).Ls, LM, M, N, ...
                              sc(c).Lg, P, A, sc(c).tag) && allok;
        end
        fprintf('\n>>> OVERALL: %s\n', ternary(allok,'ALL SCENARIOS PASSED','SOMETHING FAILED'));
        if nargout > 0, ok = allok; end
        return
      end

      L        = M*N;
      variants = {'CP','ZP','RCP','RZP'};
      z        = exp(1j*2*pi/L);
      gamma_g  = 1 + Lg/M;
      Q        = (M+Lg)*N;
      gs       = gs_table(g, ell, kap, Ls, LM, Q, L);
      Lset     = find(any(abs(gs) > 1e-13, 2)).' - 1;    % delay taps that are actually nonzero
      lmax     = max(Lset);
      onGrid   = all(abs(ell-round(ell)) < 1e-12) && all(abs(kap-round(kap)) < 1e-12);

      Pi  = circshift(eye(N), -1, 2);                    % circular left shift of columns
      Dm  = diag(exp(-1j*2*pi*(0:N-1)/N));               % the phase diagonal matrix from (4.76)
      tol = 1e-9;  ok = true;

      fprintf('\n########## VERIFY: %s ##########\n', tag);
      fprintf('  M=%d  N=%d  Lg=%d  gamma_g=%.4g  l in {%s}  on-grid=%d\n', ...
              M, N, Lg, gamma_g, num2str(Lset), onGrid);

      %% ===== T1a/b: variant-independent algebraic identities =====
      ok = rep('T1a  P orthogonal',  norm(P.'*P - eye(L),'fro'), tol, ok);
      ok = rep('T1b  A unitary',     norm(A'*A  - eye(L),'fro'), tol, ok);

      rng(0);
      s = (randn(L,1) + 1j*randn(L,1))/sqrt(2);
      x = (randn(L,1) + 1j*randn(L,1))/sqrt(2);

      for j = 1:4
        v      = variants{j};
        perBlk = any(strcmp(v, {'CP','ZP'}));
        G      = build_G(v, gs, M, N, L, Lg);
        Ht     = P.' * G * P;
        H      = A * Ht * A';
        fprintf('  ---- %s-OTFS ----\n', v);

        %% ===== T1c/d: conservation across the three domains =====
        nG = norm(G,'fro');
        ok = rep('T1c  ||G||=||Ht||=||H||', ...
                 max(abs(nG-norm(Ht,'fro')), abs(nG-norm(H,'fro'))), tol*nG, ok);
        ok = rep('T1d  nnz(G)==nnz(Ht)', abs(nnzt(G)-nnzt(Ht)), 0.5, ok);

        %% ===== T2: time-domain assembly vs. the book's element-wise formulas =====
        r_ref = zeros(L,1);
        for n = 0:N-1
          for m = 0:M-1
            q = m + n*M;
            if perBlk, qc = m + n*(M+Lg); else, qc = q; end
            for l = Lset
              val = gs(l+1, qc+1);
              switch v
                case 'CP'                                            % (4.91)
                  r_ref(q+1) = r_ref(q+1) + val*s(n*M + mod(m-l,M) + 1);
                case 'ZP'                                            % (4.109)
                  if m >= l, r_ref(q+1) = r_ref(q+1) + val*s(n*M + (m-l) + 1); end
                case 'RCP'                                           % (4.82)
                  r_ref(q+1) = r_ref(q+1) + val*s(mod(q-l,L) + 1);
                case 'RZP'                                           % (4.64)
                  if q-l >= 0, r_ref(q+1) = r_ref(q+1) + val*s(q-l+1); end
              end
            end
          end
        end
        ok = rep('T2   G*s == book eqs', norm(G*s - r_ref), tol*norm(r_ref), ok);

        %% ===== T3: block structure vs. Table 4.2 =====
        e_ge = 0;  e_lt = 0;
        for m = 0:M-1
          for l = Lset
            Kt = Ht(m*N+(1:N), mod(m-l,M)*N+(1:N));
            K  = H (m*N+(1:N), mod(m-l,M)*N+(1:N));
            if m >= l
              e_ge = max([e_ge, resid_diag(Kt), resid_circ(K)]);      % diag / circ
            else
              switch v
                case 'CP'                                            % same as m>=l
                  e_lt = max([e_lt, resid_diag(Kt), resid_circ(K)]);
                case 'ZP'                                            % zero matrix
                  e_lt = max([e_lt, norm(Kt,'fro'), norm(K,'fro')]);
                case {'RCP','RZP'}                                   % diag*Pi / circ*D
                  e_lt = max([e_lt, resid_diag(Kt*Pi.'), resid_circ(K*Dm')]);
              end
            end
          end
        end
        ok = rep('T3a  m>=l : diag & circ', e_ge, tol, ok);
        ok = rep('T3b  m<l  : Table 4.2',   e_lt, tol, ok);

        %% ===== T4/T5: closed form only available on-grid with Lg=0 =====
        if onGrid && Lg == 0
          % ---- T4: the phase table from Table 4.3 ----
          X  = reshape(x,    N, M).';                                % x(m*N+n)=X[m,n]
          Y  = reshape(H*x,  N, M).';
          Yr = zeros(M,N);
          for i = 1:numel(g)
            li = round(ell(i));  ki = round(kap(i));
            for m = 0:M-1
              for n = 0:N-1
                Xv = X(mod(m-li,M)+1, mod(n-ki,N)+1);
                if m >= li
                  a = z^(ki*(m-li));
                else
                  switch v
                    case 'CP',  a = z^(ki*(m-li));
                    case 'ZP',  a = 0;
                    case 'RCP', a = z^(ki*mod(m-li,M)) * exp(-1j*2*pi*n/N);
                    case 'RZP', a = (N-1)/N * z^(ki*(m-li)) * exp(-1j*2*pi*mod(n-ki,N)/N);
                  end
                end
                Yr(m+1,n+1) = Yr(m+1,n+1) + g(i)*a*Xv;
              end
            end
          end
          e4 = norm(Y-Yr,'fro')/norm(Y,'fro');
          if strcmp(v,'RZP')
            fprintf('    %-22s %.3e  (4.122) is an O(1/N) approximation\n', 'T4   Table 4.3', e4);
          else
            ok = rep('T4   Table 4.3 exact', e4, 1e-10, ok);
          end

          % ---- T5: closed-form nnz ----
          nBlk = M*(lmax+1);                       % total number of (m,l) blocks for l in 0..lmax
          nLt  = lmax*(lmax+1)/2;                  % number of blocks among those with m<l
          nl   = arrayfun(@(l) sum(round(ell)==l), Lset);   % |K_l|
          Sset = sum(nl);                          % S = sum_l |K_l|
          switch v
            case 'CP',  pG = M*(lmax+1)*N;         pH = nBlk*Sset/numel(Lset)*N;
            case 'ZP',  pG = N*sum(M-(0:lmax));    pH = (nBlk-nLt)*Sset/numel(Lset)*N;
            case 'RCP', pG = L*(lmax+1);           pH = nBlk*Sset/numel(Lset)*N;
            case 'RZP', pG = sum(L-(0:lmax));      pH = (nBlk-nLt)*Sset/numel(Lset)*N + nLt*N^2;
          end
          ok = rep('T5a  nnz(G) closed form', abs(nnzt(G)-pG), 0.5, ok);
          ok = rep('T5b  nnz(H) closed form', abs(nnzt(H)-pH), 0.5, ok);
        end

        %% ===== T6: Doppler grid test =====
        kg   = kap*gamma_g;
        kgOK = all(abs(kg-round(kg)) < 1e-12);
        fprintf('    T6   kappa*gamma_g = [%s]  integer=%d  nnz(H)=%d\n', ...
                num2str(kg,'%.4g '), kgOK, nnzt(H));
        if ~kgOK
          fprintf('         -> off the Doppler grid: K_{m,l} stays full, no 1-1/N collapse\n');
        end
      end

      fprintf('  ==> %s\n', ternary(ok,'ALL CHECKS PASSED','SOME CHECKS FAILED'));
    end

    %% ---------------- helper functions ----------------
    function n = nnzt(X)
    % nnz with a relative threshold. Must use the same threshold as
    % otfs_variants_spy.m's mask, otherwise the FFT floating-point dust
    % (~1e-17) in H = A*Ht*A' would get counted as nonzero.
      mx = max(abs(X(:)));
      if mx == 0, n = 0; return; end
      n = nnz(abs(X) > 1e-9*mx);
    end

    function e = resid_diag(K)
      n = norm(K,'fro');  if n < 1e-14, e = 0; return; end
      e = norm(K - diag(diag(K)), 'fro') / n;
    end

    function e = resid_circ(K)
      n = norm(K,'fro');  if n < 1e-14, e = 0; return; end
      e = norm(K - circshift(K,[1 1]), 'fro') / n;
    end

    function ok = rep(name, err, tol, ok)
      if err <= tol
        fprintf('    %-22s PASS      (%.2e)\n', name, err);
      else
        fprintf('    %-22s **FAIL**  (%.2e > %.2e)\n', name, err, tol);
        ok = false;
      end
    end

    function s = ternary(c,a,b)
      if c, s = a; else, s = b; end
    end

    %% ---------------- the three subfunctions below are identical to otfs_variants_spy.m ----------------
    %  They are duplicated here so that otfs_verify.m can run standalone
    %  (MATLAB local functions don't cross files). If otfs_variants_spy.m's
    %  version is ever updated, remember to update this copy too.

    function gs = gs_table(g, ell, kap, Ls, LM, Q, NM)
    % discrete time-varying tap g^s[l,q], book Eq. (4.6)
    %   gs(l+1,q+1) = sum_i g_i * sinc(l-ell_i) * z^{kap_i (q-l)},  z = exp(j2pi/NM)
    % NM is always M*N (the Doppler normalization base; it does not change with the guard)
    % Q is the number of samples to evaluate; variants with a per-block guard need Q = (M+Lg)*N > NM
    % Ls is the sinc truncation half-width; on-grid, taking 0.5 degenerates to (4.8)
      gs = zeros(LM+1, Q);
      q  = 0:Q-1;
      for i = 1:numel(g)
        for l = 0:LM
          if abs(l - ell(i)) > Ls, continue; end
          w = sinc_n(l - ell(i));
          if abs(w) < 1e-13, continue; end
          gs(l+1,:) = gs(l+1,:) + g(i)*w*exp(1j*2*pi*kap(i)*(q-l)/NM);
        end
      end
    end

    function G = build_G(v, gs, M, N, L, Lg)
    % assemble the time-domain channel matrix G per variant, (4.38)/(4.83)/(4.93)/(4.109)
    %   q  : matrix row index with the guard removed, always m + n*M
    %   qc : physical sample time. CP/ZP insert the guard between blocks, so qc = m + n*(M+Lg)
      G = zeros(L,L);
      active = find(any(abs(gs) > 1e-13, 2)).' - 1;
      perBlockGuard = any(strcmp(v, {'CP','ZP'}));
      for n = 0:N-1
        for m = 0:M-1
          q  = m + n*M;
          if perBlockGuard, qc = m + n*(M+Lg); else, qc = q; end
          for l = active
            val = gs(l+1, qc+1);
            switch v
              case 'CP'
                cc = n*M + mod(m-l,M) + 1;
                G(q+1, cc) = G(q+1, cc) + val;
              case 'ZP'
                if m >= l
                  cc = n*M + (m-l) + 1;
                  G(q+1, cc) = G(q+1, cc) + val;
                end
              case 'RCP'
                cc = mod(q-l, L) + 1;
                G(q+1, cc) = G(q+1, cc) + val;
              case 'RZP'
                if q - l >= 0
                  G(q+1, q-l+1) = G(q+1, q-l+1) + val;
                end
            end
          end
        end
      end
    end

    function y = sinc_n(x)
    % normalized sinc: sin(pi x)/(pi x), equal to 1 at x=0 (avoids needing the Signal Processing Toolbox)
      y = ones(size(x));
      k = (x ~= 0);
      y(k) = sin(pi*x(k)) ./ (pi*x(k));
    end
