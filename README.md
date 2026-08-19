# OTFS Chapter 4 — MATLAB Verification Suite (Sim-1 ~ Sim-9)

MATLAB implementation and numerical verification of the OTFS system model
described in Chapter 4 ("Delay-Doppler Communications: Principles and
Applications", Yi Hong, Tharaj Thaj, Emanuele Viterbo, Elsevier 2022),
Sections 4.1–4.4.

The goal is **not** BER simulation. It is a layer-by-layer check that every
matrix/vector identity in the book actually holds numerically — modulator,
demodulator, and all four channel-representation domains — before any
detection or performance study is attempted.

Full simulation plan and equation-to-code cross-reference:
[`chap4_list.md`](chap4_list.md).

---

## What this covers

| Stage | Sim items | Question being answered |
|---|---|---|
| **Stage 1** | Sim-1 – Sim-4 | Are the modulator and demodulator self-consistent, with no channel at all? |
| **Stage 2** | Sim-5 – Sim-9 | Do the time, time-frequency, delay-time, and delay-Doppler domain representations of an **on-grid** multipath channel agree with each other? |
| **Extra** | Sim-17-style, on/off-grid | Do the closed-form delay-Doppler relation (Table 4.3) and the four frame variants (CP/ZP/RCP/RZP) agree with an independently-built off-grid convolution ground truth, on-grid *and* off-grid? See [`two_path_verification/`](two_path_verification) and [`OTFS_variants/`](OTFS_variants) below. |

Everything is verified by comparing independently-written implementations
of the same equation and checking the relative error against a numerical
tolerance (`1e-10` – `1e-12`), not by visual inspection alone.

> Sim-10 and beyond (symbol-level theory validation, ideal vs. rectangular
> pulse, fractional delay/Doppler, sparsity analysis, …) already have
> working code in [`Stage2/`](Stage2) and [`Stage3/`](Stage3), but are
> outside the scope of this README revision.

---

## Repository layout

```
chap4_matlab/
├── chap4_list.md      # Full derivation + simulation spec (equations ↔ code ↔ pass criteria)
├── Stage1/             # Sim-1 ~ Sim-4: modulator/demodulator self-consistency
│   ├── otfs_params.m
│   ├── unitary_dft.m
│   ├── gen_dd_symbols.m
│   ├── gen_perm_matrix.m
│   ├── otfs_modulate.m
│   ├── otfs_demodulate.m
│   ├── sim1_perfect_reconstruction.m
│   ├── sim2_modulation_paths.m
│   ├── sim3_permutation_matrix.m
│   ├── sim4_parseval.m
│   └── run_stage1.m
├── Stage2/             # Sim-5 ~ Sim-9 (+ Sim-10 ~ Sim-13): on-grid channel, 4-domain I/O
│   ├── gen_channel_taps.m
│   ├── gen_gs.m
│   ├── gen_G.m
│   ├── apply_channel_conv.m
│   ├── get_block.m
│   ├── build_channel_matrices.m
│   ├── sim5_time_domain_convolution.m
│   ├── sim6_block_wise_time_domain.m
│   ├── sim7_time_frequency_domain.m
│   ├── sim8_delay_time_domain.m
│   ├── sim9_delay_doppler_domain.m
│   └── run_stage2.m
├── two_path_verification/   # Closed-form (Table 4.3) vs. off-grid convolution ground truth
│   ├── dirichlet_N.m          # Closed-form periodic Dirichlet kernel D_N(x)
│   ├── delay_kernel.m         # Pluggable delay kernel p(x): none/zoh/tri/sinc
│   ├── otfs_build_G.m         # Ground-truth G: CP/ZP/RCP/RZP, fractional delay/Doppler taps, selectable kernel
│   ├── otfs_dd_closed_form.m  # Y from the closed form: 'integer' (Table 4.3) or 'fractional'
│   └── otfs_verify.m          # Driver: Test 1, 2a-i/ii/iii, 2d, 3, 4; logs console output to output.txt
└── OTFS_variants/      # CP/ZP/RCP/RZP structural comparison (Table 4.2/4.3/4.4, §4.5)
    ├── otfs_verify.m          # T1-T6 self-verification of G / H̃ / H per variant
    └── otfs_variants_spy.m    # spy-style structure plots of G, H̃, H for all 4 variants
```

Stage 2 depends on Stage 1 (`addpath('../Stage1')` inside each script), so
both folders must stay siblings under `chap4_matlab/`.

---

## Requirements

- MATLAB (no toolboxes required — `unitary_dft.m` replaces `dftmtx`, and
  all DFTs are unitary by construction).
- Tested with `M = 8, N = 6` (small, structure-revealing size used in the
  book's figures) and larger sizes (e.g. `M = 32, N = 16`).

---

## Quick start

```matlab
cd chap4_matlab/Stage1
T1 = run_stage1();              % Sim-1 ~ Sim-4, table of relative errors
T1 = run_stage1(otfs_params(), true);   % also plot spy(P), Fig. 4.10, etc.

cd ../Stage2
T2 = run_stage2();              % Sim-5 ~ Sim-9 (and Sim-10~13), table of relative errors
T2 = run_stage2(otfs_params(), true);   % also plot |G|, |Ȟ|, |H̃|, |H|, ICI/ISI curves
```

Each `run_stageX` returns a table (or struct array as fallback) that can be
pasted directly into a report:

```
                            Test                             RelError      Tol      Result
    ____________________________________________________    __________    _____    ________

    {'Sim-1 IDZT/DZT perfect reconstruction'           }    2.2061e-15    1e-12    {'PASS'}
    {'Sim-2 modulation / demodulation path equivalence'}    1.0234e-15    1e-12    {'PASS'}
    {'Sim-3 permutation matrix P'                      }    1.2399e-16    1e-12    {'PASS'}
    {'Sim-4 Parseval / energy conservation'            }    6.9485e-15    1e-12    {'PASS'}
```

```
    {'Sim-5 time-domain channel: G*s vs. direct convolution'}    1.0203e-16    1e-12    {'PASS'}
    {'Sim-6 block-wise time-domain relation (4.39)'         }    2.4004e-17    1e-12    {'PASS'}
    {'Sim-7 time-frequency domain I/O relation (4.4.2)'     }    9.7077e-16    1e-10    {'PASS'}
    {'Sim-8 delay-time domain I/O relation (4.4.3)'         }     3.278e-17    1e-10    {'PASS'}
    {'Sim-9 delay-Doppler domain I/O relation (4.4.4)'      }    1.3173e-15    1e-10    {'PASS'}
```

All items pass with relative errors several orders of magnitude below the
required tolerance.

---

## Stage 1 — Modulator/demodulator self-consistency (Sim-1 ~ Sim-4)

Before any channel is introduced, the transmit/receive chain has to agree
with itself across every equivalent formulation given in the book.

| Sim | Checks | Theory |
|---|---|---|
| **Sim-1** | Perfect reconstruction: with no channel and no noise, `Y == X` | 16 combinations of 4 modulator implementations (`isfft`/`idzt`/`elem`/`vec`) × 4 demodulator implementations (`wigner`/`dzt`/`elem`/`vec`) |
| **Sim-2** | ISFFT+Heisenberg path (Fig. 4.1) ≡ IDZT path (Fig. 4.3) | Eq. (4.16)–(4.30) |
| **Sim-3** | Row-column interleaver `P`: two independently-built constructions match exactly, and `PᵀP = I_NM` | Eq. (4.32)–(4.34) |
| **Sim-4** | Parseval / energy conservation across DD, time, TF and delay-time domains | Unitarity of `F_M`, `F_N`, `P` |

Key implementation notes:
- All DFTs are **unitary** (`F_K = e^{-j2πab/K}/√K}`); mixing in raw
  `fft`/`ifft` normalization is the single most common source of failure.
- Two *different* vectorizations coexist — `s` (delay index fast) and `x`
  (Doppler index fast). The permutation matrix `P` exists solely to convert
  between them.
- Every modulator/demodulator output is a normalized `NM×1` column vector
  (the book's reference code is inconsistent about row vs. column vectors).

See [`Stage1/README.md`](Stage1/README.md) for full details.

---

## Stage 2 — On-grid channel, four-domain I/O relations (Sim-5 ~ Sim-9)

With the modulator/demodulator verified, an on-grid (integer delay,
integer Doppler tap) multipath channel is added, and the same physical
signal is checked to produce consistent results in all four domains the
book derives it in.

| Sim | Domain | Checks | Theory |
|---|---|---|---|
| **Sim-5** | Time | `r = G·s` (matrix) vs. direct tap-delay-line convolution | Eq. (4.8), (4.36), (4.38) |
| **Sim-6** | Time (block-wise) | Only the diagonal block and first sub-diagonal block of `G` are non-zero; block-wise reconstruction matches `G·s` | Eq. (4.39) |
| **Sim-7** | Time-frequency | `y̌ = Ȟx̌`, `Ȟ = (I_N⊗F_M)G(I_N⊗F_M^†)`; ICI/ISI energy decomposition vs. Doppler spread | Eq. (4.40)–(4.49) |
| **Sim-8** | Delay-time | `ỹ = H̃x̃`, `H̃ = PᵀGP`; structure of each sub-block `K̃_{m,l}` | Eq. (4.52)–(4.56), (4.73) |
| **Sim-9** | Delay-Doppler | `y = Hx`, `H = (I_M⊗F_N)H̃(I_M⊗F_N^†)`; circulant structure of each sub-block `K_{m,l}` | Eq. (4.58)–(4.61), (4.77) |

**Notable finding** (documented in full in
[`Stage2/README.md`](Stage2/README.md)): the book's short-form summary of
Sim-8/Sim-9 states that `K̃_{m,l}` is always diagonal and `K_{m,l}` is
always circulant. Numerically, this only holds when `m ≥ l`. When the
delay tap forces a wrap into the *previous* time block (`m < l`), the true
structure is a shifted diagonal / a product of a circulant matrix and a
diagonal phase matrix — which is **not**, in general, circulant. This
edge case is consistent with, and predicted by, the frame-boundary term
already present in the book's symbol-level equations (4.120)–(4.121); it
was simply omitted from the short verbal description of Sim-8/Sim-9. The
end-to-end I/O relations (`ỹ = H̃x̃`, `y = Hx`) are unaffected and verified
to `~1e-15`.

Two channel boundary conventions are supported side by side via `gen_G.m`:
- **RZP** (Reduced Zero Padding): `q − l < 0` is zero-padded — matches the
  book's "single isolated frame" assumption used for Sim-5 ~ Sim-9.
- **RCP** (Reduced Cyclic Prefix): `q − l < 0` wraps around modulo `NM` —
  used by the symbol-level theory generator (Sim-10 onward).

See [`Stage2/README.md`](Stage2/README.md) for full details, including the
corrected `H̃ = PᵀGP` convention (the book's Appendix C reference code has
`H̃ = P·G·Pᵀ`, which is inconsistent with its own `x̃ = Pᵀs` definition).

---

## Two-path verification — closed-form vs. convolution ground truth (`two_path_verification/`)

Standalone scripts (no dependency on Stage 1/2) that check the book's symbol-level
closed-form delay-Doppler relation (Table 4.3, and its fractional generalization)
against an independently-built ground-truth path: modulate → apply the off-grid
channel matrix `G` → demodulate.

```
Path A (ground truth):  X --IDZT--> s --G (off-grid, otfs_build_G)--> r --DZT--> Y_A
Path B (model):         X --closed-form H_DD (Table 4.3) or matrix identity (4.60)--> Y_B
```

| File | Role |
|---|---|
| `dirichlet_N.m` | Closed-form periodic Dirichlet kernel `D_N(x) = Σ_{p=0}^{N-1} e^{j2πxp/N}`, with the on-grid limit handled as `N·1{x ≡ 0 (mod N)}` (not a bare `x==0` test — negative multiples of `N` matter too) |
| `delay_kernel.m` | Pluggable delay-domain pulse kernel `p(x)`: `'none'` (integer delay only, no interpolation), `'zoh'` (rect Tx + ideal sampler), `'tri'` (rect Tx + matched rect Rx, 2 taps), `'sinc'` (ideal band-limited Tx/Rx pair, the original hard-coded default) |
| `otfs_build_G.m` | Builds the time-domain channel matrix `G` (ground truth) for CP/ZP/RCP/RZP, off-grid delay/Doppler, kernel selectable via `p.kernel` (default `'tri'`) |
| `otfs_dd_closed_form.m` | Computes `Y` directly from the closed-form relation: `'integer'` mode (Table 4.3 / (4.118)–(4.122)) or `'fractional'` mode (periodic-Dirichlet expansion, (4.78)/(4.105), CP/ZP only), same `p.kernel` selection |
| `otfs_verify.m` | Driver script running the tests below; also mirrors everything printed to the console into `output.txt` (own `fopen`/`fclose` + `try/catch`, not `diary()` — see note below) |

| Test | Checks | Theory |
|---|---|---|
| **Test 1** | Integer delay/Doppler taps, all four variants: convolution path vs. matrix identity (4.60) vs. symbol-wise Table 4.3 | Table 4.3, Eq. (4.60) |
| **Test 2a-i** | Doppler-only sweep, `ℓ` locked to integers (`KERNEL='none'`, so no delay interpolation is exercised at all): NMSE of the integer-tap approximation vs. the fractional closed form vs. a closed-form leakage prediction `2 − 2·Re{D_N(ε)}/N`, plus a single-path diagnostic curve | Eq. (4.105) |
| **Test 2a-ii** | Delay-only sweep, Doppler locked to integers: 3×3 NMSE matrix crossing 3 ground-truth pulses (`sinc`/`tri`/`zoh`) against 3 model kernels — only the matched (diagonal) kernel should be exact | Eq. (4.6) generalized |
| **Test 2a-ii(b)** | Guard-length sweep (`L_g` = 6…30) for the `sinc` kernel: the error floor (~ −12 dB) does **not** improve with more guard length | — |
| **Test 2a-iii** | Phase convention: `z^{κ(m−ℓ_i)}` (true delay) vs. the book's `z^{κ(m−l)}` (tap index), kernel fixed to `'tri'` | Eq. (4.105)/(4.118) |
| **Test 2d** | LS channel estimation from `T=8` pilot frames, Monte Carlo (200 trials): integer model (`P` unknowns/delay-bin) vs. fractional model (`P·N` unknowns/delay-bin), MSE vs. pilot SNR, −10…40 dB. This is the one test where the "cruder" integer model legitimately wins (at low SNR / small offset) — a bias-variance trade-off, not a ground-truth mismatch | — |
| **Test 3** | Time-frequency domain: exact identity `y̌ = Ȟx̌` (always holds) vs. the ideal-pulse single-tap approximation (biorthogonality loss) | Eq. (4.13)–(4.14), (4.42) |
| **Test 4** | Doppler-spread vector `ν_{m,l}[k]` at integer vs. fractional offset, reproducing Fig. 4.16 | Eq. (4.79)–(4.80) |

**Notable findings:**

- **Ground-truth phase bug (the significant one).** `otfs_build_G.m` computed the
  Doppler-tilt phase as `z^{κ_i(q−l)}`, using the discretized tap index `l` instead
  of the true (possibly fractional) delay `ℓ_i`. Physically the phase comes
  straight from `e^{j2πν_i(qT_s − τ_i)}` — it does not depend on which tap `l` the
  delay kernel happens to be evaluated at. For integer delay (`l ≡ ℓ_i`) this is
  invisible, which is why it went undetected through Test 1 and the original
  Test 2 (both integer-delay-only). For fractional delay it is a ~30 dB error —
  and because it lived in the *ground truth* `G`, not just a model, it silently
  made fractional-delay model comparisons meaningless until fixed (now
  `z^{κ_i(q−ℓ_i)}`, matching the already-corrected phase in
  `otfs_dd_closed_form.m`'s `'fractional'` mode).
- **Kernel tap-range off-by-one.** An earlier hand-derived tap-range table used
  `l = floor(ℓ_i)` for the `'zoh'` kernel; since `zoh`'s support is `x ∈ [0,1)`,
  the correct tap is `l = ceil(ℓ_i)`, and the floor version silently zeroed the
  channel out (`G ≡ 0`, surfacing as `NaN`/`Inf` NMSE). Fixed by no longer
  hand-deriving floor/ceil offsets per kernel: `kernel_taps()` now deliberately
  scans a couple of taps wider than the exact support and lets the kernel
  evaluation zero out anything outside it.
- **Declared pulse shape vs. implemented kernel.** The system model declares
  rectangular pulse shaping, but the delay-domain interpolation kernel was
  hard-coded to `sinc`. `p.kernel` now parameterizes this explicitly (default
  changed to `'tri'`, matching rect Tx + matched rect Rx).
- **`sinc` kernel error floor under CP/ZP.** `sinc` is non-causal (infinite,
  two-sided support), which the block-circular CP/ZP structure cannot represent
  exactly; this leaves a ~ −12 dB error floor that does **not** improve with
  more guard length `L_g` (Test 2a-ii(b)) — a structural limitation of the
  `sinc` kernel under a block-circular frame, not a bug.


---

## OTFS frame variants — CP/ZP/RCP/RZP structural comparison (`OTFS_variants/`)

Standalone scripts (no dependency on Stage 1/2) extending the on-grid/off-grid
analysis to all four frame-boundary conventions defined in §4.5, reproducing
Table 4.2 / 4.3 / 4.4 and Fig. 4.8, 4.11, 4.12, 4.14, 4.15, 4.17–4.19.

| File | Role |
|---|---|
| `otfs_verify.m` | Six-part self-verification (T1–T6, see below) of `G`, `H̃ = PᵀGP`, `H = (I_M⊗F_N)H̃(I_M⊗F_N^†)` for all four variants, run on three built-in scenarios (on-grid `Lg=0`, on-grid `Lg=2`, off-grid `Lg=5`) |
| `otfs_variants_spy.m` | Renders the `spy`-style support-set (or dB-magnitude) structure of `G`, `H̃`, `H` side by side for CP/ZP/RCP/RZP, plus a plot of which entries are lost when the Doppler tap is rounded to an integer |

| Test | Checks | Theory |
|---|---|---|
| **T1** | Variant-independent algebraic identities: `P` is a permutation (`PᵀP=I`), `A=I_M⊗F_N` is unitary, Frobenius norm and nnz are preserved across `G → H̃ → H` | — |
| **T2** | `G·s` matches the book's element-wise assembly for each variant | Eq. (4.64) RZP, (4.82) RCP, (4.91) CP, (4.109) ZP |
| **T3** | Block structure of `K̃_{m,l}` / `K_{m,l}` against Table 4.2, both for `m ≥ l` and the `m < l` corner case | Table 4.2 |
| **T4** | `H·x` matches the Table 4.3 closed-form phase relation (on-grid, `Lg=0` only) | Table 4.3 / (4.118)–(4.122) |
| **T5** | Closed-form nnz counts for `G` and `H` match hand-derived formulas (on-grid, `Lg=0` only) | — |
| **T6** | Grid test: is `κ·γ_g` an integer? Predicts whether `K_{m,l}` collapses to the `1−1/N` sparse pattern or stays full | Eq. (4.107) |

**Notable finding:** rev. 2 of `otfs_variants_spy.m` corrects an earlier assumption
that `γ_g = 1` always holds. CP/ZP sample `g^s` at `q = m + n(M+Lg)` per
(4.93)/(4.109), while RCP/RZP sample at `q = m + nM` per (4.38)/(4.83) — so for
CP/ZP with `Lg > 0`, the matrix row index and the physical sampling instant
diverge. The correct grid condition from (4.107) is `κ·γ_g ∈ ℤ`, not `κ ∈ ℤ`:
e.g. with `Lg=2, M=8`, `γ_g=1.25`, so an integer-Doppler physical channel is
still off-grid for CP/ZP. A `USE_GAMMA` switch reproduces the earlier
(incorrect) `γ_g=1` behavior for comparison.

---

## Reference

Yi Hong, Tharaj Thaj, Emanuele Viterbo,
*Delay-Doppler Communications: Principles and Applications*, Elsevier, 2022 — Chapter 4.
Equation numbers throughout this repository refer to this book.
