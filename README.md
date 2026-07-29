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
└── Stage2/             # Sim-5 ~ Sim-9 (+ Sim-10 ~ Sim-13): on-grid channel, 4-domain I/O
    ├── gen_channel_taps.m
    ├── gen_gs.m
    ├── gen_G.m
    ├── apply_channel_conv.m
    ├── get_block.m
    ├── build_channel_matrices.m
    ├── sim5_time_domain_convolution.m
    ├── sim6_block_wise_time_domain.m
    ├── sim7_time_frequency_domain.m
    ├── sim8_delay_time_domain.m
    ├── sim9_delay_doppler_domain.m
    └── run_stage2.m
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

## Reference

Yi Hong, Tharaj Thaj, Emanuele Viterbo,
*Delay-Doppler Communications: Principles and Applications*, Elsevier, 2022 — Chapter 4.
Equation numbers throughout this repository refer to this book.
