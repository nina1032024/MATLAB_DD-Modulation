# Chapter 4 OTFS 模擬 — 階段三：off-grid（分數延遲／分數都卜勒）

本資料夾實作 `chap4_list.md` 中 **階段三（Sim-14 ~ Sim-16）**：把階段二的 on-grid（整數
延遲、整數都卜勒 tap）假設放寬，改用分數（非整數）tap，驗證能量洩漏到相鄰 delay/Doppler
bin 的理論式，並分析通道稀疏度。

---

## 檔案總覽

### 輔助函式

| 檔案 | 說明 |
|---|---|
| `sinc_fn.m` | 自行實作的 `sinc(x)=sin(πx)/(πx)`，不依賴 Signal Processing Toolbox 的 `sinc` |
| `periodic_sinc.m` | 週期 sinc（Dirichlet 核）`𝒟_N(x)=sin(πx)/sin(πx/N)·e^{jπx(N-1)/N}`（(4.79)），含 `x` 為 `N` 的整數倍時的極限修正（`𝒟_N(0)=N`） |
| `gen_gs_frac.m` | 把 Stage 2 的 `gen_gs.m` 推廣到**分數延遲** `l_i`，用 `sinc(l-l_i)` 內插取代 `δ[l-l_i]`（(4.6)），整數 `l_i` 時完全退化回 Stage 2 的結果 |
| `effective_path_count.m` | 「有效路徑數 `S`」：累積能量達到門檻（預設 99%）所需的最少項數，Sim-15、Sim-16 共用 |

### 模擬項目

| 檔案 | 對應 Sim |
|---|---|
| `sim14_fractional_doppler_leakage.m` | **Sim-14**：分數都卜勒洩漏驗證，重現 Fig. 4.16 |
| `sim15_fractional_delay_leakage.m` | **Sim-15**：分數延遲洩漏驗證，`sinc` 內插模型 |
| `sim16_sparsity_analysis.m` | **Sim-16**：通道稀疏度 `S` 分析（`S` vs. `N`、`S` vs. 分數程度） |
| `run_stage3.m` | **總驅動程式**：依序執行 Sim-14~16，輸出誤差總表 |

> 沿用 Stage 1/2 的 `addpath` 慣例：每支程式在檔案開頭自動
> `addpath('../Stage1')` 與 `addpath('../Stage2')`，直接重用
> `otfs_params/gen_dd_symbols/gen_perm_matrix/otfs_modulate/otfs_demodulate/relerr`
> （Stage 1）以及 `gen_channel_taps/gen_gs/gen_G/apply_channel_conv/get_block/
> build_channel_matrices`（Stage 2），資料夾需維持
> `chap4_matlab/{Stage1,Stage2,Stage3}` 同一層結構。

---

## 如何執行

```matlab
cd chap4_matlab/Stage3

T = run_stage3();                       % M=8, N=6，驗證用尺寸，無圖
T = run_stage3(otfs_params(), true);    % 附上所有圖（Fig.4.16 重現、洩漏曲線、稀疏度曲線）
```

**執行結果（Sim-14、Sim-15 為 gated 檢定，皆通過；Sim-16 為觀察性分析）：**

```
                                                Test                                                  RelError      Tol      Result
    _____________________________________________________________________________________________    __________    _____    ________

    {'Sim-14 fractional Doppler leakage (4.79/4.80)'                                            }    3.2062e-15    1e-10    {'PASS'}
    {'Sim-15 fractional delay leakage (4.6)'                                                    }     1.855e-16    1e-10    {'PASS'}
    {'Sim-16 channel sparsity S analysis (informational, no numeric threshold in chap4_list.md)'}             0      Inf    {'PASS'}
```

---

## 各 Sim 的重點說明

### Sim-14　分數都卜勒的洩漏驗證

**理論推導**（`periodic_sinc.m` 的檔頭有完整推導，重點摘要如下）：單一路徑、整數延遲 `l`、
可為分數的都卜勒 tap `k_i`：

```
ν̃_{m,l}[n] = g_i · z^{k_i(m-l)} · e^{j2π k_i n/N}                         (由 4.73 展開)
```

對 `n` 做正向 DFT（(4.77)）得到一個等比級數，其封閉解正是週期 sinc / Dirichlet 核：

```
ν_{m,l}[k] = (1/N) · g_i · z^{k_i(m-l)} · 𝒟_N(k_i - k)                    (4.79)(4.80)
```

當 `k_i` 為整數時，`𝒟_N(k_i-k)` 在 `k=k_i` 處取極限值 `N`（`𝒟_N(0)=N`），其餘 `k` 全為 0，
恰好退化成 Stage 2 的「單一 Doppler bin」結果；`k_i` 為分數時則洩漏到所有 `N` 個 bin。

**做法**：直接由 `gen_gs.m`（Stage 2 的函式，對分數 `k_i` 本來就適用，因為 `z^{k_i(q-l)}`
本身沒有要求 `k_i` 是整數）取出 `ν̃_{m,l}[n]`，用數值 DFT 算出 `ν_{m,l}[k]`，與 (4.80) 解析式
比對；`κ=2`（整數）與 `κ=2.5`（分數）兩種情況都驗證，誤差 `~1e-15`。

**附圖**：重現 Fig. 4.16——上圖 `κ=2` 能量完全集中在一個 bin，下圖 `κ=2.5` 洩漏到所有
bin，並疊上 `𝒟_N(κ-k)/N` 的虛線包絡。

### Sim-15　分數延遲的洩漏驗證

**理論**：`g^s[l,q] = Σ_i g_i z^{k_i(q-l)} sinc(l-l_i)`（(4.6)），用 `sinc` 內插取代整數 tap
的 `δ[l-l_i]`，實作於 `gen_gs_frac.m`；當 `l_i` 為整數時 `sinc(整數)=0`（`sinc(0)=1` 除外），
**完全退化**回 Stage 2 的 `gen_gs.m`。

**做法與結果**：
1. **自我一致性**（gated，`<1e-10`）：`G_frac*s` vs. 直接摺積（沿用 Stage 2 的
   `apply_channel_conv.m`，不需修改就能處理分數延遲），誤差 `1.855e-16`。
2. **回歸測試**（gated）：強制 `l_i` 為整數，`gen_gs_frac.m` 應與 Stage 2 的 `gen_gs.m`
   完全一致——結果誤差 **精確為 0**，且延伸出去的其餘列（`sinc(非零整數)`）殘留能量也是 0，
   確認推廣沒有破壞 Stage 2 的正確性。
3. **截斷誤差曲線**（觀察性）：`capturedFraction(L_range)=Σ_{l=0}^{Lrange-1} sinc(l-l_i)^2`，
   驗證 `chap4_list.md` 所說「`l` 的範圍要擴大到 `l_max` 以外，否則會截斷洩漏能量」。
4. **有效路徑數 `S` vs. 分數程度**（觀察性）：分數偏移從 `0` 掃到 `0.5`，`S`（99% 能量門檻）
   從 `1` 增加到 `25`，確認 `S≥P`（此處 `P=1`）。

**除錯發現：截斷曲線會停在 `~0.95`，不是趨近 `1`**——`sinc` 基底的完備性等式
`Σ_{l=-∞}^{∞} sinc(l-l_i)^2 = 1` 是對**所有整數**（含負的）成立，但延遲 tap 依照 Sec. A.3
的定義必須是**因果**的（`l≥0`），所以只累加非負 `l` 天生就補不滿 `1`——這**不是**程式的
截斷誤差，而是物理上「負延遲不存在」造成的真實能量缺口。實測 `l_i=1.4` 時，`captured
fraction` 從 `Lrange=3` 的 `0.8741` 一路增加到 `Lrange=30` 的 `0.9495` 後幾乎不再成長，
缺口 `~0.05` 正好對應被排除在外的負延遲側 sinc 尾巴（`sinc(-2.4)^2, sinc(-3.4)^2,...`
數量級吻合）。程式中的註解與 `fprintf` 訊息已依此修正，避免誤導讀者以為只要 `Lrange`
夠大就能逼近 `1`。

### Sim-16　通道稀疏度 `S` 分析

`chap4_list.md` 對這項只給「觀察重點」，沒有數值準則，因此本程式做以下**明確的設計選擇**
（詳見 `sim16_sparsity_analysis.m` 檔頭註解）：

- **`S` vs. `N`**：固定 `M`，掃描 `N=6,12,24`，且**沿用同一組分數 tap 數值**（不是由物理
  車速重新推導），單純觀察「`S` 相對於 `NM` 的比例」如何隨頻框維度變化。實測：

  | `N` | `NM` | 平均 `S` | `S/NM` |
  |---|---|---|---|
  | 6 | 48 | 31.1 | 0.648 |
  | 12 | 96 | 53.27 | 0.555 |
  | 24 | 192 | 70.19 | 0.366 |

  `S` 本身持續增加，但 **`S/NM`（相對稀疏度）持續下降**——這正是支撐第 6 章低複雜度偵測
  演算法（複雜度隨 `S` 而非 `NM` 成長）的關鍵觀察。

- **`S` vs. 分數程度**：固定 `N=6`，把預設的 `l_i=[0,1,2], k_i=[0,1,-2]` 都加上同一個分數
  偏移（`0→0.5`），平均 `S` 從 `4.75` 增加到 `30.65`，與 Sim-15 用「單路徑 sinc 係數」量到
  的趨勢一致（互相印證：一個從解析係數算，一個從完整 `H` 矩陣算）。

---

## 設計上的重點提醒

1. **`l_i` 分數、`k_i` 整數/分數皆可，但 `l` 本身（TDL 的列索引）恆為非負整數**：分數只發生
   在「路徑的真實延遲/都卜勒值」上，`gen_gs_frac.m` 內部仍是對整數 `l`（`0…Lrange-1`）做
   `sinc` 加權求和，而不是對 `l` 本身做內插。
2. **`gen_G.m`、`apply_channel_conv.m`、`build_channel_matrices.m` 完全不需修改**：這三個
   Stage 2 函式都是照 `gs` 矩陣的實際列數（`size(gs,1)-1`）運作，所以不管 `gs` 是 Stage 2
   的整數版（`gen_gs.m`）還是 Stage 3 的分數版（`gen_gs_frac.m`），都能直接重用——這是
   Stage 2 設計時預留的彈性，Stage 3 不必重造輪子。
3. **`sinc_fn.m` 自行實作**：避免依賴 Signal Processing Toolbox 的 `sinc`，與 Stage 1
   `unitary_dft.m` 不依賴 `dftmtx` 的作法一致。
4. **`periodic_sinc.m` 的奇異點處理**：`𝒟_N(x)` 在 `x` 為 `N` 的整數倍時是 `0/0`，程式用
   L'Hôpital 極限直接給出精確值 `N`，而不是依賴浮點數的自然消去（後者在 `sin(π·整數)`
   有 `~1e-16` 的浮點殘差，雖然量級夠小不影響 `1e-10` 門檻，但用極限公式更嚴謹）。

---

## 下一步

階段三（off-grid）驗證完成後，依 `chap4_list.md` 的建議執行順序，接下來是「階段四（選做）」：
Sim-17（RZP/RCP/CP/ZP 四種 `G` 建構方式比較）與 Sim-18（on-grid + LMMSE 偵測、初步 BER）。
