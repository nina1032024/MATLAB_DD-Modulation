# Chapter 4 OTFS 模擬 — 階段二：on-grid 通道，四個 domain 的 I/O 關係驗證

本資料夾實作 `chap4_list.md` 中 **階段二（Sim-5 ~ Sim-13）**：在階段一確認 modulator/demodulator
沒問題之後，加入 on-grid（整數延遲、整數都卜勒 tap）的多路徑通道，逐一驗證時域、時頻域、
延遲-時間域、延遲-都卜勒域四個 domain 的 I/O 關係，最後用符號級理論式做正式驗收。

> **補充內容**：`appendixC_extra/` 子資料夾額外驗證了 Appendix C 程式碼中先前 Stage1~3
> 都沒實際跑過的部分——3GPP 標準通道（EPA/EVA/ETU）、移動速度→都卜勒參數，以及附錄 E.4
> 要求的「r 的 4 種算法 × Y 的 3 種算法」整合誤差表。詳見該資料夾的 README.md。

---

## 檔案總覽

### 通道相關輔助函式

| 檔案 | 說明 |
|---|---|
| `gen_channel_taps.m` | 產生 on-grid 通道 `(g_i, l_i, k_i)`，預設 `l_i=[0 1 2], k_i=[0 1 -2]`，等功率、隨機相位 |
| `gen_gs.m` | 由 `(g_i,l_i,k_i)` 依 (4.8) 建立時域擴散函數 `g^s[l,q]`（`(l_max+1)×NM` 矩陣） |
| `gen_G.m` | 由 `g^s[l,q]` 依 (4.38) 建立時域通道矩陣 `G`（`NM×NM`），支援 `'RZP'`／`'RCP'` 兩種邊界處理 |
| `apply_channel_conv.m` | 獨立於 `G` 的純 for-loop 時變摺積，與 `gen_G.m` 互相驗證用（Sim-5） |
| `get_block.m` | 通用的區塊矩陣切割函式（`G/Ȟ` 用 `M`-block；`H̃/H` 用 `N`-block） |
| `build_channel_matrices.m` | 由 `G` 建立 `Ȟ`（時頻域）、`H̃`（延遲-時間域）、`H`（延遲-都卜勒域） |
| `theory_dd_relation.m` | Sim-10 用的符號級理論產生器，實作 (4.120)(4.121)（矩形脈波、RCP） |
| `ideal_pulse_relation.m` | Sim-11 用的理想脈波理論式，實作 (4.15)（純 2D 循環摺積） |
| `add_awgn.m` | 產生複數 AWGN，`sigma_w^2` 可由 `SNR=E_s/N_0` 換算 |

### 模擬項目

| 檔案 | 對應 Sim |
|---|---|
| `sim5_time_domain_convolution.m` | **Sim-5**：矩陣 `G*s` vs. 直接摺積 |
| `sim6_block_wise_time_domain.m` | **Sim-6**：Block-wise 時域關係 (4.39) |
| `sim7_time_frequency_domain.m` | **Sim-7**：時頻域 I/O 關係、ICI/ISI 能量分析 |
| `sim8_delay_time_domain.m` | **Sim-8**：延遲-時間域 I/O 關係、`K̃_{m,l}` 結構 |
| `sim9_delay_doppler_domain.m` | **Sim-9**：延遲-都卜勒域 I/O 關係、`K_{m,l}` circulant 性 |
| `sim10_symbol_level_theory.m` | **Sim-10**：符號級理論式 vs. 模擬（最終驗收） |
| `sim11_ideal_vs_rectangular_pulse.m` | **Sim-11**：理想脈波 vs. 矩形脈波比較 |
| `sim12_single_path_shift_test.m` | **Sim-12**：單路徑位移測試 |
| `sim13_noise_statistics.m` | **Sim-13**：雜訊統計一致性 |
| `run_stage2.m` | **總驅動程式**：依序執行 Sim-5~13，輸出誤差總表 |

> 各檔案內部用 `addpath(fullfile(thisDir,'..','Stage1'))` 自動載入階段一的
> `otfs_params/gen_dd_symbols/gen_perm_matrix/otfs_modulate/otfs_demodulate/relerr`，
> 所以 Stage1、Stage2 資料夾需維持同一層（`chap4_matlab/Stage1`, `chap4_matlab/Stage2`）。

---

## 如何執行

```matlab
cd chap4_matlab/Stage2

T = run_stage2();                       % M=8, N=6，驗證用尺寸，無圖
T = run_stage2(otfs_params(), true);    % 附上所有圖（|G|,|Ȟ|,|H̃|,|H|,ICI/ISI 曲線等）
```

**最終結果（9 項全數通過）：**

```
                                           Test                                            RelError      Tol      Result
    __________________________________________________________________________________    __________    _____    ________

    {'Sim-5 time-domain channel: G*s vs. direct convolution'                         }    1.0203e-16    1e-12    {'PASS'}
    {'Sim-6 block-wise time-domain relation (4.39)'                                  }    2.4004e-17    1e-12    {'PASS'}
    {'Sim-7 time-frequency domain I/O relation (4.4.2)'                              }    9.7077e-16    1e-10    {'PASS'}
    {'Sim-8 delay-time domain I/O relation (4.4.3)'                                  }     3.278e-17    1e-10    {'PASS'}
    {'Sim-9 delay-Doppler domain I/O relation (4.4.4)'                               }    1.3173e-15    1e-10    {'PASS'}
    {'Sim-10 symbol-level DD-domain theory vs. simulation (final Stage-2 validation)'}    1.6994e-15    1e-10    {'PASS'}
    {'Sim-11 ideal vs. rectangular pulse (4.2)'                                      }             0    1e-10    {'PASS'}
    {'Sim-12 single-path delta shift test'                                           }    2.4488e-15    1e-10    {'PASS'}
    {'Sim-13 noise statistics across domains'                                        }     0.0089794     0.01    {'PASS'}
```

---

## 各 Sim 的重點說明

### Sim-5　時域通道兩種實作一致（4.4.1）
矩陣形式 `r=G*s` vs. 純 for-loop 的時變摺積 `r[q]=Σ_l g^s[l,q]s[q-l]`，兩者互相獨立實作
（`gen_G.m` 與 `apply_channel_conv.m` 沒有共用程式碼），比對到 `1e-16` 量級。此處與後續
Sim-6~9 皆使用 **RZP**（Reduced Zero Padding：`q-l<0` 直接補零，對應 4.4 節「單一 frame
孤立分析」的假設）。附圖重現 Fig. 4.8 的三條次對角線結構（`l_max+1=3` 條）。

### Sim-6　Block-wise 時域關係（4.39）
把 `G` 切成 `N` 個 `M×M` 區塊（區塊索引＝時間槽 `n`），驗證只有對角塊 `G_{n,0}` 與第一條
次對角塊 `G_{n,1}` 非零（因 `l_max<M`），並用 `r_n=G_{n,0}s_n+G_{n,1}s_{n-1}` 逐塊重建 `r`，
與整體 `G*s` 比對。

### Sim-7　時頻域 I/O 關係（4.4.2）
驗證 `y̌=Ȟx̌`，並用同一組 `G_{n,0}/G_{n,1}` 記號算出 `Ȟ_{n,0}=F_M G_{n,0}F_M^†`（自身塊，其
非對角能量即 ICI）與 `Ȟ_{n,1}=F_M G_{n,1}F_M^†`（前一塊，其總能量即 ISI），並掃描
`k_max=0,1,2` 觀察 ICI/ISI 隨都卜勒展寬增加的趨勢。

### Sim-8　延遲-時間域 I/O 關係（4.4.3）
驗證 `ỹ=H̃x̃`（`H̃=P^T G P`，見下方「修正」說明），並檢查每個子塊 `K̃_{m,l}` 的完整結構
（詳見下方「除錯發現」）。

### Sim-9　延遲-都卜勒域 I/O 關係（4.4.4）
驗證 `y=Hx`，並檢查每個子塊 `K_{m,l}` 是否為 circulant（詳見下方「除錯發現」）。

### Sim-10　符號級 DD 域理論式驗證（最終驗收）
純 for-loop 理論產生器（`theory_dd_relation.m`，實作 (4.120)(4.121)，**使用 RCP** 變體，對應
`chap4_list.md` Sec. E.1 code 16 的對照）直接由 `X` 與 `(g_i,l_i,k_i)` 算出 `Y_theory`，與
「`X→s→r=G_RCP·s→Y_sim`」的完整模擬路徑比對，誤差 `~1e-15`。

### Sim-11　理想脈波 vs. 矩形脈波比較（4.2 節）
比較 `ideal_pulse_relation.m`（(4.15) 純 2D 循環摺積，無相位）與矩形脈波的 `theory_dd_relation.m`。
`chap4_list.md` 本身沒有給這項嚴格的數值門檻，只有觀察重點，因此本程式回報三種情況（詳見下方）。

### Sim-12　單路徑位移測試
`P=1,g_1=1`，`X` 只放一個 delta，分別測純延遲、純都卜勒、兩者皆有三種情況，同時用理論式與
完整通道模擬（`G` 走 RCP）兩條獨立路徑驗證 peak 位置、大小、相位、以及其餘位置能量是否為 0。

### Sim-13　雜訊統計一致性
產生 `1e5` 組複數高斯雜訊，推過 `w→w̌(TF)→w̃(delay-time)→z(DD)` 四個么正變換，檢查每個
domain 的逐點變異數（門檻 `<1%`）與非對角相關（統計量，隨樣本數收斂，門檻 `<5%`）。

---

## 除錯發現：`chap4_list.md` 摘要與實際結構不完全一致的三處

這正是 `chap4_list.md` Sec. H.9 提到的精神：「模擬是用來檢查自己的理解對不對」。
第一次跑 `run_stage2()` 時，Sim-8、Sim-9、Sim-11 都 FAIL，追查後發現問題出在
**「`m<l_i`（延遲 tap 超出目前 delay bin）時，會有一個跨 frame-boundary 的邊界修正」**
這件事，在清單的精簡描述中被省略了。以下記錄問題與修正，供撰寫報告時參考。

### 1. Sim-8：`K̃_{m,l}` 並非對所有 `m` 都是「對角矩陣」

- **原始假設（依 chap4_list.md 字面描述）**：`K̃_{m,l}` 一律是對角矩陣，對角元素
  `ν̃_{m,l}[n] = g^s[l, m+nM]`（(4.73)）。
- **第一次執行結果**：off-diagonal 最大 Frobenius norm 高達 `1.291`（遠超過 `1e-12` 門檻），
  對角線比對誤差 `nu_err_max ≈ 1.0`（完全不對）。
- **追查過程**：把 `H̃=P^T G P` 的元素展開，`H̃` 在 (區塊 `m`, 區塊 `m'`) 的 `(n,n')` 元素
  剛好等於 `G[m+nM, m'+n'M]`。當 `m≥l`（沒有跨越 block 邊界）時，`n'=n` 對角成立；但當
  `m<l` 時，`q-l` 會落到「前一個時間 block」（`n'=n-1`），使得非零元素其實出現在
  **偏移一格的次對角線**（`(n,n-1)`），而非主對角線；在 RZP 模式下第一列（`n=0`）甚至
  完全為 0（因為沒有「block -1」可借用）。這與 Sim-10 已驗證過的 (4.121) 中 `m<l_i` 分支
  （借用前一個時間 block 的符號）是同一個物理效應。
- **修正**：`sim8_delay_time_domain.m` 改成針對每個 `(m,l)` 建立「完整的解析預期矩陣」
  （`m≥l` 時填對角線，`m<l` 時填偏移對角線，RZP 邊界填 0），再與實際的 `K̃_{m,l}` 比對，
  而不是單純檢查「非對角是否為 0」。修正後全部通過（`~1e-17`）。

### 2. Sim-9：`K_{m,l}` 並非對所有 `m` 都是 circulant

- **原始假設**：`K_{m,l}=F_N K̃_{m,l} F_N^†` 一律是 circulant（因為「對角矩陣經么正 DFT
  相似變換 → circulant」）。
- **問題根源**：這個性質只有在 `K̃_{m,l}` 是**真正的對角矩陣**時才成立（即 Sim-8 發現的
  `m≥l` 情況）。當 `m<l` 時 `K̃_{m,l}` 其實是「循環位移矩陣 × 對角矩陣」，可以證明
  `F_N·(shift·diag)·F_N^† = Λ_shift·(F_N·diag·F_N^†)`，也就是「對角矩陣 × circulant 矩陣」的
  乘積，一般而言**不是** circulant（除非 `Λ_shift` 是常數，但都卜勒相移對應的 shift 一般不是）。
- **修正**：`sim9_delay_doppler_domain.m` 把 circulant 檢查拆成兩組：`m≥l` 的「安全」區塊
  （嚴格門檻 `<1e-10`，全部通過）與 `m<l` 的「邊界」區塊（僅回報數值，不設門檻，因為理論上
  就是不會是 circulant）。I/O 關係本身（`y=Hx`）完全不受影響，一路都在 `1e-15` 量級。

### 3. Sim-11：「`k_i=0` 時兩者必然相同」的說法不夠精確

- **原始假設**：`chap4_list.md` 寫「當所有 `k_i=0` 時兩者應完全相同 → 可作為 sanity check」。
- **問題**：矩形脈波公式 (4.121) 在 `m<l_i` 分支多出的 `e^{-j2πn/N}` 邊界相位，其來源是
  **frame 邊界的時間位移本身**（借用前一個 time block 的符號後，經過 DZT 解調產生的線性相位），
  與都卜勒 `k_i` 完全無關；即使 `k_i=0`，只要 `l_max>0`（延遲 tap 造成跨 block 借用），這個
  相位依然存在，`Y_rect` 與 `Y_ideal` 仍會有非零差異。
- **驗證**：程式內第一次的「case 3」用 `l_i=0, k_i=2`（想測試「沒有延遲跨界」），但因為
  `k_i≠0`，矩形脈波公式仍保留 `z^{k_i(m-l_i)}` 相位，與理想脈波（無相位）本來就不該相同，
  這其實是設計錯誤，不是理論錯誤。
- **修正**：把「必然完全相同」的情境改為 `l_i=0 且 k_i=0`（真正的平凡情況：沒有延遲、沒有
  都卜勒，通道退化成純量增益），這時兩式都退化成 `Y=g_1·X`，相對誤差為 **精確的 0**。
  一般的「`k_i=0` 但 `l_max>0`」情況（case 2）改標記為「僅供觀察」，不再視為必須通過的檢定。

**共同結論**：三個問題的根源相同——矩形脈波（非雙正交）在 `m<l_i`（延遲 tap 造成需要借用
前一個 OTFS block 的符號）時，會產生一個額外的邊界相位/位移，這件事在 (4.120)(4.121) 的公式
本身已經正確描述（Sim-10 對這個完整公式的驗證從頭到尾都通過），但 `chap4_list.md` 對 Sim-8/9/11
的**文字摘要**把這個邊界情形省略了，只描述了「主要／安全」的情況。程式中已把正確的完整結構
寫清楚並加上詳細註解，可直接作為報告中「發現並修正教材摘要與正文推導不一致之處」的例子
（對應 `chap4_list.md` Sec. H 第 9 點的精神）。

---

## 設計上的重點提醒

1. **RZP vs. RCP**：Sim-5~9 使用 **RZP**（對應 4.4 節「單一 frame 孤立分析」假設，`gen_G.m`
   的 `q-l<0` 直接補零）；Sim-10~12 使用 **RCP**（`q-l<0` 時以 `mod NM` 環繞），因為
   (4.118)-(4.122) 的簡潔公式正是在 RCP 假設下推導出來的（對應 `chap4_list.md` Sec. E.1
   code 16 → Sim-10 的對照）。兩種變體都由同一份 `gen_G.m`/`apply_channel_conv.m` 透過
   `variant` 參數切換，不需要重複程式碼。
2. **`k_i` 一律用「原始（可為負）整數」而非 `mod N`**：`z=e^{j2π/(NM)}` 的週期是 `NM`，不是
   `N`，所以 `g^s[l,q]`、`G`、`Ȟ`、`H̃`、`H` 全部都必須用**帶正負號的原始 `k_i`**；只有在
   對 DD 符號矩陣 `X` 做索引（例如 `X[...,[n-k_i]_N]`）時才需要 `mod N`。`gen_channel_taps.m`
   的註解特別強調這點。
3. **`H̃ = P^T G P`，不是 `P G P^T`**：修正 `chap4_list.md` Sec. E.3.2 指出的 Appendix-C code 9
   錯誤（書上參考碼寫成 `H_tilde=P*G*P.'`，但由 `x̃=P^T s` 反推應為 `H̃=P^T G P`）。
   `build_channel_matrices.m` 直接採用修正後的版本，並在註解中說明推導。
4. **`get_block.m` 有兩種不同的切塊大小**：`G`/`Ȟ` 用區塊大小 `M`（區塊索引＝時間 `n`）；
   `H̃`/`H` 用區塊大小 `N`（區塊索引＝延遲 `m`）。混用會直接導致錯誤的子塊比對，
   這也是撰寫 Sim-6~9 時最容易搞混的地方。

---

## 下一步

階段二全數通過後，即可依 `chap4_list.md` 進入「階段三：off-grid（分數延遲／分數都卜勒）」
（Sim-14 ~ Sim-16），驗證分數 tap 造成的都卜勒/延遲洩漏，以及通道稀疏度 `S` 的變化。
