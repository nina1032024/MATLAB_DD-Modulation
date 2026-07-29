# Chapter 4 OTFS 模擬 — 階段一：無通道自我一致性檢查

本資料夾實作 `chap4_list.md` 中 **階段一（Sim-1 ~ Sim-4）** 的程式碼：在加入通道之前，
先確認 modulator／demodulator 本身的實作正確、彼此等價、且沒有寫錯正規化係數。

---

## 檔案總覽

| 檔案 | 說明 |
|---|---|
| `otfs_params.m` | 產生系統參數 struct（`M, N, Δf, T, fc, F_M, F_N, tol, ...`），對應清單 B.1 |
| `unitary_dft.m` | 建立么正（unitary）DFT 矩陣 `F_K`，不依賴 Signal Processing Toolbox 的 `dftmtx` |
| `relerr.m` | 共用的相對誤差函式 `norm(a-b)/norm(b)`，`b≈0` 時退回絕對誤差 |
| `gen_dd_symbols.m` | 產生 DD 域符號矩陣 `X (M×N)` 與其向量化 `x`（QPSK/16QAM/BPSK，`E_s=1`） |
| `gen_perm_matrix.m` | 建立 row-column interleaver 置換矩陣 `P`，**用兩種獨立方法建構**以互相驗證 |
| `otfs_modulate.m` | 調變器：`X → s`，內建 **4 種等價實作**（`isfft` / `idzt` / `elem` / `vec`） |
| `otfs_demodulate.m` | 解調器：`r → Y`，內建 **4 種等價實作**（`wigner` / `dzt` / `elem` / `vec`） |
| `sim1_perfect_reconstruction.m` | **Sim-1**：無通道無雜訊時 `Y == X` |
| `sim2_modulation_paths.m` | **Sim-2**：ISFFT+Heisenberg 路徑 ≡ IDZT 路徑（含元素式、向量式） |
| `sim3_permutation_matrix.m` | **Sim-3**：`P` 的正確性（`PᵀP=I`）與 Fig. 4.10 示意圖 |
| `sim4_parseval.m` | **Sim-4**：Parseval / 各域能量守恆 |
| `run_stage1.m` | **總驅動程式**：依序執行 Sim-1~Sim-4，輸出誤差總表 |

---

## 如何執行

```matlab
cd chap4_matlab/Stage1

% 最快的檢查（M=8, N=6，驗證用尺寸，無圖）
T = run_stage1();

% 附上所有圖（spy(P)、Fig.4.10 示意圖、|X| vs |Y| 等）
T = run_stage1(otfs_params(), true);

% 換成效能用尺寸
T = run_stage1(otfs_params('M', 32, 'N', 16));
```

`run_stage1` 會回傳一張表（`table`，若無 Statistics/base 支援則退回 struct array），
可直接複製貼到進度報告：

```
                            Test                             RelError      Tol      Result
    ____________________________________________________    __________    _____    ________

    {'Sim-1 IDZT/DZT perfect reconstruction'           }    2.2061e-15    1e-12    {'PASS'}
    {'Sim-2 modulation / demodulation path equivalence'}    1.0234e-15    1e-12    {'PASS'}
    {'Sim-3 permutation matrix P'                      }    1.2399e-16    1e-12    {'PASS'}
    {'Sim-4 Parseval / energy conservation'            }    6.9485e-15    1e-12    {'PASS'}
```

四項全數通過（相對誤差 ~1e-15 量級，遠低於 `1e-12` 門檻），代表 modulator/demodulator
的實作彼此一致，可以放心進入階段二（加入通道）。

也可以單獨執行任一個 Sim（例如只想看 Sim-3 的圖）：

```matlab
p = otfs_params();
r3 = sim3_permutation_matrix(p, true);   % 第二參數 = 是否畫圖
```

---

## 各 Sim 的重點說明

### Sim-1　IDZT/DZT 完美重建
- **理論**：無通道、無雜訊時 `r = s`，理應 `Y = X`。
- **做法**：對調變器的 4 種實作 × 解調器的 4 種實作（共 16 種組合）全部跑一次
  `X → s → r=s → Y`，避免只測到「剛好互相抵消」的單一組合。
- **通過準則**：`norm(Y-X,'fro') < 1e-12`。

### Sim-2　兩條調變/解調路徑等價
- **理論**：Fig. 4.1（ISFFT + Heisenberg）≡ Fig. 4.3（IDZT）；接收端 (4.24)(4.27)(4.30) 亦同。
- **做法**：以 IDZT／DZT 為基準，比對 ISFFT 路徑、元素式迴圈、向量式（含 `P`）三種寫法，
  並額外驗證 `X_tf = F_M X F_N'`、`F_Mᵀ F_M = I` 等中間恆等式。
- **通過準則**：所有路徑對比基準的相對誤差 `< 1e-12`。

### Sim-3　置換矩陣 P 正確性
- **理論**：`(4.33)(4.34)` 的區塊建構法 vs.「先寫成 N×M 矩陣、再逐列讀出」的直覺建構法，
  兩者必須完全相同；且 `PᵀP = I_NM`（正交矩陣）。
- **做法**：`gen_perm_matrix.m` 內部用兩種獨立方法各建一次 `P`，直接相減應為 **精確 0**
  （不是誤差很小，是完全相等，因為兩者都是 0/1 矩陣）。同時驗證 `(4.32)` 與 `(4.53)` 兩個
  向量恆等式。
- **附圖**：`spy(P)`（在跑 `run_stage1(p, true)` 用的模擬尺寸下）以及 **M=3, N=4** 的
  Fig. 4.10 row-column interleaver 示意圖（左圖：DD 域寫入順序；右圖：時域讀出順序）。

### Sim-4　Parseval / 功率檢查
- **理論**：`ISFFT`、`IDZT`、`Heisenberg`（矩形脈波）、`P` 全部是么正變換，因此
  `‖x‖² = ‖s‖² = ‖x̌‖² = ‖x̃‖²`（`E_s=1` 正規化下應等於 `NM`）。
- **做法**：把同一個 frame 表示成 DD 域向量 `x`、時域向量 `s`、時頻域向量 `x̌ (4.41)`、
  延遲-時間域向量 `x̃ (4.53)`，逐一量測能量；並直接檢查 `F_M, F_N, P` 及組合算子的么正性。
- **通過準則**：所有能量之間的相對誤差 `< 1e-12`。

---

## 設計上的重點提醒（對應 `chap4_list.md` 附錄 E.3 的地雷）

1. **么正 DFT**：一律使用 `unitary_dft.m`（`F_K = exp(-j2πab/K)/√K`），不要混用未正規化的
   `dftmtx`／`fft`，否則 Sim-1 會直接 fail。
2. **兩種向量化方向不同**：`x`（`x[n+mN]=X[m,n]`，Doppler 快）與 `s`（`s[m+nM]`，delay 快）
   不是同一種 vec，`P` 存在的唯一理由就是接住這個差異。`otfs_modulate.m` / `otfs_demodulate.m`
   內的 `'vec'` method 就是照這個順序寫的。
3. **輸出一律為 NM×1 行向量**：不像書上 Appendix C code 3 Method 1 會不小心產生列向量
   （導致後面要寫 `G*s.'`），這裡所有調變器輸出都統一成行向量，避免日後串接通道矩陣時出錯。
4. **`z = e^{j2π/(NM)}`**：分母是 `NM` 不是 `N`，已在 `otfs_params.m` 中正確設定（下一階段
   通道模擬會用到）。

---

## 下一步

階段一全數通過後，即可依 `chap4_list.md` 的「階段二：on-grid 通道」繼續實作
`gen_channel_taps.m`、`gen_gs.m`、`gen_G.m`（Sim-5 ~ Sim-13），驗證時域／時頻域／
延遲-時間域／延遲-都卜勒域四個 domain 的 I/O 關係。
