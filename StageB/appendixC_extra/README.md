# Appendix C 補充驗證（Stage 2 子資料夾）

本資料夾是對照書本 **Appendix C（`matlab_code.pdf`）程式碼** 與 `chap4_list.md` 附錄 E 對照表後，
額外補上的驗證項目——這些內容**不在** `chap4_list.md` 原本編號的 Sim-1~16 清單中，而是
Appendix C 程式碼裡「有程式、但先前 Stage1~3 都沒有實際跑過/驗證過」的部分。因為都是直接
沿用 Stage 2 已驗證過的 on-grid 通道機制（`gen_gs.m`、`gen_G.m`、`build_channel_matrices.m`
等），所以獨立成 `Stage2/appendixC_extra/` 子資料夾，跟 Sim-5~13 的官方項目分開放，避免混淆。

---

## 為什麼要補這些

對照 `chap4_list.md` 附錄 E.1 的「Code ↔ Sim 對照表」，可以發現有幾段 Appendix C 程式碼
**只對應到「B 參數表」或「B.2 通道參數」，沒有對應到任何一個編號 Sim**，代表先前的 Stage1~3
從未真正執行/測試過它們：

| Appendix C 程式碼 | 內容 | 先前狀態 |
|---|---|---|
| **Code 4** | 由最大車速算 `ν_max`, `k_max` | 從未實作成獨立函式 |
| **Code 5+6** | EPA/EVA/ETU 3GPP 標準通道 | 從未使用過；Stage2/3 全部只用合成 taps |
| **Code 3 Method 3** | `s=kron(Fn',Im)*P*x`（另一種運算子順序） | Sim-2 只測過 Method 1/2，沒測過這個順序 |
| **Code 12 Method 3** | `Y=(P.')*kron(Fn,Im)*r`（另一種順序） | 同上，解調端也沒測過 |
| **Code 10（4 種 r 算法）+ Code 12（3 種 Y 算法）整合成一張表** | chap4_list.md 附錄 E.4 明確要求的「誤差總表」 | Sim-5/8/9 只各自驗證單一 domain 關係，從未把 4×3 種算法整理成同一張表 |

以上就是本資料夾要補的內容。

---

## 檔案總覽

| 檔案 | 對應 Appendix C | 說明 |
|---|---|---|
| `mobility_params.m` | Code 4 | 由最大車速（km/h）算 `ν_max`（Hz）與正規化 `k_max`，並檢查 `k_max<N/2`（Sec. H.8 的走離失真門檻） |
| `gen_3gpp_channel.m` | Code 5 + 6 | 產生 EPA/EVA/ETU 標準通道：dB→線性 PDP 正規化、Rayleigh 複數增益、delay 量化到 grid（`round`）、Jakes 頻譜的都卜勒 tap |
| `sim_appendixC_channel_models.m` | Code 5+6 驗證 | 驗證 3GPP 通道模型，含 Sec. H.5/E.3.5 提到的「delay 解析度不足導致路徑塌陷」現象 |
| `sim_appendixC_full_consistency.m` | Code 3+8+9+10+12 | 附錄 E.4 要求的「r 的 4 種算法 × Y 的 3 種算法」整合誤差表，並新增 Method-3 運算子順序的交叉驗證 |

> 因為現在資料夾多了一層（`Stage2/appendixC_extra/`），這兩支 Sim 檔案開頭的 `addpath`
> 多加了一層 `..`：`addpath('../../Stage1')` 找 Stage1 函式，`addpath('..')` 找 Stage2 本身的
> `gen_channel_taps/gen_gs/gen_G/build_channel_matrices/get_block/apply_channel_conv`。

---

## 如何執行

```matlab
cd chap4_matlab/Stage2/appendixC_extra

r1 = sim_appendixC_full_consistency();     % 4×3 誤差總表 + Method-3 順序驗證
r2 = sim_appendixC_channel_models();       % 3GPP 標準通道驗證
```

（兩支都預設 `doPlot=true`，會各自跳出圖：`r1` 顯示 4×4/3×3 誤差矩陣的 heatmap，`r2` 顯示
distinct delay tap 數量 vs. `M` 的曲線圖。)

---

## 各項重點說明

### `sim_appendixC_full_consistency.m`——附錄 E.4「最小可跑驗證腳本」的誤差總表

實測結果（皆通過，最大跨方法誤差 `~1e-15`）：

```
modulation Method 3 (kron(Fn',Im)*P*x, NEW) vs Method 1 (idzt) : 1.240e-16

r (code 8-10) 兩兩相對誤差：
                M1:TDL-conv   M2:G*s        M3:via-Htilde M4:via-H
M1:TDL-conv     0             1.02e-16      1.02e-16      1.01e-15
M2:G*s          1.02e-16      0             3.28e-17      1.02e-15
M3:via-Htilde   1.02e-16      3.28e-17      0             1.02e-15
M4:via-H        1.01e-15      1.02e-15      1.02e-15      0

demod Method 3 ((P.')*kron(Fn,Im)*r, NEW) vs Method 1 (dzt) : 1.152e-16

Y (code 12) 兩兩相對誤差：
                        M1:dzt      M2:vec(P after)  M3:vec(P before)[NEW]
M1:dzt                  0           1.15e-16         1.15e-16
M2:vec(P after)         1.15e-16    0                0
M3:vec(P before)[NEW]   1.15e-16    0                0

-> 全部跨方法最大誤差 = 1.017e-15  -> PASS
```

**這張表直接對應 `chap4_list.md` 附錄 E.4 的要求**：「你手上就有 `r` 的 4 種算法 × `Y` 的 3
種算法 + 1 個理論值（Sim-10 已驗證，這裡不重複），任兩兩相減都應該 `<1e-10`」。四種 `r` 算法
與三種 `Y` 算法皆使用本專案已修正過的慣例（**不是**書本原始程式碼）：

- **`r` Method 1**（TDL 摺積）已修正 Appendix C 的 off-by-one（書上 `for ell=0:(delay_spread-1)`
  少算最後一個 tap；本專案 `apply_channel_conv.m` 從一開始就是 `ell=0:delay_spread`）。
- **`r` Method 3/4**（經 `H̃`/`H`）使用修正後的 `H̃=P^T·G·P`（書上 Code 9 第一行誤植為
  `P*G*P.'`；本專案 `build_channel_matrices.m` 已修正，見 Stage2 README 的除錯發現）。
- 所有向量統一為 **NM×1 直行向量**（書上 Code 3 Method 1 產生 1×NM 列向量，導致後面要寫
  `G*s.'`；本專案從 Stage 1 起就統一成行向量，見 Stage1 README）。

**新加入的驗證**：Code 3 / Code 12 的「Method 3」使用了跟 Method 2 不同的運算子順序
（`kron(Fn',Im)*P*x` vs. `P*kron(Im,Fn')*x`；解調端同理）。這兩種順序在數學上會透過
「置換矩陣與 Kronecker 積可交換」的恆等式（跟混合基 FFT 分解用的 stride permutation 是同一
類恆等式）保證相等，但先前 Sim-2/Sim-7 都沒有真的測試過這個特定排列方式——現在補上，
確認 `P` 矩陣在這個不同的運算順序下依然正確。

**誠實的獨立性說明**：4 種 `r` 算法中，Method 1（TDL 摺積）與 Method 2（`G*s`）是用**互相獨立**
的程式碼寫出來的（這正是 Sim-5 本來就驗證的事）；但 Method 3、4（經 `H̃`、`H`）在 `H̃`、`H`
的建構方式底下，數學上是被**迫**要跟 Method 2 一致的（`build_channel_matrices.m` 的推導可以
直接證明），所以它們是「建構正確性的迴歸測試」，不是獨立的物理驗證。不過把全部 7 種算法
整理進同一張表，正是老師（`chap4_list.md` 作者）要的整合誤差表格式，這件事本身就是先前
沒做過的。

### `sim_appendixC_channel_models.m`——3GPP 標準通道模型驗證

Appendix C 的 Code 5/6/7 提供三種通道 taps 來源：EPA/EVA/ETU（真實 3GPP 標準）或合成通道
（`chap4_list.md` 建議的驗證階段做法）。**Stage 2/3 從頭到尾都只用合成通道**（`gen_channel_taps.m`
的 `l_i=[0,1,2], k_i=[0,1,-2]`），從未真正跑過 EPA/EVA/ETU——這是本檔案要補的。

四項檢查：

1. **PDP 正規化**（gated，精確為 0）：`sum(pdp_linear)==1`，這是套用隨機衰落前、確定性的
   功率分配，理論上就該精確等於 1。
2. **Rayleigh 增益的蒙地卡羅驗證**（觀察性）：`g_i = sqrt(pdp_i)*CN(0,1)`，用 5000 次實現估計
   `E|g_i|^2` 應收斂到 `pdp_linear(i)`（統計量，隨樣本數以 `~1/√5000≈1.4%` 的速度收斂，
   不是機器精度等級，因此不設嚴格門檻）。
3. **delay 解析度不足導致路徑塌陷**（觀察性，重現 `chap4_list.md` Sec. H.5/E.3.5 的具體數字）：
   掃描 `M=[8,32,64,128,256,512]`，計算量化後**相異** `l_i` 的個數。**實測 `EVA @ M=64` 得到
   `l_i=[0 0 0 0 0 1 1 2 2]`，與 `chap4_list.md` 文中所寫的數字完全吻合**——這是先前從未
   真正驗證過、只是照抄書上文字的一個結論，現在用程式親自算出來確認無誤。`M` 增加到
   `512`（書本 Code 1 的預設值）後，9 條路徑的相異 `l_i` 數量會顯著回升，印證「要看到完整
   延遲解析度就要用大 `M`」的建議。
4. **整合測試**（gated）：把 `EVA @ M=64` 這個真實通道，餵進 Stage 2 已驗證過的完整管線
   （`gen_gs`→`gen_G`→`build_channel_matrices`），重跑 Sim-5 風格的自我一致性檢查
   （`G*s` vs. 直接摺積）與 Sim-9 風格的 circulant 檢查（僅 `m>=l` 的「安全」子塊）。
   通過的話，代表 Stage 2 的機制不只對合成 taps 有效，對真實的標準化通道模型一樣成立。

---

## 設計上的重點提醒

1. **`k_max` 的預設值刻意選 500 km/h**（`chap4_list.md` Sec. B.2 建議值），不是 Appendix C
   Code 4 範例的 100 km/h——因為 Stage 2/3 的 `otfs_params.m`/驗證慣例都是照著 `chap4_list.md`
   走，這裡保持一致。若要重現書本 Code 4 的確切數字，呼叫時傳入
   `mobility_params(p, 100)` 即可。
2. **`gen_3gpp_channel.m` 產生的 `k_i` 預設是分數（Jakes 頻譜）**，這正是 `chap4_list.md`
   Sec. E.3.4 提醒的「Code 6/7 預設是 off-grid」。這裡刻意保留分數 `k_i`（不像 Stage 2 的
   `gen_channel_taps.m` 那樣可選擇 `grid='on'` 取整數），因為 `gen_gs.m` 本來就不要求 `k_i`
   是整數（只有 `l_i` 需要是整數用來做列索引），所以可以直接驗證「真實、分數都卜勒的標準
   通道」通過 Stage 2 管線，不需要額外改寫任何函式。
3. **`sum(pdp_linear)` 是精確等於 1 的確定性檢查，`E|g_i|^2` 才是統計量**：兩者的誤差門檻
   應該分開設計（前者可以要求機器精度，後者只能要求統計收斂），程式裡明確把兩者分開報告，
   避免用同一個門檻誤導判讀。

---

## 下一步

若還要更完整對照 Appendix C，剩下的 Code 13/14（LMMSE 偵測）與 Code 17/18（CP-OTFS/ZP-OTFS
的 `G` 矩陣）分別對應 `chap4_list.md` 的 **Sim-17、Sim-18**——這兩項屬於「階段四（選做）」，
不在本次「補上 Stage1~3 沒做過的模擬」範圍內，之後若要做可以再另外建立 Stage4 資料夾。
