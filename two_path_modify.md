# Result 2a — 零假設推導與修訂後的 Test 2 規格

**對象**：`0805_W3_weeklyReport` P.38 (slide 38), Test 2 / Eq. (4.105) vs (4.118)
**日期**：2026-08-19

---

## 0. 結論先講

1. **P.38 的圖是對的，而且和你自己 slide 36 寫的 "Expected: Integer model degrades; fractional stays exact" 完全一致。**
2. **sinc 不可能是這張圖的誤差來源。** Test 2 的測試通道是 `ℓ = [0, 2, 5]`，全部是整數延遲，所以 `sinc(l − ℓᵢ) = δ[l − ℓᵢ]` 是**恆等式**，band-limited 假設在這個實驗裡沒有被使用到。我用完全不含任何內插（連 sinc 都沒出現）的 ground truth 重跑，得到和你圖上一模一樣的曲線（§1.1）。
3. **公式中真正造成 integer model 崩潰的是 Doppler 方向的 Dirichlet kernel `𝒟_N(κγ_g − k)`，它來自一個有限等比級數，是恆等式，不是近似**，跟 band-limit 無關。
4. 因此**不存在任何「不經過假設的推導」可以讓 integer 的 NMSE 小於 fractional**——這是數學上不可能，不是建模選擇問題（證明見 §3）。fractional model 的誤差恆等於 0，而 NMSE ≥ 0。
5. 但**你的直覺不是沒有道理**，只是被放在錯的指標上。integer model 真正會贏的地方是 model order / 估計變異數 / 複雜度，需要換一個實驗（§5，新的 Test 2d，有 bias–variance crossover）。
6. 順帶查出報告裡**四個真正該修的問題** F1–F4（§4），其中 F2（rect pulse 卻用 sinc kernel）和 F4（phase convention）是實質錯誤。

---

## 1. 先排除 sinc

### 1.1 零假設 ground truth 重現 P.38

我用 waveform-level 的 ground truth，**完全沒有用到內插、band-limit、sinc**（因為 ℓᵢ 是整數，`s(qTs − τᵢ) = s[q − ℓᵢ]` 精確成立）：

```
r_phys[p] = Σᵢ gᵢ · z^{κᵢ(p − ℓᵢ)} · s_tx[p − ℓᵢ]
```

參數與報告相同：`M=16, N=8, Lg=6, γ_g=1.375, g=[0.8, 0.5+0.2j, 0.3], ℓ=[0,2,5], κγ_g=[1,−2,3]+offset`。

| fractional offset | Fractional (4.105) | Integer (4.118) | 你報告的圖 |
|---|---|---|---|
| 0.000 | −296.9 dB | −296.9 dB | ≈ −289 dB |
| 0.005 | −293.8 dB | **−35.6 dB** | ≈ −35 dB |
| 0.050 | −293.5 dB | **−15.6 dB** | ≈ −15 dB |
| 0.100 | −295.6 dB | **−9.6 dB** | ≈ −9 dB |
| 0.200 | −296.1 dB | **−3.8 dB** | ≈ −4 dB |
| 0.300 | −295.4 dB | **−0.7 dB** | ≈ −1 dB |
| 0.500 | −292.5 dB | **+2.0 dB** | +3.22 dB |

完全重現。**結論：把 sinc 拿掉，圖一模一樣。** 所以 sinc 不是 Result 2a 的誤差來源。

### 1.2 sinc 在哪裡才有影響

sinc 只在 **fractional delay** 時才會作用。Result 2a 掃的是 **fractional Doppler**，兩者是正交的兩個維度。你手寫紙上 M2 的綠色結論

```
ν_{m,l}[k] = (1/N) Σᵢ gᵢ z^{κᵢ(m−l)} · sinc(l − ℓᵢ) · 𝒟_N(κᵢγ_g − k)
                                        └── delay 維 ──┘   └── Doppler 維 ──┘
```

是兩個 kernel 的乘積。Test 2 固定 `sinc(l−ℓᵢ) = δ`，只掃 `𝒟_N`。而且更關鍵的是：**Path A (ground truth)、(4.105)、(4.118) 三者用的是同一個 `sinc(l−ℓᵢ)`**，它是 common-mode 因子，換成任何其他 kernel 三條路徑會同步改變，NMSE 差值不變。

---

## 2. 從源頭的推導（每一步標明是「恆等」還是「假設」）

### S0. 連續時間物理（恆等）

```
r(t) = Σᵢ gᵢ · e^{j2πνᵢ(t − τᵢ)} · s(t − τᵢ)
s(t) = Σ_q s_tx[q] · ψ(t/Tₛ − q)          ψ = Tx pulse（任意，不假設 band-limited）
```

### S1. 接收濾波 + 取樣（恆等，給定 φ）

以 Rx filter `φ` 匹配、在 `t = q'Tₛ` 取樣：

```
r[q'] = Σᵢ gᵢ Σ_q s_tx[q] · e^{j2πνᵢ(q' − q)Tₛ} · ρᵢ(q' − q)

ρᵢ(u) ≜ ∫ ψ(x) φ*(x + u − ℓᵢ) e^{j2πνᵢ x Tₛ} dx        ← 精確的 composite kernel
```

`ρᵢ(u)` 是 Tx/Rx pulse 對的 **Doppler-tilted cross-ambiguity function**。這是「不經過任何假設」能寫到的最一般形式。注意它**同時依賴 τᵢ 和 νᵢ**。

> **A1（唯一必要的近似）**：`|νᵢ| Tₛ ≪ 1` 時 `e^{j2πνᵢxTₛ} ≈ 1`，於是
> `ρᵢ(u) → p(u − ℓᵢ) ≜ A_{ψφ}(u − ℓᵢ)`（零 Doppler 的 pulse 互相關）。
> 本系統 `ν_max Tₛ = 5.6 kHz × 4.17 µs = 0.023`，A1 造成的誤差 ~ −33 dB 以下，可接受。
> 這是**唯一**一個和 Doppler 有關的近似，而且它在 (4.105) 和 (4.118) 裡是共用的。

> **A2（pulse shape，這才是 sinc 的來源）**：
> - `ψ = φ = sinc`（理想 band-limited Tx pulse **且** 理想 band-limited Rx filter）→ `p(x) = sinc(x)`
> - `ψ = φ = rect(Tₛ)`（矩形脈波 + integrate-and-dump matched filter）→ `p(x) = tri(x) = max(0, 1−|x|)`，**只有 2 個 tap**
> - `ψ = rect, φ = δ`（矩形脈波 + 理想取樣器，無接收濾波）→ `p(x) = 1{0 ≤ x < 1}`，**1 個 tap**（等於把 ℓᵢ 直接無條件捨去）
>
> **sinc 完全是 A2 的產物。** 一般式寫成 `p(l − ℓᵢ)`，不預設任何形狀。

### S2. Delay-time channel（恆等，給定 p）

CP/ZP-OTFS 的物理取樣點 `q = m + n(M + L_g) = m + nMγ_g`：

```
h_{m,l}[n] = Σᵢ gᵢ · p(l − ℓᵢ) · z^{κᵢ(m − ℓᵢ)} · e^{j2πκᵢγ_g n/N}
                                        ↑
                    注意是 ℓᵢ（真實延遲），不是 l（tap 索引）
```

> **A3（書上 (4.6)/(4.91) 的隱藏近似）**：把相位寫成 `z^{κᵢ(m − l)}` 等於用 tap 索引 `l` 取代真實延遲 `ℓᵢ`。當 `ℓᵢ ∈ ℤ` 時 `p` 是 delta，`l ≡ ℓᵢ`，兩者相同 → **Test 2 完全不受影響**。但 fractional delay 時這會造成約 **−30 dB 的誤差地板**（實測見 §4 F4）。

### S3. DZT → Doppler 維（**純恆等式，這是重點**）

```
ν_{m,l}[k] ≜ (1/N) Σ_{p=0}^{N−1} h_{m,l}[p] e^{−j2πkp/N}
           = (1/N) Σᵢ gᵢ p(l − ℓᵢ) z^{κᵢ(m − ℓᵢ)} · Σ_{p=0}^{N−1} e^{j2π(κᵢγ_g − k)p/N}
           = (1/N) Σᵢ gᵢ p(l − ℓᵢ) z^{κᵢ(m − ℓᵢ)} · 𝒟_N(κᵢγ_g − k)
```

其中

```
𝒟_N(x) ≜ Σ_{p=0}^{N−1} e^{j2πxp/N} = (1 − e^{j2πx}) / (1 − e^{j2πx/N})
       = [sin(πx) / sin(πx/N)] · e^{jπx(N−1)/N}
```

**這是一個有限等比級數的封閉式，對任意實數 x 恆成立。** 沒有 band-limit、沒有截斷、沒有 DTFT/CTFT 的取捨、沒有任何近似。你手寫紙上紅字寫的「找 sinc 來源 / 可能有 band-limited / 用 DTFT」——這個疑慮對 **delay 維的 `sinc(l−ℓᵢ)`** 是成立的（就是 A2），但對 **Doppler 維的 `𝒟_N`** 不成立，`𝒟_N` 沒有任何可以拆掉的假設。

也因為它是恆等式，Parseval 直接成立：

```
Σ_{k=0}^{N−1} |𝒟_N(x − k)|² = N · Σ_p |e^{j2πxp/N}|² = N²    （對任意 x）
```

→ 能量守恆，這也解釋了你 Result 2b 「Total energy is conserved」的觀察。

### S4. Integer-tap model = 對上式做截斷（**近似**）

(4.118) 做的是

```
𝒟_N(κᵢγ_g − k)  ⟶  N · δ[k − round(κᵢγ_g)]
```

即：把整個 Dirichlet kernel 換成一根 delta。這不是另一個模型，**這是同一個模型刪掉 N−1 項**。

---

## 3. 為什麼 NMSE 排序在數學上不可能反轉

令 `ε = κᵢγ_g − round(κᵢγ_g)` 為分數偏移。對單一路徑：

- 訊號能量：`|𝒟_N|² 總和 = N²`
- integer model 的誤差能量：`|N − 𝒟_N(ε)|² + Σ_{k≠0} |𝒟_N(ε − k)|² = |N − 𝒟_N(ε)|² + N² − |𝒟_N(ε)|²`

化簡後得到**封閉式**：

```
NMSE_int(ε) = 2 − (2/N)·Re{𝒟_N(ε)} = 2 − 2·cos(πε(N−1)/N)·sin(πε)/(N·sin(πε/N))
```

小 ε 展開：**NMSE_int(ε) ≈ (πε)²·(1 − 1/N²)/... ≈ (πε)²**，即每 offset ×10 → NMSE +20 dB。

驗證（單路徑，N=8）：

| ε | 實測 | 封閉式 `2 − 2Re{𝒟_N}/N` |
|---|---|---|
| 0.005 | −35.60 dB | −35.69 dB |
| 0.020 | −23.56 dB | −23.65 dB |
| 0.050 | −15.62 dB | −15.71 dB |
| 0.100 | −9.65 dB | −9.74 dB |
| 0.200 | −3.86 dB | −3.95 dB |
| 0.300 | −0.73 dB | −0.80 dB |
| 0.500 | **+2.46 dB** | **+2.43 dB** |

誤差 < 0.1 dB。ε = 0.5 時 `𝒟_N(0.5)` 的實部剛好是 1（N=8），得 NMSE = 2 − 0.25 = 1.75 → **2.43 dB**，這就是你 slide 上「NMSE ≈ 2 → 3.01 dB」那段的精確版本（3.01 dB 是把 `Re{𝒟_N}` 當 0 的粗略上界，實際是 2.43 dB，多路徑時 3.22 dB）。

**不可能反轉的證明**：
S3 的 `ν_{m,l}[k]` 對任意實數 `κᵢγ_g` 恆等於 ground truth，所以 `NMSE_frac ≡ 0`。任何模型 M 的 `NMSE(M) ≥ 0`。因此
```
NMSE_int(ε) ≥ 0 = NMSE_frac,  等號僅在 ε = 0
```
沒有任何重新推導能改變這件事。要讓 integer 贏，唯一的辦法是**把 ground truth 也換成 integer-tap 通道**——那就是拿模型去比自己，不是驗證。

**「沒有 leakage 所以誤差小」的直覺錯在哪**：leakage 是**通道的物理事實**（`κᵢγ_g ∉ ℤ` 時能量真的散在所有 N 個 bin），不是模型可以選擇要不要產生的東西。integer model 不是「沒有 leakage」，是「**忽略了確實存在的 leakage**」。忽略 ≠ 消除。

---

## 4. 報告中真正需要修的四個問題

### F1 — Result 2a 的敘事錯位（文件問題）
P.33「Why sinc? …this is the delay-domain leakage」被放在解釋 Test 2 的脈絡，但 Test 2 的 `ℓ = [0,2,5]` 全是整數，sinc 在該測試中是 no-op。
**修正**：把 slide 標題改成 **"Result 2a — Doppler-only: integer-tap model collapses off-grid"**，並在圖上註明 `ℓ ∈ ℤ (sinc ≡ δ), 只掃 Doppler`。

### F2 — 系統模型宣告 rect pulse，實作卻用 sinc kernel（實質不一致）
P.32 明寫「Discrete baseband OTFS frame, **rectangular pulse shaping**」，但 `otfs_build_G.m` 用 `sinc` 內插（Eq. 4.6）。這兩者對應**不同的 A2**：

| Tx/Rx pulse | 正確 kernel `p(x)` | tap 數 | 是否 causal |
|---|---|---|---|
| ideal BL / ideal BL | `sinc(x)` | ∞（要截斷） | ✗ 非因果 |
| rect / matched rect | `tri(x)` | **2** | ✓ |
| rect / ideal sampler | `1{0≤x<1}` | **1** | ✓ |

實測（fractional delay `ℓ = [0.4, 2.4, 4.4]`, Doppler offset 0.3）：

| 真實 physics 的 pulse | 用「對應 kernel + 正確相位」的模型 | 用「sinc kernel」的模型 |
|---|---|---|
| sinc | −11.9 dB | −11.9 dB |
| **tri (rect pulse)** | **−295.5 dB（精確）** | **−6.2 dB（完全錯）** |
| **zoh (rect pulse)** | **−295.9 dB（精確）** | **−0.4 dB（完全錯）** |

→ 如果系統真的是 rect pulse shaping，現在的 `G` 是錯的，誤差高達 −6 dB。
**修正**：`otfs_build_G.m` 增加 `KERNEL` 參數（`'sinc' | 'tri' | 'zoh'`），且與 slide 上宣告的 pulse shaping 一致。

### F3 — sinc kernel 在 CP/ZP block 結構下有結構性誤差地板（新發現）
sinc 是**非因果且無限長**的，會產生 `l < 0` 的 tap（`ℓᵢ = 0.4` 時 tap 從 `l = −8` 開始）。CP-OTFS 的 block 結構無法表示負延遲：`l = −1` 在區塊邊界會取到**下一個 block 的 CP**，而不是循環回捲。結果是一個**不隨 L_g 改善**的誤差地板：

| L_g | 6 | 10 | 14 | 20 | 30 |
|---|---|---|---|---|---|
| NMSE（精確 sinc 模型 vs sinc ground truth） | −11.9 | −11.0 | −11.3 | −11.7 | −11.1 dB |

（把負 tap 用 mod-M 回捲只改善到 −18 dB，仍不精確。）

→ P.28 寫的「`L_g ≥ l_max + L_sinc`，否則會有 error floor」**不足以**解決問題；根本原因是非因果性，不是 guard 長度。
**修正**：改用 compact-support 且 causal 的 kernel（`tri`），或在 kernel 上加 bulk delay `L_sinc` 使其因果化。用 `tri` 時上表全部變成 −295 dB。

### F4 — 相位慣例 `z^{κ(m−l)}` vs `z^{κ(m−ℓᵢ)}`（實質錯誤，僅 fractional delay 時顯現）

| 真實 pulse | 正確相位 `z^{κ(m−ℓᵢ)}` | 書上相位 `z^{κ(m−l)}` |
|---|---|---|
| tri | −295.5 dB | **−31.3 dB** |
| zoh | −295.9 dB | **−29.3 dB** |

→ 這和你 README 裡已經記錄的「(4.118) phase convention 修正」是同一族問題的另一面。整數延遲時看不出來，fractional delay 時會鎖死在 −30 dB。
**修正**：`otfs_dd_closed_form.m` 的 `'fractional'` 模式，相位一律用真實 `ℓᵢ`。

---

## 5. 修訂後的 Test 2 規格

### 設計原則
原本的 Test 2 把「delay 維近似」和「Doppler 維近似」混在一個標題底下，導致敘事錯位（F1）。新規格把兩個維度**正交拆開**，並新增一個**能讓 integer model 正當獲勝**的測試（2d）。

---

#### **Test 2a-i — Doppler-only sweep**（取代現有 Result 2a）

| 項目 | 內容 |
|---|---|
| 通道 | `ℓ = [0, 2, 5]`（整數，鎖定）；`κγ_g = [1, −2, 3] + ε`，`ε ∈ [0, 0.5]` |
| Ground truth | `r[p] = Σᵢ gᵢ z^{κᵢ(p−ℓᵢ)} s_tx[p−ℓᵢ]`（**不含任何內插**） |
| 比較對象 | (4.105) Dirichlet 展開 / (4.118) integer-tap |
| 額外曲線 | **疊上封閉式 `2 − 2Re{𝒟_N(ε)}/N`（黑色虛線）** |
| 預期 | frac ≈ −290 dB 平坦；int 貼合封閉式，斜率 +20 dB/decade；ε=0.5 時 +2.4~3.2 dB |
| 通過條件 | `NMSE_frac < −250 dB` ∀ε；`|NMSE_int − closed_form| < 1 dB` |
| Slide 文案 | 標題加 "Doppler-only"；註明 `sinc ≡ δ, not exercised here` |

> 這一版把「圖為什麼長這樣」從觀察升級為**有解析預測的驗證**，比現在的版本強。

---

#### **Test 2a-ii — Delay-only sweep + kernel identifiability**（新增，處理 F2/F3）

| 項目 | 內容 |
|---|---|
| 通道 | `κγ_g = [1, −2, 3]`（整數，鎖定）；`ℓ = [0, 2, 5] + δ`，`δ ∈ [0, 0.5]` |
| Ground truth | 三種 pulse 各跑一次：`ψ ∈ {sinc, tri, zoh}` |
| 比較對象 | 三種 kernel 的模型（3×3 矩陣） |
| 預期 | **對角線（kernel 匹配）→ −290 dB；離對角線 → −6 ~ −12 dB** |
| 通過條件 | 對角線 `< −250 dB`（`tri`/`zoh`）；`sinc` 對角線因非因果性只到 ≈ −12 dB，需在報告中明列為 known limitation |
| 產出 | 3×3 NMSE 熱圖 + 「pulse shaping 決定 kernel」的一句話結論 |

---

#### **Test 2a-iii — Phase convention**（新增，處理 F4）

| 項目 | 內容 |
|---|---|
| 通道 | `ℓ = [0.4, 2.4, 4.4]`，`κγ_g = [1, −2, 3] + 0.3` |
| 比較 | `z^{κᵢ(m−ℓᵢ)}` vs `z^{κᵢ(m−l)}`，kernel 固定 `tri` |
| 預期 | 前者 −295 dB，後者 ≈ −31 dB |
| 通過條件 | 兩者差距 > 200 dB（證明 `ℓᵢ` 版本才是正確的） |

---

#### **Test 2d — 何時 integer model 才真的比較好（bias–variance）**（新增，回應你的直覺）

這是**唯一**會出現「integer 贏」的正當實驗，而且它問的是一個真正有意義的工程問題。

| 項目 | 內容 |
|---|---|
| 場景 | Channel estimation，不是 model verification |
| 通道 | `ℓ ∈ ℤ`，`κγ_g = [1,−2,3] + ε`，`ε ∈ {0.02, 0.05, 0.1, 0.3}` |
| 兩個估計器 | (a) **integer model**：每路徑估 1 個複係數 → `P` 個參數<br>(b) **fractional model**：每路徑估 `N` 個 Doppler bin → `P·N` 個參數 |
| 掃描 | pilot SNR，−10 ~ 40 dB |
| 指標 | 總 MSE（= bias² + variance），對通道係數或等效 `H_DD` |
| 預測 | bias²(int) ≈ `(πε)²`（常數，與 SNR 無關）<br>var(frac)/var(int) ≈ `N`（多 `N` 倍參數）<br>→ **crossover 在 `(πε)² ≈ (N−1)·σ²/E_p` 處** |
| 預期 | 低 SNR / 小 ε：**integer 贏**（你的直覺在此成立）<br>高 SNR / 大 ε：fractional 贏 |
| 產出 | MSE vs SNR 曲線，標出 crossover SNR |

**估計的 crossover（N=8, P=3）**：
- ε = 0.02（bias² = −24 dB）→ crossover ≈ 24 dB SNR
- ε = 0.05（bias² = −16 dB）→ crossover ≈ 16 dB SNR
- ε = 0.30（bias² = −0.8 dB）→ 幾乎全 SNR 範圍 fractional 勝

→ 這才是「integer taps 比較好」的正確陳述：**不是模型比較準，而是參數比較少、估計變異數比較小。**

---

#### **Test 2e — 複雜度 / 稀疏度對照表**（建議，一頁 slide）

| 指標 | Integer | Fractional |
|---|---|---|
| `ν_{m,l}[k]` 非零項數 | `P` | `P·N` |
| `H_DD` nnz | `O(P·MN)` | `O(P·MN·N)` |
| 每 symbol 等化複雜度 | `O(P)` | `O(P·N)` |
| 模型偏差 | `2 − 2Re{𝒟_N(ε)}/N` | `0` |

---

### 程式碼異動清單

| 檔案 | 異動 |
|---|---|
| `otfs_build_G.m` | 新增 `KERNEL ∈ {'sinc','tri','zoh'}`（預設改 `'tri'` 以符合 rect pulse 宣告）；新增 `'none'` 模式（整數延遲，跳過內插）供 Test 2a-i 用 |
| `otfs_dd_closed_form.m` | `'fractional'` 模式相位改用 `z^{κᵢ(m−ℓᵢ)}`；kernel 參數化；`𝒟_N` 抽成獨立函式並在 `x∈ℤ` 加 `mod N` 保護 |
| `otfs_verify.m` | Test 2 拆成 2a-i / 2a-ii / 2a-iii；新增 2d（需要噪聲與 pilot 設定）與 2e |
| 新增 `dirichlet_N.m` | `𝒟_N(x) = sin(πx)/sin(πx/N)·exp(jπx(N−1)/N)`，整數處回傳 `N·1{x≡0 mod N}` |
| 新增 `nmse_int_closed_form.m` | `2 − 2·real(dirichlet_N(eps,N))/N`，供 2a-i 疊圖 |

---

## 6. 一句話總結

> Result 2a 沒有問題，`𝒟_N` 是恆等式所以 fractional model 永遠精確、integer model 永遠有 `2 − 2Re{𝒟_N(ε)}/N` 的偏差，排序不可能反轉；sinc 只影響 delay 維、在本測試中是 no-op，但它確實隱藏了三個真正的問題（rect pulse 卻用 sinc kernel、sinc 非因果造成 −12 dB 地板、相位用 `l` 而非 `ℓᵢ` 造成 −30 dB 地板），這三個要修；而「integer 比較好」的直覺要放到 Test 2d 的 bias–variance 框架下才成立。