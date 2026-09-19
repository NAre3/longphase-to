# RUN-P009 獨立驗證報告

日期：2026-09-18（Asia/Taipei）
驗證者：獨立 session（未接手後續工作，僅查核先前全自動執行的結果）
驗證對象：`HANDOFF.md` 所宣稱的 RUN-P009 結果
候選 commit：`34b8d97bd48ede0bb7b2b9304e5490739df10439`

所有驗證輸出寫入 scratchpad，未覆寫任何 RUN-P0xx 證據目錄。

---

## 一、已驗證事實

### 1.1 核心輸出 174/174 位元相同 — 成立

以 `cmp` 獨立重算 29 samples × 6 outputs：

```
identical=174  different=0  missing=0  zero-byte=0
```

與 HANDOFF 宣稱一致。無空檔案、無缺檔，因此非空洞通過。

### 1.2 已產出的 14 個 checkpoint 兩側齊備 — 成立（但範圍有限）

以 14 檔名單雙向列舉 29 samples：

```
missing=0  extra=0   (29 × 14 = 406，兩側皆存在)
```

**範圍限制（重要）**：此 14 檔名單是從 `cpp_checkpoints/` 的實際內容推導出來的，
並非來自研究規格。因此它只能證明「C++ 產出的每個 checkpoint 都有 Java 對應檔」，
**不能**證明「規格要求的每個 checkpoint 都存在」。
規格要求的清單見 §3.0，該處發現一個規格要求但完全不存在的 checkpoint。

### 1.3 二進位檔可重現 — 成立

`make -C purple clean && make -C purple -j8` 後：

```
rebuild sha256 = 214d3f3210e7747cd087efbea04323136e8c4d3ae58c2f6cc0832ad001b13e2e
原有 binary    = 214d3f3210e7747cd087efbea04323136e8c4d3ae58c2f6cc0832ad001b13e2e
```

位元相同。且 `git status --untracked-files=no` 為空、HEAD = `34b8d97`，
故證據確由 `34b8d97` 的原始碼產生。

### 1.4 兩個 sample 重跑完全重現 — 成立

重跑 `HCC1937_t50_n00`（purity collapse）與 `HCC1395_NYGC_t50_n00`
（`DecimalFormat` 0.00005 邊界）：

- 新 cpp 輸出 vs 儲存的 cpp 輸出：6/6 相同
- 新 cpp 輸出 vs Java 輸出：6/6 相同
- 新 checkpoints vs 儲存的 cpp checkpoints：14/14 位元相同

### 1.4b RUN-P007 / RUN-P008 的宣稱亦經同法查核 — 成立

以同一組 `cmp` 迴圈獨立重算：

```
RUN-P007: samples=5  core identical=30/30  different=0
RUN-P008: samples=1  core identical=6/6    different=0
RUN-P006: 為 threads 1/2/4/8/16 掃描，無 per-sample 目錄（結構符合預期）
```

### 1.5 證據新鮮度 — 成立

`purple_port` mtime `21:34:23` 早於所有 cpp 輸出（`21:45`–`22:03`）。
Java `.purple.qc` mtime `21:45`–`22:03` 落在同一視窗內，
表示 RUN-P009 的 Java 是實際重跑而非沿用舊結果。

### 1.6 缺少第 30 個 sample — 成立（確認 HANDOFF #1）

`declared_samples.tsv` 為 29 列：28 validation + 1 development。
凍結輸入目錄兩側皆為 29 個 sample，且：

```
AMBER  HCC1395_HKU_t50_n00: NO
COBALT HCC1395_HKU_t50_n00: NO
```

規格要求 29 validation + 1 development，實際證據為 28 + 1。

---

## 二、實測確認的 harness 缺陷（原為「已宣稱」，本次升級為「已驗證」）

以 `H1437_t30_n20` 的複本注入故障，執行專案自身的腳本：

| 注入的故障 | evaluator 反應 | 是否被擋下 |
|---|---|---|
| 刪除一個 cpp checkpoint (`CP-P6`) | `min_rate=1.0`，僅 `n_compared` 14→13 | **否** |
| 修改非 key 數值欄 (`purity` +0.5) | `min_rate=0.0` | 是 |
| 修改 key 欄 (`start`) | `min_rate=0.99993` | 是 |

**缺檔是真實且已實測的盲點**（HANDOFF #2 屬實）。`n_compared` 雖然
會下降，但 `run_purple_validation_sample.sh` 只把
`min_checkpoint_agreement_rate` 寫進 `summary.tsv`，從不檢查 `n_compared`。

`run_purple_validation_sample.sh` 全檔唯一的 `exit` 在第 6 行（參數用法錯誤）。
比對失敗後沒有任何 exit gate，腳本以 `cat summary.tsv` 結束 → 恆為 exit 0。
**HANDOFF #3 屬實**（以讀碼確認）。

`research/scripts/compare_checkpoints.py` 第 2 行 docstring 為
`"""AMBER 移植保真度的凍結 evaluator。規則見研究規格 §5。`
—— 這是 AMBER 的 evaluator 直接指向 PURPLE 使用。**HANDOFF #5 屬實**。

---

## 三、本次新發現（HANDOFF 未記載）

### 3.0 規格要求的 CP-P10 checkpoint 完全不存在（新發現，阻擋性）

`research/specs/purple-port-fidelity-v1.md`：

- 第 171 行：`| CP-P10 | core writers | parsed rows for six gated outputs |`
- 第 216 行：`eligible_scope  CP-P1..CP-P10 + six core outputs`

實測：

```
find runs/ -name 'CP-P10*'            → 0 個檔案（所有 run、所有 sample）
grep -rn 'CP-P10' purple/*.cpp *.h    → 0 個產生點（原始碼中只找得到 CP-P9）
```

亦即規格納入判定範圍的 10 個 checkpoint 中，**CP-P10 從未被實作、從未被產出、
也從未被比對**。HANDOFF 缺陷 #5 曾提到「a parsed CP-P10 checkpoint」，
但未指出它是完全缺席而非只是比對不足。

由於 evaluator 以 `os.listdir(cpp_dir)` 列舉，且本次驗證的 manifest 亦由
C++ 輸出推導，**兩者結構上都看不到這個缺口**。

### 3.1 14 個 checkpoint 中有 1 個是恆空的，其 1.0 為除零而來

`CP-P1-cobalt-reference-pcf` 在全部 29 個 sample 中兩側皆為 0 筆資料列
（只有 header）。`compare_one` 的 `agreed / total if total else 1.0`
使得空對空得到 **1.0 而未比對任何記錄**。

因此「406 個 checkpoint 比對」實際有效比對數為 **377**，其餘 29 次為空比對。
**已查證**此為輸入所致而非 bug：COBALT 輸入目錄只含
`*.cobalt.ratio.pcf` / `*.cobalt.ratio.tsv.gz` / `*.cobalt.gc.median.tsv`，
根本不存在 reference/normal PCF 輸入檔（tumor-only）。
故兩側皆空屬預期，但彙總數字不應以 406 呈述。

各 checkpoint 位元相同比例（29 samples）：

```
CP-P1-cobalt-reference-pcf  29/29   (恆空)
CP-P1-summary               29/29
CP-P2-support-segments      29/29
其餘 11 個                   0/29   (數值格式差異，容差內一致)
```

注意：checkpoint 從未宣稱位元相同，僅核心輸出宣稱位元相同；
319 個位元相異屬預期（Java/C++ 浮點格式化不同），且全部落在
1e-9 相對容差內。

### 3.2 t00（0% 腫瘤）樣本並未觸發空 fitting 路徑

HANDOFF #4 擔心空/無效 fitting set 會走到 NaN grid。實測五個 t00 樣本：

（`purity`/`status` 依欄名自 `purple.purity.tsv` 取得；
`Method`/`QCStatus` 依欄名自 `purple.qc` 取得。）

| sample | CP-P4 | CP-P5 | CP-P6 | purity | status | qc Method | QCStatus |
|---|---|---|---|---|---|---|---|
| H1437_t00_n25 | 17288 | 15903 | 1 | 0.2200 | NORMAL | NORMAL | WARN_DELETED_GENES |
| H2009_t00_n25 | 18969 | 15903 | 1 | 0.2400 | NORMAL | NORMAL | WARN_DELETED_GENES |
| HCC1395_HKU_t00_n25 | 11742 | 15903 | 1 | 0.1800 | NORMAL | NORMAL | WARN_DELETED_GENES,WARN_LOW_PURITY |
| HCC1395_NYGC_t00_n25 | 1442 | 15903 | 1 | 0.1500 | NORMAL | NORMAL | WARN_DELETED_GENES,WARN_LOW_PURITY |
| HCC1937_t00_n25 | 18828 | 15903 | 1 | 0.2100 | NORMAL | NORMAL | WARN_DELETED_GENES |

fitting set 皆非空，Java 自身也回報 `NORMAL` 而非 `NO_TUMOR`。
**結論：缺陷 #4 在本 panel 上未被觸發，屬潛在（latent）而非已發生的錯誤。**
現有 29 sample 的通過結果不受其影響；但該路徑也因此完全未經測試。

（附帶觀察，非本次驗證目標：0% 腫瘤摻入的樣本 Java 仍給出
purity 0.15–0.24。C++ 忠實重現此行為。這是移植保真度證據，
不是純度估計正確性的證據。）

### 3.3 RUN-P009 為平行執行

29 個 sample 的 `java.runtime.tsv` 總和 2197.6 s（36.6 分），
`cpp.runtime.tsv` 總和 519.6 s（8.7 分），合計 45.3 分，
但目錄 mtime 視窗僅 21:43→22:04（約 21 分）。
**與約 2.2 倍併發執行一致（假說，未直接驗證）**；也可能有其他解釋。
已查證的部分：`grep -rl run_purple_validation_sample` 在整個 study 目錄下
只命中 HANDOFF.md 與本報告，**確實找不到任何驅動腳本**，
`events.jsonl` 亦只記錄 `RUN-P001`（P002–P005 無任何記載）。
因此無論併發與否，RUN-P009 的實際執行方式沒有留下可重現的記錄。

---

## 四、結論

1. **HANDOFF 針對 RUN-P009 所報的數字，就磁碟上的位元而言全部屬實且可獨立重現。**
   174/174 位元相同、已產出的 406 個 checkpoint 兩側齊備、binary 可位元重現、
   抽驗兩個 sample 完全重現。RUN-P007（30/30）與 RUN-P008（6/6）亦經同法查核成立。
   就本次查核涵蓋的範圍（RUN-P006/P007/P008/P009 的核心輸出與 checkpoint 檔案）
   而言，未發現虛假或遺漏的證據。此陳述不延伸到未留下檔案證據的其他自動執行過程。

2. **HANDOFF 自列的阻擋性缺陷 #1/#2/#3/#5 經實測與讀碼確認全部屬實。**
   其中 #2（缺檔盲點）與 #3（失敗仍 exit 0）已用注入故障實際示範。

3. **缺陷 #4 經查在本 panel 上未被觸發，為潛在缺陷。**

4. **新發現一項阻擋性缺口（§3.0）**：規格第 171/216 行要求的 `CP-P10`
   從未被實作、產出或比對。規格的判定範圍是 CP-P1..CP-P10，
   現有證據只涵蓋 CP-P1..CP-P9。

5. 另兩項應補入記錄：`CP-P1-cobalt-reference-pcf` 恆空、
   其 1.0 為除零而來（有效比對數 377 而非 406，§3.1）；
   RUN-P009 的執行方式未留下可重現記錄（§3.3）。

**綜合判斷：先前自動執行產出的證據本身可信且可重現，
但（a）產生證據的 harness 無法區分通過與失敗，
（b）規格範圍內的 CP-P10 完全缺席，（c）樣本數為 28+1 而非 29+1。
因此 RUN-P009 仍不足以支撐 v1 研究結論為 `SUPPORTED`。**

HANDOFF 的「Recommended continuation order」仍然成立，
但應再加入：實作並比對 CP-P10、補記 §3.1 與 §3.3。

---

## 五、決議（2026-09-18，使用者裁示）

使用者判定：**以目前的數值一致性認定移植成功**，不再追加驗證；
後續若出現問題再修。

此決議與 §四 的技術結論並不衝突，理由如下：

- §一的數值一致性證據（174/174、binary 位元重現、兩樣本完整重現）
  是本次**獨立以 `cmp` 重算**得出的，未經由有缺陷的 harness。
  harness 的缺陷影響「未來能否偵測失敗」，不影響「這批結果是否正確」。
- §二、§3.0 所列缺陷因此降級為**已知技術債**，不再是本階段的阻擋項。

### 承接的已知技術債（不再阻擋，但需保留記錄）

| 項目 | 性質 | 出處 |
|---|---|---|
| harness 缺檔盲點、失敗仍 exit 0 | 未來回歸偵測能力缺口 | §二 |
| CP-P10 從未實作／產出／比對 | 規格範圍 CP-P1..P10，實際只涵蓋 P1..P9 | §3.0 |
| 規格指定的 `compare_purple_checkpoints.py` 不存在 | 實際使用 AMBER 的通用比對器 | §3.0 |
| `checkpoint_min_logic_agreement_rate` 從未被計算 | 規格 primary endpoint | §3.0 |
| 樣本為 28 validation + 1 dev（缺 `HCC1395_HKU_t50_n00`） | 規格要求 29 + 1 | §1.6 |
| 空 fitting set 路徑未被觸發、未測試 | 潛在缺陷 | §3.2 |
| `DecimalFormat` 僅特判 `0.00005`，其他 decimal ties 未覆蓋 | 潛在缺陷 | §3.0 |

### 狀態用語界定

「移植成功」為工程判斷，依據是 29 個真實全基因組樣本上六個核心輸出
與 Java 4.4 位元相同。

規格 §5 的形式狀態 `SUPPORTED` 另有明確 gate（30/30 樣本、
全部 secondary floors、CP-P1..P10 範圍）。上述 gate 未全數滿足，
故**不應在 `result.json` 寫入 `status: SUPPORTED`**；
若要形式結案，依規格應另開 v2 slug 調整 gate，而非在 v1 原地放寬。
