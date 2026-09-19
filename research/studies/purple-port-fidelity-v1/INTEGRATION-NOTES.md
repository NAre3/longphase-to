# PURPLE 整合到 longphase-to：實作與驗證記錄

日期：2026-09-19（Asia/Taipei）
分支：`research/purple-port-fidelity-v1/run-p001`
起點：`34b8d97`（凍結候選）

使用者裁示（2026-09-18）：
- 移植依數值一致性認定成功，缺陷降為技術債（見 `VERIFICATION-2026-09-18.md` §五）
- 整合採**記憶體接手**，中間檔案預設不落地，要檔案再加 CLI 參數
- `amber_port` / `cobalt_port` / `purple_port` 保留，作為驗證行為未改變的依據
- `-ref_genome`「沒必要就不用放」
- purity 參數可設化、CLI 命名細節，暫不處理

---

## 一、提交

| commit | 內容 |
|---|---|
| `770bfe7` | 抽出 `PurplePipeline`（`loadInputs` / `runFromInputs`） |
| `474184b` | PCF 區間→位置邏輯共用（`buildPcfPositions`） |
| `bcb2a80` | 記憶體 adapter（`PurpleInputAdapter`） |
| `14f028b` | 三支 Makefile 各自 `build/` 物件目錄 |
| `6339506` | 染色體長度改由呼叫端提供 |
| `07a9719` | AMBER/COBALT postscan 回傳結果 + `writeStageOutputs` |
| `ead86a8` | 接線進 longphase-to + CLI |

---

## 二、已驗證事實

### 2.1 `purple_port` 行為全程未變

每一次改動後都重跑 `HCC1937_t50_n00` 與 `HCC1395_NYGC_t50_n00`：

```
core outputs vs Java (RUN-P009)   6/6  位元相同
checkpoints  vs 凍結 C++ 輸出      14/14 位元相同
```

四次改動（`770bfe7`、`474184b`、`6339506`、Makefile 隔離）後皆如此。

### 2.2 建置物件隔離生效

三支 target 的物件改放各自 `build/`、`build/common/`，`common/` 保持乾淨。
交錯建置實測：

```
make -C purple clean && make -C purple
  amber_port  sha256 前16碼 0a1901e7a0a229f2 -> 0a1901e7a0a229f2  未受影響
  cobalt_port sha256 前16碼 d15341a67ad3919b -> d15341a67ad3919b  未受影響
```

### 2.3 AMBER postscan 改動未改變行為

`amber_port` 重跑 `HCC1937_t50_n00`（2m13s，rc=0），對凍結 Java 參考：

```
amber.qc          位元相同
amber.baf.pcf     位元相同
amber.baf.tsv.gz  解壓後 sha256 相同（22520851e67ab3e3e4d651740fec417cd5eefbc7c8bc54261035f4b6421b3fbf）
```

`.tsv.gz` 壓縮檔差 **1 個位元組**：byte 10 = gzip 的 OS 欄位（`0xff` unknown vs `0x03` Unix），
由 zlib 寫入，與內容無關。RUN-I005 的既有比對腳本對 gz 檔同樣以 `zcat` 比對，
故此處理方式與既有作法一致。

掃描階段的中間數字亦相符：`noise floor applied: 580618 of 622586 retained`，
而凍結參考的 BAF 筆數正是 580,618。

### 2.4 COBALT postscan 改動未改變行為

`cobalt_port` 重跑 `HCC1937_t50_n00`（5m13s，rc=0），對凍結 C++ 輸出：

```
cobalt.ratio.tsv.gz    位元相同
cobalt.ratio.pcf       位元相同
cobalt.gc.median.tsv   位元相同
```

### 2.5 `-ref_genome` 可省略的前提

`readChromosomeLengths` 只讀 `.fai` 的前兩欄，從不讀序列。
BAM header `@SQ` 提供同一份資訊（`cobalt/CobaltPipeline.cpp:37-54` 的 `loadChromosomes`）。
實測 `HCC1937_t50_n00`：

```
比對 195 條 contig  長度相同 195  不同 0
```

### 2.6 round-trip 對檔案來源值恆等（必要條件）

`PurpleInputAdapter` 的四位小數 round-trip，對凍結檔案中的實際數值：

```
COBALT tumorGCRatio      3,088,257 個值   非恆等 0
COBALT tumorReadDepth    3,088,257 個值   非恆等 0
COBALT referenceGCRatio  3,088,257 個值   非恆等 0
AMBER  tumorBAF            580,618 個值   非恆等 0
AMBER  normalBAF           580,618 個值   非恆等 0
```

### 2.7 整合版可建可執行

`longphase-to` 連入全部 9 個 PURPLE 物件，CLI 與參數檢查正常：

```
--ensembl_data_dir=DIR   Ensembl data cache. enables PURPLE; requires the AMBER and COBALT options too.
--purple-output-dir=DIR  write the six PURPLE core outputs.
--purple-sample=NAME     sample id used for the PURPLE output file names.

缺 AMBER/COBALT 時：
phase: --ensembl_data_dir enables PURPLE, which needs both the AMBER (--amber-loci)
       and COBALT (--cobalt-gc-profile) options.
```

---

## 三、端到端驗證（已通過）

**「記憶體接手等同寫檔再讀回」已由實測確立。**

`scratchpad/e2e_purple.sh`，2026-09-19 22:53–23:05，wall 723s，rc=0。
整合版 `longphase-to` 一個 process 跑完 AMBER + COBALT + PURPLE，
PURPLE 走記憶體接手，未經過五個 stage 檔案。

主要判定 —— PURPLE 六個核心輸出 vs `RUN-P009` 的 Java PURPLE 輸出：

```
purity.tsv          位元相同
purity.range.tsv    位元相同
segment.tsv         位元相同
cnv.somatic.tsv     位元相同
chromosome_arm.tsv  位元相同
qc                  位元相同
=> 6/6
```

對照 —— 上游 stage 輸出 vs 凍結候選（同一次執行落地的中間檔案）：

```
=> stage 6/6
```

上游 6/6 排除了「差異來自共用掃描」的可能，因此 PURPLE 的 6/6 直接支持
adapter 的等價主張，而非被上游的一致性掩蓋。

執行過程中的中間量亦逐項相符：

| 量 | 整合版 | 凍結候選 |
|---|---|---|
| AMBER PCF 分段輸入值數 | 580,618 | 580,618 |
| COBALT 分段 arms / totalCount / gamma | 43 / 2,533,912 / 100 | 43 / 2,533,912 / 100 |
| PURPLE BAF / ratio 筆數 | 580,618 / 3,088,257 | 580,618 / 3,088,257 |
| purity/ploidy 候選數 | 15,903 | 15,903 |

註：共用掃描回報 `consumed 5689345 reads`，獨立 `amber_port` 為 `6247452`。
這是 design.md D3 已記載的語義差異（前者為通過 slicer filter 的 read 數，
後者為 per-region 造訪次數），該欄位不進任何 checkpoint 或 stage 輸出。

判讀基礎：整合研究 **RUN-I005** 已確立共用掃描的 stage 輸出與獨立 port
逐位元組相同（`stage_output_byte_identity_rate = 1.0`，60/60），
且 `HCC1937_t50_n00` 正是其 10 個樣本之一。推理鏈：

```
共用掃描 AMBER ≡ amber_port ≡ Java AMBER = RUN-P009 的 PURPLE 輸入
```

### 跨樣本重複（2026-09-20）

單樣本結論已由 10 個樣本的重複支持。`scratchpad/e2e_batch.sh`，
2026-09-19 23:25 – 2026-09-20 01:50。

| sample | secs | PURPLE | 上游 |
|---|---|---|---|
| HCC1937_t50_n00 | 723 | 6/6 | 6/6 |
| HCC1395_NYGC_t50_n00 | 1066 | 6/6 | 6/6 |
| H1437_t30_n20 | 1235 | 6/6 | 6/6 |
| HCC1395_HKU_t10_n40 | 986 | 6/6 | 6/6 |
| HCC1395_HKU_t40_n10 | 914 | 6/6 | 6/6 |
| HCC1395_NYGC_t7_n29 | 820 | 6/6 | 6/6 |
| H1437_t10_n40 | 1288 | 6/6 | 6/6 |
| H2009_t8_n30 | 809 | 6/6 | 6/6 |
| H2009_t50_n00 | 789 | 6/6 | 6/6 |
| HCC1937_t6_n26 | 790 | 6/6 | 6/6 |

**10/10 樣本 PURPLE 6/6，且 10/10 上游 6/6。** 合計 60 個核心輸出、60 個 stage 輸出，
無一例外。

涵蓋的條件：

| 維度 | 範圍 |
|---|---|
| 細胞株 | H1437、H2009、HCC1395_HKU、HCC1395_NYGC、HCC1937（全 5 條） |
| 性別 | MALE 6、FEMALE 4 |
| purity | 0.20 – 0.93 |
| ploidy | 1.46 – 4.65 |

性別兩端都涵蓋，這對 `34b8d97` 的男性 X/Y baseline 與女性 Y 排除修正是必要的；
purity 兩端（0.20 的低腫瘤比例、0.93 的 purity collapse 案例）亦都走到。

因此可陳述為：**在上述 10 個樣本涵蓋的條件下，PURPLE 的記憶體接手與寫檔再讀回
產生逐位元組相同的六個核心輸出。**

仍未涵蓋（沿用 v1 規格的既有範圍限制，非本次新增）：
- t00（0% 腫瘤摻入）樣本未納入本批次——它們不在 RUN-I005 的 10 組內
- 非 `chr` 前綴的染色體命名
- gamma=1000 的上游輸入（見 `VERIFICATION-2026-09-18.md` 與 AMBER spec A2b）

## 四、過程中的一個錯誤（已修正）

在此 worktree 建 htslib 時跳過了 `./configure`，直接 `make -C htslib lib-static`。
它「成功」產出一個壞掉的 `libhts.a`（7.55 MB，正確為 8.78 MB），連結得過但
`amber_port` 執行時 **segfault**，且崩在 BAM 掃描階段，表面上完全像是剛改的程式碼有 bug。

已照 `autoreconf -i && ./configure` 重建，並更新對應的記憶條目。
判斷方法：`ls htslib/config.mk` 不存在即表示該份 htslib 未 configure 過。

`purple_port` 不連 htslib，不受此問題影響；`amber_port` / `cobalt_port` 受影響。
