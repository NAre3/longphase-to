# AMBER tumor-only 的 C++ 移植

對照對象：hmftools tag `amber-v4.3` 的 `AmberApplication` tumor-only 路徑。

> **`research/` 與 `tools/` 底下的路徑是工作站本機的,不隨本 repo 發佈。**
> 保真度研究的行為契約、math provenance、逐 run 驗證證據放在
> `research/studies/purple-port-amber-fidelity-v1/`,診斷用的 CDF 一致性工具放在
> `amber/tools/`——兩者都已刻意排除在版控之外:前者含數十 GB 的 checkpoint dump,
> 後者是驗證器具而非產品程式碼,皆屬過程材料。本檔與原始碼註解中對這些路徑的引用
> 是指向本機記錄,clone 下來不會有這些檔案。

本階段（保真度驗證）**刻意不共用 LongPhase-TO 既有的 BAM 掃描，也不接進主程式的建置**，
以便單獨編、單獨跑、單獨與 Java 逐值對照。整合是第二階段的事。

```
make                       # 需先在 ../htslib 建好 libhts.a（autoreconf -i && ./configure && make lib-static）
./amber_port -loci <AmberGermlineSites.38.tsv.gz> \
             -tumor_only_excluded_bed <tumorOnlyExcludedSnp.38.bed> \
             -tumor_bam <tumor.bam> \
             -cpdump_dir <dir>
```

`-tumor_bam` 未給定時只跑到 CP-A2。`-debug_only_chr` 只是迭代時縮短週期用的 harness 旗標，
**不是移植的行為**（AMBER 的 `-specific_chr` 並不限制 loci），正式記錄的執行不得帶它。

BAM 讀取用的是 LongPhase-TO in-tree 的 htslib 1.16，與行為對照表引用的版本一致。

`tumorOnlyExcludedSnp.38.bed` 必須是**從對照用的那個 jar 內抽出**的那一份，
不能用 HMF bundle 裡的其他 bed——兩者若不同版，比對就不是在比實作。

## 目前涵蓋的 checkpoint

| Checkpoint | 內容 | 對應 Java |
| --- | --- | --- |
| CP-A1 | germline site 載入 | `AmberApplication.loadAmberSites()` → `AmberSitesFile.loadFile()` |
| CP-A2 | tumor-only blacklist 過濾 | `AmberApplication.hetLociTumorOnly()` |
| CP-A2b | region 切分（minGap 4000） | `BamEvidenceReader.populateTaskQueue()` |
| CP-A3 | 七個 per-locus 計數器 | `BamEvidenceReader.processBam()` → `RegionTask` → `PositionEvidenceChecker` |
| CP-A4 | IndelCount == 0 的保留集合 | `TumorAnalysis.tumorBAFAndContamination()` |
| CP-A5 | 四道 filter + 排序（含 idx） | `AmberApplication.runTumorOnly()` |
| CP-A6 | noise floor / contamination | `TumorOnlyPurityAnalysis`（amber.purity 套件） |
| CP-A6b | noise floor 的套用 | `AmberApplication.runTumorOnly()` |
| CP-A7 | AmberBAF 轉換 | `AmberUtils.fromTumorBaf()` |
| CP-A8 | 分段輸入（per-arm 陣列）與 penalty | `PerArmSegmenter` 建構 |
| CP-A9 | 分段結果 | `BAFSegmenter.writeSegments` |
| （stage） | `amber.baf.tsv.gz`、`amber.qc`、`amber.baf.pcf` | `ResultsWriter.persistBAF/persistQC` |

## 移植時逐條對齊的行為

- 非 human chromosome（含 MT 與 alt contig）在載入時即跳過；判定為去 `chr` 前綴後
  1-22 或 X 或 Y，其餘皆否
- `SnpCheck` 欄以 Java `Boolean.parseBoolean` 語義解析：僅不分大小寫的 `true` 為真
- `Frequency` 為選用欄位；bundle 版本沒有這一欄，此時頻率以 0 代入
- BED 為 0-based 半開區間，載入時 `start+1` 轉 1-based 含端點
- blacklist 為逐區間線性掃描並在命中時 break（區間僅 32 個），與 Java 相同；
  不改成區間樹——保真度階段以行為一致為先
- dump 的排序必須是 **stable**：Java 用 `List.sort`（TimSort，穩定），
  同一 `(chromosome, position)` 的多筆記錄要保留讀入順序

## CP-A2b / CP-A3 移植時逐條對齊的行為

完整依據見 `research/studies/purple-port-amber-fidelity-v1/runs/RUN-004/behaviour-contract.md`，
以下只列最容易寫錯的幾項：

- read position 一律採 htsjdk 約定：1-based、**含 soft clip、不含 hard clip**，
  等於 `bam_get_seq()` 索引加一。alignment block 每條 read 只建一次（對應 htsjdk 快取
  `mAlignmentBlocks`），因為 `getBaseQuality` 對每個 (read, 位點) 都會被呼叫
- `alignmentEnd` **不使用 `bam_endpos()`**：htslib 在 reference length 為 0 時強制 rlen=1，
  htsjdk 則得到 `alignmentStart - 1`，兩者差 1
- `readDepth` 在 filter 判定之後、return 之前無條件遞增 → 被品質擋掉的 read 仍計入
- MAPQ 與 BASEQ 兩個 filter 各自獨立判定，同一條 read 可同時計入兩個計數器
- 位點落在 deletion 上（`readIndex < 0`）時 Ref/Alt 與 `indelCount` **三者都不動**
- `getBaseQuality` 在 deletion 上往右掃到第一個有 read position 的參考位置，
  掃到 `alignmentEnd` 仍無則回傳 0
- slicer 只以四個 flag 過濾（unmapped/secondary/supplementary/duplicate），
  **不做 MAPQ 過濾**——`BamSlicer(0, ...)` 的最低品質 0 等於不過濾
- region 切分條件為嚴格小於：`end + 4000 == pos` 併入同一個 region
- `RegionTask` 的提前中止（`haltProcessing`）在本設定下不可能觸發，故未實作，理由見對照表 §2.4
- Java 端遇到非 ACGTN 的 IUPAC 碼會丟例外，C++ 改為計數並輸出，讓這個已知的不對稱點可觀察

## CP-A4 / CP-A5 移植時逐條對齊的行為

- **排序用的是染色體 rank 的數值序，不是字串序**（`ContigComparator` → `HumanChromosome.chromosomeRank`：
  1-22 取數值、X=23、Y=24、MT/M=25、其餘 26）。照字串排會讓 chr10 跑到 chr2 前面，
  72 萬筆的 idx 全錯。CP-A5 的驗收規則明文要求 idx 逐筆相同
- CP-A5 用 `CpDump::write` 而非 `writeSorted`——**順序本身就是被比對的對象**
- tumor-only 的最低深度是 `DEFAULT_TUMOR_ONLY_MIN_DEPTH` **25**，不是 tumor/normal 的 8
  （`AmberConfig.java:135-146` 依 ReferenceIds 是否為空分支）
- `aboveQualFilter` 的分母是 ReadDepth，三個 filtered 計數器相加後**嚴格小於** 0.15 才通過
- CP-A3 是全部位點，CP-A4 才是第一次縮減（634 萬 → 560 萬），CP-A5 再縮到 72 萬

## CP-A6 / CP-A6b：數值一致性是這一段的核心

AMBER 的 peak 捕捉判定是拿 binomial CDF 去比 0.16 / 0.84 兩個門檻（`CandidatePeak.java:88-93`）。
門檻判定是布林值，CDF 的偏差有機會翻轉某個點的歸屬，因此這一段的數值行為要單獨驗證。

- `CommonsMath.{h,cpp}`：commons-math3 3.6.1 的 `BinomialDistribution.cumulativeProbability`
  → `Beta.regularizedBeta`（連分數，epsilon 1E-14）。連分數展開保留照抄，因為
  **正規化不完全 beta 函數 C++ 標準庫沒有**（`std::beta` 是完全 beta，且此處
  a 小、b 可達 1000，會下溢為 0，不能替代）。其餘一律走標準庫：
  `log` / `log1p` / `exp` / `lgamma`。
- **編譯帶 `-ffp-contract=off`**：`a*b+c` 若被融合成 FMA 會少一次捨入。
  這不再是與 Java 對齊的契約，而是讓**本程式自己的輸出**不隨編譯器與 `-march` 飄移
  ——PURPLE 直接消費這些檔案。實測：baseline `x86-64` 下此旗標無作用（ISA 無 FMA 指令），
  但加上 `-march=native` 且允許收縮時，61366 組 CDF 有 5988 組改變（2026-09-06）。

> **【2026-09-06 就地修正】** 本節原本寫「必須連數值實作一起複製」「不能用 `std::log`
> 等系統版本」，並以 `FastMath.{h,cpp}` + `FastMathTables.inc`（自 jar bytecode 反射倒出的
> 查表）與 commons-math 的 `Gamma.logGamma1p` 高精度 `logBeta` 機制達成逐位元相同。
> **該敘述已不成立，相關檔案已刪除。** 逐位元相同是當時達成的結果，不是規格要求的門檻
> ——研究規格 §5 對 double 欄位宣告的容差本來就是相對 1e-9。
>
> 改用標準庫後的實測（依據見 `research/studies/purple-port-amber-fidelity-v1/math_provenance/`）：
>
> | 量 | 值 |
> |---|---|
> | CDF 與 Java 逐位元不同 | 15021 / 61366（24.5%） |
> | 最大相對誤差 | 4.11e-11（n=1000, k=9），在規格容差 1e-9 之內 |
> | 跨越 0.16 / 0.84 門檻的筆數 | 0 |
> | 全體 CDF 值距門檻最近的距離 | 6.49e-05（比最大偏差大六個數量級） |
> | 30 組樣本的 noise floor 輸出 | 逐位元組相同（含每個 grid level 的 `%.17g` score） |
>
> **未涵蓋**：以上為 30 組 ONT 樣本的實測，不是對所有輸入的證明。
> 若日後樣本的 CDF 值落到距門檻 1e-10 以內，判定仍可能翻轉。

驗證工具：`tools/CdfConformance.java`（以 jar 內的 commons-math 產生基準）+
`tools/cdf_conformance.cpp`（逐位元比對）。基準檔 `java_cdf.tsv` 已存於
`research/studies/purple-port-amber-fidelity-v1/math_provenance/`。
**其角色已從「必須為零」的驗收門檻改為偏差的量測儀器。**

其餘逐條對齊的行為：

- 著絲點座標與免疫排除區間皆取自 jar 內的資源檔／原始碼，非人工輸入
- `RegionsFilter` 的狀態式掃描照抄（只檢查第一個 end >= position 的區間），
  不改成「檢查所有區間」——後者雖在不重疊區間下等價，但那是額外假設
- `LocalMaximaFinder` 末端的 `maxima.remove(firstNonZero)` 是「移除第一個相等的元素」，
  不是「移除索引 0」
- gnomad 頻率：loci 檔無 `Frequency` 欄，全部為 0 → 檢查退化為「兩個 band 皆非空」
- `Doubles.greaterOrEqual` 是 `value - reference > -1e-10`，不是 `>=`

### 開發用 harness（非移植行為）

`tools/noisefloor_from_cp_a5` 直接載入 CP-A5.tsv 只跑 noise floor 這一段（22 秒），
避免每次迭代都重掃 BAM（21 分鐘）。與 `-debug_only_chr` 同性質，正式記錄的執行仍為完整流程。

## CP-A7 與 stage 輸出

```
./amber_port ... -tumor <sampleId> -output_dir <dir>
```
兩者同時給定時才寫出 `amber.baf.tsv.gz` 與 `amber.qc`。

逐條對齊的行為：

- **`tumorBAF` 的分母是 `Alt+Ref`，不是 `ReadDepth`**（`AmberUtils.java:60-61`）；
  `tumorDepth` 才是 `ReadDepth`
- tumor-only 沒有 normal 樣本，normal 的三個計數皆為 0
  → `normalBAF = 0/0 = NaN`，`normalDepth = 0`。CP-A7 直接寫 NaN；
  寫檔時 `AmberBAFFile.toString` 以 `Doubles.isFinite` 判斷，非有限值輸出字面 **"0"**
  （`AmberBAFFile.java:84-86`）
- `amber.baf.tsv.gz` 為 4 位小數（`DecimalFormat("0.0000")`）。
  Java 用 HALF_EVEN、C++ 用 `printf %.4f`——**兩者不保證一致**，
  故先以參考端既有輸出實測：70 萬列零差異（見 RUN-007 的 IMPLEMENTATION checkpoint）
- `amber.qc` 的 `QCStatus` 由 `Doubles.greaterThan`（epsilon 1e-10）判定，不是 `>`；
  contamination 為 0 → PASS。tumor-only 下 `ConsanguinityProportion` 為 0、
  `UniparentalDisomy` 為 `NONE`
- `persistBAF` 內另含 PCF 分段（`ResultsWriter.java:43-51`），屬 EXP-010，此處不實作

## CP-A8 / CP-A9：PCF 分段

**tumor-only 也會做分段**——`runTumorOnly()` 本身沒呼叫，但它呼叫的 `persistBAF()` 內有
（`ResultsWriter.java:43-51`）。gamma 硬編碼 100.0，AMBER 4.3 無 CLI 可調。

演算法：每個染色體臂各跑一次最小成本分段的動態規劃（O(n²)，最大的臂 5.3 萬點，
全基因體約 90 億次內層迴圈，C++ 單執行緒約 19 秒）。penalty 逐臂計算：
寬度 51 的移動中位數 → 殘差的 MAD → 平方乘 gamma。

逐條對齊的行為：

- **段的 `MeanRatio` 取的是 rawValues 該段的平均**（`Doubles.mean`），
  **不是** `PiecewiseConstantFit.means` 裡那個已四捨五入到三位小數的值——
  兩者在程式裡是不同的東西（`ChromosomeArmSegments.java:27`）
- 動態規劃在成本相同時選哪個切點，由 `cost < minCost`（嚴格小於）與
  `Double.MAX_VALUE` 的初值決定。改成 `<=` 會選到不同切點
- cumulative sum 必須照原順序累加；改用其他求和方式會動到最後一位
- `Doubles.mean` 也是依序累加後再除，不可改寫
- `WindowedMedian` 的 `getMedian()` 回傳 maxHeap 頂端，即視窗內第 ceil(w/2) 小的值。
  視窗為奇數且不超過資料長度時等於真中位數；此處以「第 k 小」實作，
  與 heap 機制無關但輸出相同
- `.pcf` 的 `MeanRatio` 用 `DecimalFormat("#.####")`：最多四位小數、**去掉尾端的 0**。
  本資料 9309 段中有 981 段不足四位，該分支確實被觸發
- arm 的排序為 `ChrArm.compareTo`：先比染色體的 enum 序（1..22, X, Y），再比 arm（P < Q）

### uniform penalty 分支：刻意不實作

`totalCount < 100000` 時 `PerArmSegmenter` 走 uniform 分支，其 `allRatios` 的組裝順序
取自 `HashMap.keySet()`，**跨語言不保證重現**（spec U5）。本研究的 dev 樣本
`totalCount = 701544`，走 per-arm-gamma，不會進該分支。

C++ 端在偵測到該分支時**明確拋錯**，而非用某個自訂順序默默算出結果——
若日後樣本的 BAF site 數低於 10 萬，這個問題必須重新處理，屆時應該要看到失敗而不是看到數字。

## 選用輸出（CLI 控制，預設關閉）

```
-write_tumor_data   寫出 <sample>.amber.tumor.raw.tsv.gz
-write_version      寫出 amber.version
```

**PURPLE 只讀三個檔**：`amber.qc`、`amber.baf.tsv.gz`、`amber.baf.pcf`
（`purple/src/main/java/com/hartwig/hmftools/purple/AmberData.java`，缺一即拋 `ParseException`）。
上面兩個選用輸出 PURPLE 不讀取，純屬稽核／除錯用途，故預設關閉。

**已實證不影響計算**，兩條獨立證據：

1. 產生 Java 參考端時本就帶著 `-write_tumor_data`，而 C++ 端在未實作此輸出的情況下，
   30 組樣本的十一個 checkpoint 與三個 stage 輸出仍逐位元組相同。
2. 同一支 C++ 在開與不開旗標下執行，十三個 dump 完全相同。

### 兩者的內容

- `amber.tumor.raw.tsv.gz` **就是 CP-A5 的 rawData**，只是換上 AMBER 的欄名
  （`RefCount`／`AltCount` 而非 `refSupport`／`altSupport`，且無 `idx` 欄）。
  已對 Java 的輸出逐位元組驗證（724,166 列）。
- `amber.version` **刻意不照抄 Java 的 `version=4.3`**：那會讓這個檔看起來像 hmftools 的產出。
  此處記錄的是「本移植重現的是哪個 AMBER 版本」。Java 的 `build.date` 是 jar 的建置時間，
  本就無法重現。
