# AMBER tumor-only 的 C++ 移植

對照對象：hmftools tag `amber-v4.3` 的 `AmberApplication` tumor-only 路徑。

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
門檻判定是布林值，**CDF 的最後一個位元不同就可能翻轉某個點的歸屬**，因此這一段不能只
「算出同一個數學函數」，必須連數值實作一起複製：

- `CommonsMath.{h,cpp}`：commons-math3 3.6.1 的 `BinomialDistribution.cumulativeProbability`
  → `Beta.regularizedBeta`（連分數，epsilon 1E-14）→ `logBeta` → `Gamma.logGamma1p`。
  Gamma 的 36 個常數由原始碼機械抽出，未經人工轉錄。
- `FastMath.{h,cpp}` + `FastMathTables.inc`：commons-math 的 `FastMath.log/log1p/exp`。
  **不能用 `std::log` 等系統版本**——實測 FastMath 與 `java.lang.Math` 本身就不逐位元相同
  （20 萬組取樣中 log1p 差 13142 組、exp 差 457 組、log 差 21 組），
  改用系統版本會讓 CDF 出現 1~32 ulp 偏差（61366 組中 215 組）。
  查表以反射自 `amber_v4.3.jar` 內的 bytecode 倒出，不是抄上游原始碼。
- **編譯必須帶 `-ffp-contract=off`**：`a*b+c` 若被融合成 FMA，中間結果少一次捨入，
  結果就與 Java 不同。此旗標已寫進 Makefile。

驗證：`tools/CdfConformance.java`（以 jar 內的 commons-math 產生基準）+
`tools/cdf_conformance.cpp`（逐位元比對）→ **61366 組零不一致**。

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
