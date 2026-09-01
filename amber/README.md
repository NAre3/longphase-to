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
