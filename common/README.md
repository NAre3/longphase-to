# common/ — AMBER 與 COBALT 共用的移植核心（`namespace lp`）

這裡放的是**不知道自己在處理什麼生物量**的那一層。AMBER 餵它 BAF、COBALT 餵它 log2 ratio，
兩邊拿到的是同一份 `hmf-common` 語義。

| 檔案 | 對應的 hmftools 來源 |
|---|---|
| `Segmentation.{h,cpp}` | `common.segmentation.copynumber` 的 `Runmed`／`Gamma`／`Segmenter`（最小成本分段 DP），以及 `ChrArmLocator` 的著絲點判定 |
| `Centromeres38.inc` | `RefGenomeCoordinates.COORDS_38`，由 jar 的 `refgenome/centromeres.38.tsv` 原樣抽出 |
| `HumanChromosome.{h,cpp}` | `common.genome.chromosome.HumanChromosome` |
| `ChrBaseRegion.h` | `common.region.ChrBaseRegion` |
| `SamRecordView.{h,cpp}` | htsjdk `SAMRecord` 的座標語義 |
| `CpDump.{h,cpp}` | checkpoint dump 基礎設施（非 hmftools，用於與 Java 參考輸出逐點對照） |

**這一層不含 AMBER 或 COBALT 專屬的東西。** `segmentBafs`、`writeSegmentsFile`、
CP-A8／CP-A9 的 dump 留在 `amber/Segmentation.{h,cpp}`；COBALT 的對應物將留在 `cobalt/`。

`hmf-common/.../segmentation/` 的 19 個 Java 原始碼檔在 `amber-v4.3` 與 `cobalt-v3.0`
兩個 tag 上 git blob SHA 逐一相同，因此兩邊共用同一份 C++ 移植是有依據的，不是便宜行事。

## 這次抽出沒有改變任何行為（已驗證）

抽出前後各跑一次 AMBER，與抽出前的產物比對：13 個 checkpoint dump（CP-A1 到 CP-A9，
含各 summary 與 diagnostic）以及 3 個 stage 輸出（`amber.qc`、`amber.baf.pcf`、`amber.baf.tsv.gz`，
解壓後）全部逐位元組相同。

**驗證的是資料流的每一個切點，不只是最後的結果。** 這 13 個切點就是當初移植 AMBER 時
按資料流切出來的；只比對 stage 輸出不足以證明中間邏輯沒有被改動。

## 日後改動這一層的規矩

共用核心被兩個工具依賴，改它就可能同時動到兩邊。**每次改完都重跑一次上面的比對**
（13 個 dump 加 3 個 stage 輸出），逐位元組相同才算沒有回歸。
這是目前最便宜的護欄——AMBER 的 checkpoint 已經是現成的、涵蓋整條資料流的回歸測試。

COBALT 需要的泛化（`valuesForSegmentation` 與 `rawValues` 分離、`WindowSegments` 分支）放在
`cobalt/Segmentation.{h,cpp}`；共用核心只多了兩個輸出中間量的多載
（`gammaSegmentPenalty` 的 `GammaTrace`、`segment` 的 `pcfMeans`），原本的版本行為不變。
