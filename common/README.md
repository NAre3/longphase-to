# common/ — AMBER 與 COBALT 共用的移植核心（`namespace lp`）

這裡放的是**不知道自己在處理什麼生物量**的那一層。AMBER 餵它 BAF、COBALT 餵它 log2 ratio，
兩邊拿到的是同一份 `hmf-common` 語義。

| 檔案 | 對應的 hmftools 來源 |
|---|---|
| `Segmentation.{h,cpp}` | `common.segmentation.copynumber` 的 `Runmed`／`Gamma`／`Segmenter`（最小成本分段 DP），以及 `ChrArmLocator` 的著絲點判定 |
| `Centromeres38.inc` | `RefGenomeCoordinates.COORDS_38`，由 jar 的 `refgenome/centromeres.38.tsv` 原樣抽出 |
| `HumanChromosome.{h,cpp}` | `common.genome.chromosome.HumanChromosome` |
| `ChrBaseRegion.h` | `common.region.ChrBaseRegion` |
| `SamRecordView.{h,cpp}` | htsjdk `SAMRecord` 的座標語義（行為對照表 §4） |
| `CpDump.{h,cpp}` | checkpoint dump 基礎設施（非 hmftools，是本研究的量測接縫） |

**這一層不含 AMBER 或 COBALT 專屬的東西。** `segmentBafs`、`writeSegmentsFile`、
CP-A8／CP-A9 的 dump 留在 `amber/Segmentation.{h,cpp}`；COBALT 的對應物將留在 `cobalt/`。

`hmf-common/.../segmentation/` 的 19 個 Java 原始碼檔在 `amber-v4.3` 與 `cobalt-v3.0`
兩個 tag 上 git blob SHA 逐一相同，因此兩邊共用同一份 C++ 移植是有依據的，不是便宜行事。

## 這次抽出沒有改變任何行為（已驗證）

抽出前後以 `HCC1395_HKU_t50_n00`（dev 樣本，`-threads 16`）各跑一次，對照
`purple-port-amber-fidelity-v1/runs/RUN-011-revalidate-stdmath` 的產物：

- **13 個 checkpoint dump 逐位元組相同**：CP-A1（6,345,996 列）、CP-A2（6,342,273）、
  CP-A2-summary、CP-A2b（22,763）、CP-A3（6,342,273）、CP-A4（5,600,784）、CP-A5（724,167）、
  CP-A6、CP-A6-diagnostic、CP-A6b（701,545）、CP-A7（701,545）、CP-A8（701,590）、CP-A9（9,310）
- **3 個 stage 輸出逐位元組相同**：`amber.qc`、`amber.baf.pcf`、`amber.baf.tsv.gz`（解壓後）

**驗證的是資料流的每一個切點，不只是最後的結果。** 這 13 個切點就是當初移植 AMBER 時
按資料流切出來的；只比對 stage 輸出不足以證明中間邏輯沒有被改動。

## 日後改動這一層的規矩

共用核心被兩個工具依賴，改它就可能同時動到兩邊。**每次改完跑一次上面那組 dev 樣本比對**
（單一樣本、13 個 dump 加 3 個 stage 輸出），逐位元組相同才算沒有回歸。
這是目前最便宜的護欄——AMBER 的 checkpoint 已經是現成的、涵蓋整條資料流的回歸測試。

COBALT 需要的泛化（`valuesForSegmentation` 與 `rawValues` 分離、`WindowSegments` 分支、
gamma 參數化）尚未做，屬 `purple-port-cobalt-fidelity-v1` 的實作票。
