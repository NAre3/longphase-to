# AMBER tumor-only 的 C++ 移植

對照對象：hmftools tag `amber-v4.3` 的 `AmberApplication` tumor-only 路徑。

本階段（保真度驗證）**刻意不共用 LongPhase-TO 既有的 BAM 掃描，也不接進主程式的建置**，
以便單獨編、單獨跑、單獨與 Java 逐值對照。整合是第二階段的事。

```
make
./amber_port -loci <AmberGermlineSites.38.tsv.gz> \
             -tumor_only_excluded_bed <tumorOnlyExcludedSnp.38.bed> \
             -cpdump_dir <dir>
```

`tumorOnlyExcludedSnp.38.bed` 必須是**從對照用的那個 jar 內抽出**的那一份，
不能用 HMF bundle 裡的其他 bed——兩者若不同版，比對就不是在比實作。

## 目前涵蓋的 checkpoint

| Checkpoint | 內容 | 對應 Java |
| --- | --- | --- |
| CP-A1 | germline site 載入 | `AmberApplication.loadAmberSites()` → `AmberSitesFile.loadFile()` |
| CP-A2 | tumor-only blacklist 過濾 | `AmberApplication.hetLociTumorOnly()` |

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
