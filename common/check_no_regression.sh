#!/bin/bash
# 共用核心（common/）的回歸護欄。
#
# 改動 common/ 之後跑這支：以 AMBER 的 dev 樣本重跑一次，對照既有的參考產物，
# 比對資料流上的每一個 checkpoint 與三個 stage 輸出。逐位元組相同才算沒有回歸。
#
# 為什麼比 stage 輸出更嚴：AMBER 的 13 個 dump 是按資料流切出來的，
# 只比最後的結果無法排除「湊出相同結果但中間邏輯已改變」。
#
#   用法：common/check_no_regression.sh [參考目錄]
#   參考目錄預設為 RUN-011-revalidate-stdmath 的 dev 樣本產物。
set -uo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
STUDY=/bip8_disk/yungjen114/longphase-to/research/studies/purple-port-amber-fidelity-v1
SAMPLE=HCC1395_HKU_t50_n00
REF=${1:-$STUDY/runs/RUN-011-revalidate-stdmath/cpp/$SAMPLE}

BAM=/big8_disk/data/HCC1395/ONT/subsample/t50_n00/HCC1395_t50_n00.bam
LOCI=/bip8_disk/yungjen114/run_Pacbio_purple/hmf_38/hmf_pipeline_resources.38_v2.3.0--2/dna/copy_number/AmberGermlineSites.38.tsv.gz
BED=$STUDY/runs/RUN-001/jar/tumorOnlyExcludedSnp.38.bed
THREADS=${THREADS:-16}

for f in "$BAM" "$LOCI" "$BED"; do
  [ -r "$f" ] || { echo "缺少輸入：$f"; exit 2; }
done
[ -d "$REF" ] || { echo "缺少參考目錄：$REF"; exit 2; }
[ -x "$ROOT/amber/amber_port" ] || { echo "請先 build：cd $ROOT/amber && make"; exit 2; }

OUT=$(mktemp -d -t common-regression-XXXXXX)
trap 'echo "產物保留於 $OUT"' EXIT
mkdir -p "$OUT/cp" "$OUT/o"

echo "[1/2] 以 $SAMPLE 重跑（threads=$THREADS）…"
"$ROOT/amber/amber_port" -loci "$LOCI" -tumor_only_excluded_bed "$BED" -tumor_bam "$BAM" \
  -min_base_quality 13 -min_map_quality 50 -tumor "$SAMPLE" \
  -output_dir "$OUT/o" -cpdump_dir "$OUT/cp" -threads "$THREADS" > "$OUT/run.log" 2>&1
rc=$?
[ $rc -eq 0 ] || { echo "amber_port 失敗 rc=$rc，見 $OUT/run.log"; exit 1; }

echo "[2/2] 逐 checkpoint 比對…"
fail=0
DUMPS="CP-A1 CP-A2 CP-A2-summary CP-A2b CP-A3 CP-A4 CP-A5 CP-A6 CP-A6-diagnostic CP-A6b CP-A7 CP-A8 CP-A9"
for f in $DUMPS; do
  a=$OUT/cp/$f.tsv; b=$REF/cp/$f.tsv
  if [ ! -f "$a" ] || [ ! -f "$b" ]; then echo "  缺檔        $f"; fail=1; continue; fi
  if cmp -s "$a" "$b"; then printf "  相同        %-18s %s 列\n" "$f" "$(wc -l < "$a")"
  else echo "  ** 不同 **  $f"; fail=1; fi
done
for f in amber.qc amber.baf.pcf; do
  if cmp -s "$OUT/o/$SAMPLE.$f" "$REF/o/$SAMPLE.$f"; then echo "  相同        $f"
  else echo "  ** 不同 **  $f"; fail=1; fi
done
if cmp -s <(zcat "$OUT/o/$SAMPLE.amber.baf.tsv.gz") <(zcat "$REF/o/$SAMPLE.amber.baf.tsv.gz"); then
  echo "  相同        amber.baf.tsv.gz（解壓後）"
else echo "  ** 不同 **  amber.baf.tsv.gz"; fail=1; fi

echo
if [ $fail -eq 0 ]; then echo "PASS：16 項產物逐位元組相同，共用核心無回歸"; else
  echo "FAIL：有產物不同。差異的 checkpoint 指出邏輯改變的位置，先看最早出現不同的那一個。"; fi
exit $fail
