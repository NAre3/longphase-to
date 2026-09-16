#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <sample> <run-directory>" >&2
    exit 2
fi

sample=$1
run_dir=$2
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
studies=/bip8_disk/yungjen114/longphase-to/research/studies
resources=/bip8_disk/yungjen114/run_Pacbio_purple/hmf_38/hmf_pipeline_resources.38_v2.3.0--2
reference=/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta
instrumented_jar="$repo_root/research/studies/purple-port-fidelity-v1/runs/RUN-P001/purple_v4.4.instrumented.jar"
amber="$studies/purple-port-amber-fidelity-v1/validation_reference/$sample/out"
cobalt="$studies/purple-port-cobalt-fidelity-v2/runs/RUN-C009/samples/$sample/cpp"

java_out="$run_dir/$sample/java"
java_cp="$run_dir/$sample/java_checkpoints"
cpp_out="$run_dir/$sample/cpp"
cpp_cp="$run_dir/$sample/cpp_checkpoints"
mkdir -p "$java_out" "$java_cp" "$cpp_out" "$cpp_cp"

if [[ ! -s "$java_cp/CP-P9-context-qc.tsv" || ! -s "$java_out/$sample.purple.qc" ]]; then
    /usr/bin/time -f '%e\t%M' -o "$run_dir/$sample/java.runtime.tsv" \
        java -Xmx32G -Dpurple.cpdump.dir="$java_cp" -jar "$instrumented_jar" \
        -tumor "$sample" -amber "$amber" -cobalt "$cobalt" \
        -gc_profile "$resources/dna/copy_number/GC_profile.1000bp.38.cnp" \
        -ref_genome "$reference" -ref_genome_version 38 \
        -ensembl_data_dir "$resources/common/ensembl_data" -threads 8 -no_charts \
        -output_dir "$java_out" >"$run_dir/$sample/java.stdout.log" 2>"$run_dir/$sample/java.stderr.log"
fi

/usr/bin/time -f '%e\t%M' -o "$run_dir/$sample/cpp.runtime.tsv" \
    "$repo_root/purple/purple_port" -tumor "$sample" -amber "$amber" -cobalt "$cobalt" \
    -ref_genome "$reference" -ensembl_data_dir "$resources/common/ensembl_data" \
    -threads 8 -cpdump_dir "$cpp_cp" -output_dir "$cpp_out" \
    >"$run_dir/$sample/cpp.stdout.log" 2>"$run_dir/$sample/cpp.stderr.log"

python3 "$repo_root/research/scripts/compare_checkpoints.py" \
    --ref-dir "$java_cp" --cpp-dir "$cpp_cp" \
    --out "$run_dir/$sample/checkpoint_compare.json" --examples 20 \
    >"$run_dir/$sample/checkpoint_compare.stdout.log" \
    2>"$run_dir/$sample/checkpoint_compare.stderr.log"

core_pass=true
: > "$run_dir/$sample/core_output_mismatches.txt"
for suffix in purity.tsv purity.range.tsv segment.tsv cnv.somatic.tsv chromosome_arm.tsv qc; do
    if ! cmp -s "$java_out/$sample.purple.$suffix" "$cpp_out/$sample.purple.$suffix"; then
        core_pass=false
        echo "$suffix" >> "$run_dir/$sample/core_output_mismatches.txt"
    fi
done

checkpoint_rate=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["min_checkpoint_agreement_rate"])' "$run_dir/$sample/checkpoint_compare.json")
printf 'sample\tcheckpoint_min_agreement\tcore_outputs_byte_identical\n%s\t%s\t%s\n' \
    "$sample" "$checkpoint_rate" "$core_pass" > "$run_dir/$sample/summary.tsv"
cat "$run_dir/$sample/summary.tsv"
