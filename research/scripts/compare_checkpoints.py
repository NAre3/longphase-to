#!/usr/bin/env python3
"""AMBER 移植保真度的凍結 evaluator。規則見研究規格 §5。

- 以記錄鍵做 outer join；只在一端出現的記錄一律計為不一致
- 整數與字串欄：完全相等
- 浮點欄：a == b，或 |a-b| <= max(1e-12, 1e-9 * max(|a|,|b|))；NaN 對 NaN 視為一致，
  ±Inf 需同號
- 每個 checkpoint 內 micro（記錄數加權），跨 checkpoint 取最小值

記錄鍵取該 checkpoint 存在的欄位組合：
  (chromosome, position, idx) / (arm, idx) / (field) / (chromosome, start, end)
"""
import argparse
import json
import math
import os
import sys

KEY_CANDIDATES = [
    ("chromosome", "position", "idx"),
    ("chromosome", "start", "end"),
    ("chromosome", "position"),
    ("arm", "idx"),
    ("field",),
]


def read_tsv(path):
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in handle]
    return header, rows


def pick_key(header):
    for candidate in KEY_CANDIDATES:
        if all(column in header for column in candidate):
            return candidate
    # 無可辨識的鍵時退化為整列，等同逐列嚴格比對
    return tuple(header)


def as_float(text):
    try:
        return float(text)
    except ValueError:
        return None


def values_agree(a, b):
    if a == b:
        return True
    fa, fb = as_float(a), as_float(b)
    if fa is None or fb is None:
        return False
    if math.isnan(fa) and math.isnan(fb):
        return True
    if math.isinf(fa) or math.isinf(fb):
        return fa == fb
    return abs(fa - fb) <= max(1e-12, 1e-9 * max(abs(fa), abs(fb)))


def compare_one(ref_path, cpp_path, examples):
    ref_header, ref_rows = read_tsv(ref_path)
    cpp_header, cpp_rows = read_tsv(cpp_path)

    if ref_header != cpp_header:
        return {
            "agreement_rate": 0.0,
            "reason": "header mismatch",
            "ref_header": ref_header,
            "cpp_header": cpp_header,
        }

    # Summary/context checkpoints contain exactly one record and intentionally
    # have no locus key.  Use the singleton itself as the record identity so
    # numeric fields are still compared with the frozen floating tolerance.
    key_columns = () if len(ref_rows) <= 1 and len(cpp_rows) <= 1 else pick_key(ref_header)
    key_index = [ref_header.index(column) for column in key_columns]

    def index_rows(rows):
        table = {}
        for row in rows:
            key = tuple(row[i] for i in key_index)
            table.setdefault(key, []).append(row)
        return table

    ref_table, cpp_table = index_rows(ref_rows), index_rows(cpp_rows)
    all_keys = set(ref_table) | set(cpp_table)

    agreed = 0
    mismatches = []
    for key in all_keys:
        a, b = ref_table.get(key), cpp_table.get(key)
        if a is None or b is None or len(a) != len(b):
            mismatches.append((key, "only-in-ref" if b is None else "only-in-cpp"))
            continue
        if all(all(values_agree(x, y) for x, y in zip(ra, rb)) for ra, rb in zip(a, b)):
            agreed += len(a)
        else:
            mismatches.append((key, "value"))

    total = max(len(ref_rows), len(cpp_rows), len(all_keys))
    for key, reason in mismatches[:examples]:
        print(f"  mismatch {key} ({reason})", file=sys.stderr)

    return {
        "key_columns": list(key_columns),
        "ref_rows": len(ref_rows),
        "cpp_rows": len(cpp_rows),
        "agreed": agreed,
        "mismatched_keys": len(mismatches),
        "agreement_rate": agreed / total if total else 1.0,
        "record_count_delta": abs(len(ref_rows) - len(cpp_rows)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref-dir", required=True)
    parser.add_argument("--cpp-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--examples", type=int, default=20)
    args = parser.parse_args()

    names = sorted(f for f in os.listdir(args.cpp_dir) if f.endswith(".tsv"))
    result = {}
    for name in names:
        ref_path = os.path.join(args.ref_dir, name)
        if not os.path.exists(ref_path):
            result[name] = {"agreement_rate": 0.0, "reason": "missing in reference"}
            continue
        print(f"comparing {name}", file=sys.stderr)
        result[name] = compare_one(ref_path, os.path.join(args.cpp_dir, name), args.examples)

    rates = [v["agreement_rate"] for v in result.values()]
    summary = {
        "per_checkpoint": result,
        "checkpoints_compared": len(result),
        "min_checkpoint_agreement_rate": min(rates) if rates else 0.0,
        "checkpoint_record_count_delta_max": max(
            (v.get("record_count_delta", 0) for v in result.values()), default=0),
    }
    with open(args.out, "w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps({k: v for k, v in summary.items() if k != "per_checkpoint"}, indent=2))


if __name__ == "__main__":
    main()
