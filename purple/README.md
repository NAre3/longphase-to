# PURPLE v4.4 tumor-only C++ port

This directory is a checkpoint-driven port of the tumor-only WGS core used by
PURPLE v4.4.

> **Paths under `research/` are workstation-local and are not distributed with
> this repository.** The frozen behaviour contract, the Java reference
> artifacts and the per-run validation evidence live in
> `research/studies/purple-port-fidelity-v1/`, which is deliberately
> gitignored: it holds tens of GB of checkpoint dumps, and the study records
> are process material rather than product code. References to those paths in
> this file and in source comments are pointers to the local record, not to
> files you will find in a clone.

The current executable implements CP-P1 through CP-P10: input parsing and
tumor-only adjustment, support segmentation, observed-region aggregation,
fitting-region selection, the full purity/ploidy score grid, tumor-only best-fit
selection, fitted-region calculations, and consolidated copy number with BAF
inference. It also reconstructs the summary/QC context, including canonical-gene
deletion counts from the Ensembl cache, LOH, polyclonality, and WGD, and writes
the six gated core outputs when `-output_dir` is supplied.
