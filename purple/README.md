# PURPLE v4.4 tumor-only C++ port

This directory is a checkpoint-driven port of the tumor-only WGS core used by
PURPLE v4.4.  The frozen behavior contract and Java reference artifacts live in
`research/studies/purple-port-fidelity-v1/runs/RUN-P001`.

The current executable implements CP-P1 through CP-P7: input parsing and
tumor-only adjustment, support segmentation, observed-region aggregation,
fitting-region selection, the full purity/ploidy score grid, tumor-only best-fit
selection, and fitted-region calculations. Later checkpoints add consolidated
copy number, QC, and the six gated writers.
