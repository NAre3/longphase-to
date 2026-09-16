# PURPLE v4.4 tumor-only C++ port

This directory is a checkpoint-driven port of the tumor-only WGS core used by
PURPLE v4.4.  The frozen behavior contract and Java reference artifacts live in
`research/studies/purple-port-fidelity-v1/runs/RUN-P001`.

The current executable implements CP-P1 through CP-P9: input parsing and
tumor-only adjustment, support segmentation, observed-region aggregation,
fitting-region selection, the full purity/ploidy score grid, tumor-only best-fit
selection, fitted-region calculations, and consolidated copy number with BAF
inference. It also reconstructs the summary/QC context, including canonical-gene
deletion counts from the Ensembl cache, LOH, polyclonality, and WGD. The final
checkpoint adds the six gated writers.
