# PURPLE v4.4 tumor-only C++ port

This directory is a checkpoint-driven port of the tumor-only WGS core used by
PURPLE v4.4.  The frozen behavior contract and Java reference artifacts live in
`research/studies/purple-port-fidelity-v1/runs/RUN-P001`.

The current executable implements the CP-P1 input boundary: AMBER BAF/QC/PCF,
COBALT ratio/PCF, tumor-only gender adjustment, and Java-compatible PCF boundary
construction.  Later checkpoints add support segmentation, fitting, copy number,
QC, and the six gated output writers.
