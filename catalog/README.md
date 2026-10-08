# Released NCBI catalog

This directory contains the result of applying TaxaAudit to 5,865
molecular TimeTree study-tree files representing 5,072 publications, using the
NCBI Taxonomy downloaded July 22, 2026.

- `supported_original_constraint_deletions.tsv` contains the 516 selected
  NCBI groups and their evidence summaries.
- `constraint_changes.json` contains the same selected taxids in a compact
  machine-readable form.
- `summary.json` records the tested universe, selected publication-count
  filter, and number of findings.
- `independent_filter_profile.tsv` reports results across candidate minimum
  publication counts.
- `provenance.json` records input vintage, model settings, and the SHA-256 of
  the released selected-group table.

The catalog identifies exact clade assertions that are recurrently
contradicted. It is not a replacement taxonomy. The input study trees are not
redistributed in this repository.
