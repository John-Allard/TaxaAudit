# Catalog of flagged NCBI taxa

This directory contains the result of applying TaxaAudit to 5,865
molecular TimeTree study-tree files representing 5,072 publications, using the
NCBI Taxonomy downloaded July 22, 2026.

- `supported_original_constraint_deletions.tsv` contains the 516 flagged
  NCBI groups and their evidence summaries.
- `constraint_changes.json` contains the same flagged taxids in a compact
  machine-readable form.
- `summary.json` records the tested universe, selected publication-count
  filter, and number of findings.
- `independent_filter_profile.tsv` reports results across candidate minimum
  publication counts.
- `provenance.json` records input vintage, model settings, and the SHA-256 of
  the released flagged-taxon table.
- `recurrence_upper_limits.tsv` contains all 26,291 represented taxa, with
  conditional 95% recurrence upper limits for the 18,109 having informative
  publications. See [column definitions and interpretation](UPPER_LIMITS.md).
- `upper_limits_provenance.json` records the upper-limit calculation settings,
  input hashes, and the released table's SHA-256.

The catalog identifies taxa with statistical evidence of recurrent contradiction
of monophyly. It is not a replacement taxonomy. Excluding a flagged taxon's exact
monophyly requirement does not remove its species or change their names. The
input study trees are not redistributed in this repository.

Filenames and machine-readable field names retain earlier terminology such as
`deletions`, `selected_groups`, and `deletion_only` for compatibility with existing
analyses. These refer to exclusion of monophyly requirements, not deletion of
species. The data files are unchanged from the original analysis.
