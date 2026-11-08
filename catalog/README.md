# Catalog of flagged NCBI taxa

This directory contains the result of applying TaxaAudit to 5,865
molecular TimeTree study-tree files representing 5,072 publications, using the
NCBI Taxonomy downloaded 22 July 2026.

- `supported_original_constraint_deletions.tsv` contains the 516 flagged
  NCBI groups and their evidence summaries.
- `constraint_changes.json` contains the same flagged taxids in a compact
  machine-readable form.
- `summary.json` records the tested universe, selected publication-count
  filter, and number of findings.
- `independent_filter_profile.tsv` reports results across candidate minimum
  publication counts.
- `provenance.json` records the input collection, taxonomy snapshot, and model
  settings.
- `recurrence_upper_limits.tsv` contains all 26,291 represented taxa, with
  conditional 95% recurrence upper limits for the 18,109 having informative
  publications. See [column definitions and interpretation](UPPER_LIMITS.md).
- `upper_limits_provenance.json` records the upper-limit calculation settings
  and input collection.

The catalog identifies taxa with statistical evidence of recurrent contradiction
of monophyly. It is not a replacement taxonomy. Targeted exclusion means omitting
a flagged taxon's exact monophyly requirement while retaining all species and
the monophyly requirements of unflagged ancestors and descendants.
The study trees are available from [TimeTree](https://timetree.org/).

Filenames and machine-readable field names retain earlier terminology such as
`deletions`, `selected_groups`, and `deletion_only` for compatibility with existing
analyses. These refer to exclusion of monophyly requirements, not deletion of
species. The data files are unchanged from the original analysis.
