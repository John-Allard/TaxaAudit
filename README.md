# TaxaAudit

**TaxaAudit: statistical identification of taxa recurrently contradicted by
published molecular phylogenies**

TaxaAudit identifies NCBI taxonomic groups whose sampled members are
recurrently contradicted as a clade across published molecular
phylogenies. It does **not** construct a replacement taxonomy. Its primary
output is a catalog of flagged taxa: NCBI groups with statistical evidence of
recurrent contradiction of monophyly.

The statistical model separates three quantities for NCBI group `g` in
publication `i`:

```text
P(contradiction) = f_ig + (1 - f_ig) d_ig rho_g
```

- `f_ig` is a cross-fitted estimate of ordinary apparent conflict expected
  from the publication, rank, and focal-group sampling.
- `d_ig` is a cross-fitted estimate of whether the publication's species
  sampling could expose a real violation.
- `rho_g` is the recurring group-specific contradiction tendency fitted across
  publications.

The repository contains the complete fitting workflow. Fitted machine-learning
objects are intentionally not committed: the background and detection models
are cross-fitted for the supplied study collection and are reproducibly rebuilt
from the input trees.

## Installation

Python 3.10 or later is recommended.

```bash
git clone https://github.com/John-Allard/TaxaAudit.git
cd TaxaAudit
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

Required Python packages are Biopython, NumPy, pandas, SciPy, and
scikit-learn.

## Run on a directory of Newick trees

Download and unpack an NCBI `taxdump.tar.gz`, or let the command do it:

```bash
taxaaudit download-taxonomy --out data/ncbi_taxonomy
```

Then run the complete workflow:

```bash
taxaaudit run \
  --trees /path/to/newick_trees \
  --taxonomy data/ncbi_taxonomy \
  --out outputs/my_scan
```

Recognized filename extensions are `.nwk`, `.newick`, `.tree`, and `.tre`.
By default, every file is treated as a separate publication. If a publication
contributed more than one tree, supply a tab-separated manifest so that it
still contributes only one result per taxonomic group:

```text
path	publication
/absolute/path/study1_tree1.nwk	PMID_123456
/absolute/path/study1_tree2.nwk	PMID_123456
/absolute/path/study2.nwk	PMID_987654
```

```bash
taxaaudit run \
  --manifest trees.tsv \
  --taxonomy data/ncbi_taxonomy \
  --out outputs/my_scan
```

Terminal labels may be bare NCBI taxids (`9606`), labels beginning
`tx<TAXID>_` (`tx9606_Homo_sapiens`), or unambiguous NCBI scientific names
with spaces or underscores. Unmapped or ambiguous tips are reported and
omitted.

The most useful outputs are:

```text
catalog/supported_original_constraint_deletions.tsv
catalog/constraint_changes.json
catalog/all_group_evidence.tsv
catalog/publication_evidence.tsv
```

The first file is the catalog of flagged taxa. Its filename retains the earlier
term `deletions` for compatibility. Excluding a flagged taxon's monophyly
requirement does **not** remove any species from a tree or dataset. See
[`catalog/README.md`](catalog/README.md) for the bundled results and file formats.

## Upper confidence limits on recurrence

After a run, calculate an upper confidence limit for every taxon with informative
publications:

```bash
taxaaudit upper-limits --run outputs/my_scan --out outputs/my_scan/upper_limits
```

This uses the fitted background and detection probabilities saved by `run` in
`catalog/publication_evidence.tsv`. It does not retrain the models or change the
catalog of flagged taxa. It writes `taxaaudit_recurrence_upper_limits.tsv`, a
summary, a settings file, and a checkpoint. Repeat the same command to resume
an interrupted calculation; changed inputs or settings require a new output
directory. Groups with no informative publications remain in the table with
missing limits, not zero limits.

For an existing collection of fitted observations, provide the inputs directly:

```bash
taxaaudit upper-limits \
  --observations publication_evidence.tsv \
  --universe all_group_evidence.tsv \
  --catalog supported_original_constraint_deletions.tsv \
  --out outputs/upper_limits
```

The observation table needs `taxid`, `name`, `rank`, `publication` (or `source`),
`observed_failure`, `background_apparent_failure_probability`,
`cross_study_detection_probability`, and `mean_sampled_members`. Each row is
one informative observation of a taxon in a publication. The universe table
needs `taxid`, `name`, and `rank`, and the flagged catalog needs `taxid`.
The first two probability columns are fitted values, not raw conflict fractions.
Older runs without `publication_evidence.tsv` need their saved fitted observations
or a new run. Upper limits with nonuniform publication weights are not supported.

The limits are approximate upper endpoints of **individual conditional 95%
confidence sets** for the recurrence parameter. They condition on the fitted
background/detection probabilities and independent publication outcomes. They
are not probabilities of monophyly or a simultaneous confidence guarantee across
all taxa. Low limits distinguish well-sampled groups with little evidence of
recurrent contradiction from groups with too little evidence to constrain it.
See the [released upper-limit data](catalog/UPPER_LIMITS.md) for interpretation,
numerical settings, and column definitions.

## Optional study weights and higher ranks

The default analysis uses no external metadata-based study weights. Its background
and detection probabilities already account for evidence within the supplied trees.
For an exploratory weighted analysis, add `--study-weights weights.tsv`, with
tab-separated `publication` and `weight` columns. Publication identifiers must
match the manifest. Values must be finite and positive; unlisted publications
receive weight one. Relative weights multiply publication log-likelihoods after
normalization to mean one within each group. The null simulations use the same
weights, and weighted runs bootstrap all observed groups rather than using an
unweighted asymptotic screen. Weights neither add replication nor change the
background and detection fits. Prespecify weights independently of focal outcomes.
This option is unit-tested but has not received the released analysis's full
end-to-end validation and was not used for the released catalog.

To scan explicit higher ranks, supply for example
`--ranks superfamily,order,class,phylum`. The default remains family and below.
Rank-specific models and calibration should be checked before treating higher-rank
results as a validated extension.

## Reproduce the released NCBI scan

The released catalog in [`catalog/`](catalog/) was generated from:

- 5,865 molecular study-tree files representing 5,072 source publications in
  the TimeTree literature collection;
- the NCBI Taxonomy downloaded **July 22, 2026**;
- NCBI taxdump SHA-256
  `9706198c129c6b70ec591009431f742649b939d5de86e102d469e6928034a1fe`;
- tests of named NCBI groups at the family level and below, including named
  unranked clades but excluding administrative placeholder groups;
- 2,000 contradicted complete-tree/group templates, stratified by rank and
  tree size, with 40 reductions per template;
- five taxonomic-region folds for production detection probabilities;
- five group-held-out folds for the background sampling model;
- a 201-point grid for `rho` from 0 through 1;
- 99,999 parametric-bootstrap replicates for groups passing the broad
  asymptotic screen (`q <= 0.30`);
- a minimum of five informative publications, selected by the publication-count
  filtering procedure;
- Benjamini-Hochberg-adjusted `q <= 0.05` as the flagging threshold.

Of 26,291 represented NCBI groups, 5,395 passed the publication-count filter
and 516 were flagged. Of the represented groups, 18,109 had at least one
informative publication. The released table includes exact NCBI identifiers,
names, ranks, evidence counts, fitted recurrence, P values, and q values.

The bundled [`catalog/recurrence_upper_limits.tsv`](catalog/recurrence_upper_limits.tsv)
adds the upper-limit results for the same 26,291 taxa: 18,109 have informative
publications and calculated limits, while 8,182 are explicitly marked as lacking
informative publications. It includes membership in the unchanged 516-taxon
flagged catalog.

The study trees themselves are not redistributed. Users with those trees can
reproduce the analysis using:

```bash
taxaaudit run \
  --manifest timetree_studies.tsv \
  --taxonomy taxonomy_cache_20260722 \
  --max-templates 2000 \
  --draws-per-template 40 \
  --detection-folds 5 \
  --background-folds 5 \
  --rho-grid-size 201 \
  --bootstrap-replicates 99999 \
  --bootstrap-screen-q 0.30 \
  --filter-min 2 \
  --filter-max 20 \
  --fdr-q 0.05 \
  --detection-seed 20260821 \
  --background-seed 20260728 \
  --bootstrap-seed 20260728 \
  --out reproduced_catalog
```

## Interpretation

A flagged taxon has evidence of recurrent contradiction of monophyly in the
sampled molecular phylogenetics literature. Flagging supports excluding that
taxon's exact monophyly requirement in applications that rely on it. It does
not identify a unique replacement group, prove a particular biological mechanism, or
establish that unflagged groups are monophyletic. Groups with insufficient
published sampling are evidence-limited rather than confirmed.

## Paper

The accompanying manuscript is in preparation. A link to the paper and its
supplementary methods will be added here when a public version is available.

## Development

Run the automated tests from the repository root:

```bash
pip install pytest
python -m pytest
```
