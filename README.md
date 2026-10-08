# TaxaAudit

**TaxaAudit: statistical identification of taxa recurrently contradicted by
published molecular phylogenies**

TaxaAudit identifies NCBI taxonomic groups whose sampled members are
recurrently contradicted as a clade across independent published molecular
phylogenies. It does **not** construct a replacement taxonomy. Its primary
output is a conservative list of exact NCBI taxids for which monophyly should
not be assumed without further evidence.

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
still contributes only one independent result per taxonomic group:

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
```

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
- an independently selected minimum of five informative publications;
- Benjamini-Hochberg false-discovery-rate control at `q <= 0.05`.

Of 26,291 represented NCBI groups, 5,395 passed the publication-count filter
and 516 were selected. The released table includes exact NCBI identifiers,
names, ranks, evidence counts, fitted recurrence, P values, and q values.

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

A selected taxid has reproducible evidence that its sampled members do not
form an exclusive clade in the molecular phylogenetics literature. Selection
supports withdrawing that exact monophyly assertion. It does not identify a
unique replacement group, prove a particular biological mechanism, or
establish that unselected groups are monophyletic. Groups with insufficient
published sampling are evidence-limited rather than confirmed.

See [`METHODS.md`](METHODS.md) for the algorithm and
[`catalog/README.md`](catalog/README.md) for the released output.
