# Recurrence upper limits

`recurrence_upper_limits.tsv` supplements the catalog of flagged taxa. It uses
the same 5,865 molecular TimeTree study trees, 5,072 source publications, and NCBI
Taxonomy downloaded on 22 July 2026. Of 26,291 represented taxa, 18,109 have
informative publications and calculated limits; 8,182 have no informative
publications and missing limits. The flag column identifies the same 516 taxa
as the flagged-taxon catalog. No additional threshold for selecting supported
taxa is imposed.

## Interpretation and calculation

These values are approximate upper endpoints of individual conditional 95%
confidence sets for the recurring-contradiction parameter, rho. They are not
probabilities that taxa are monophyletic, proportions of misplaced species, or
raw fractions of contradictory publications. A low upper limit means that
strong recurrence is inconsistent with the available evidence under the model.
An upper limit of one means the evidence does not exclude the maximum recurrence
allowed by the model; it does not estimate recurrence as one.

The calculation holds the fitted background and detection probabilities fixed.
For each proposed rho, a likelihood-ratio statistic compares it with the best
fit in [0,1]. The null distribution is calculated by enumerating all binary
publication outcomes for up to 12 publications, or by Monte Carlo for larger
records. Tests are inverted to obtain the confidence set. This uses the general
likelihood-ratio ordering principle of Feldman and Cousins (1998), adapted to
heterogeneous Bernoulli outcomes rather than their Poisson model:
<https://doi.org/10.1103/PhysRevD.57.3873>.

The calculation does not propagate uncertainty from estimating background and
detection probabilities or account for unknown dependence among publications.
The 95% level applies to each conditional confidence set, not simultaneously
to all taxa. The upper-limit MLE is optimized continuously by bisection, whereas
the flagging analysis uses a fixed grid, so their fitted rho values can differ
slightly.

## Numerical settings

The released calculation uses 4,000 Monte Carlo replicates per candidate rho,
increased to 20,000 when the initial P value lies within .015 of .05. Random
streams depend on taxid and candidate rho, with base seed 20261105. Candidate
values span [0,1] at .025 spacing, also including the fitted value. The highest
accepted grid value is refined at .0025 spacing. Finite endpoints receive .0025
upward padding. These are finite-grid and Monte Carlo approximations, not
certified continuous-domain bounds. Do not overinterpret the final decimal
places. When rho=1 is accepted, the upper endpoint is one without further search.

## Columns

| Column | Meaning |
|---|---|
| `taxid`, `name`, `rank` | Identity in the frozen NCBI taxonomy |
| `informative_publications` | Number of informative publication observations |
| `contradicted_publications` | Observations contradicting the group's monophyly |
| `monophyletic_publications` | Observations displaying the group as a clade |
| `rho_mle` | Best-fitting recurrence in [0,1] |
| `upper_recurrence_limit_95` | Approximate calibrated confidence-set upper endpoint |
| `mean_background_probability` | Mean fitted background across informative publications |
| `mean_detection_probability` | Mean fitted detection across informative publications |
| `max_sampled_members` | Maximum of the publication-level mean sampled-member counts |
| `calibration` | Exact enumeration or parametric Monte Carlo |
| `largest_MC_sample` | Largest simulation count used at any candidate rho; zero if none was needed |
| `evaluated_rho_values` | Number of candidate values evaluated |
| `endpoint_grid_padding` | Upward padding for a finite endpoint |
| `evidence_status` | Informative observations or no informative publications |
| `currently_flagged_recurrent_nonmonophyly` | Membership in the 516-taxon flagged catalog |

For taxa without informative observations, publication counts are zero and
fitted recurrence and limits are missing. No observation is not evidence of
monophyly. All rows of this released resource are complete.

## Reproduction

Use `taxaaudit upper-limits` as described in the main README. The released table
is unchanged from the completed analysis; its input hashes and settings are in
`upper_limits_provenance.json`. Input study trees and publication-level records
are not redistributed here. The software can fit those probabilities for a
user-supplied collection with `taxaaudit run`.
