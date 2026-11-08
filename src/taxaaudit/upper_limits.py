"""Conditional confidence-set upper endpoints for contradiction recurrence.

The likelihood-ratio ordering and numerical settings match the published-data
calculation. Background and detection probabilities are treated as fixed.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def log_likelihood(y, f, a, rho):
    p = np.clip(f + a * np.asarray(rho), 1e-12, 1 - 1e-12)
    return np.sum(y * np.log(p) + (1 - y) * np.log1p(-p), axis=-1)


def fit_recurrence(y, f, a):
    """Fit one or more binary records, including the parameter boundaries."""
    y = np.atleast_2d(y)
    lo = np.zeros((len(y), 1))
    hi = np.ones_like(lo)
    for _ in range(25):
        mid = (lo + hi) / 2
        p = np.clip(f + a * mid, 1e-12, 1 - 1e-12)
        score = np.sum(a * (y - p) / (p * (1 - p)), axis=1, keepdims=True)
        lo = np.where(score > 0, mid, lo)
        hi = np.where(score > 0, hi, mid)
    rho = (lo + hi) / 2
    best = log_likelihood(y, f, a, rho)
    zero = log_likelihood(y, f, a, 0)
    one = log_likelihood(y, f, a, 1)
    rho = np.where((zero >= best)[:, None], 0, rho)
    best = np.maximum(best, zero)
    rho = np.where((one >= best)[:, None], 1, rho)
    return rho[:, 0], np.maximum(best, one)


def binary_records(n):
    return ((np.arange(2**n, dtype=np.uint64)[:, None] >>
             np.arange(n, dtype=np.uint64)) & 1).astype(float)


class RecurrenceCalibrator:
    """Invert 95% likelihood-ratio tests for a single taxon's publications."""

    def __init__(self, y, f, d, seed, replicates=4000):
        self.y, self.f, self.d = (np.asarray(v, dtype=float) for v in (y, f, d))
        if (self.y.ndim != 1 or not len(self.y) or
                self.y.shape != self.f.shape or self.y.shape != self.d.shape):
            raise ValueError("y, background and detection must be equal-length nonempty vectors")
        if not np.isin(self.y, [0, 1]).all():
            raise ValueError("outcomes must be zero or one")
        if any(not np.isfinite(v).all() or ((v < 0) | (v > 1)).any()
               for v in (self.f, self.d)):
            raise ValueError("background and detection must be finite probabilities")
        if seed < 0 or replicates < 1:
            raise ValueError("seed must be nonnegative and replicates must be positive")
        self.a = (1 - self.f) * self.d
        self.seed, self.replicates = seed, replicates
        self.fitted, self.best = fit_recurrence(self.y, self.f, self.a)
        self.exact = len(self.y) <= 12
        if self.exact:
            self.sim = binary_records(len(self.y))
            _, self.sim_best = fit_recurrence(self.sim, self.f, self.a)
        self.cache = {}

    def pvalue(self, rho):
        if not 0 <= rho <= 1:
            raise ValueError("rho must lie in [0, 1]")
        if rho in self.cache:
            return self.cache[rho]
        observed = max(0, float(2 * (self.best[0] -
                       log_likelihood(self.y, self.f, self.a, rho))))
        if observed <= 1e-9:
            self.cache[rho] = (1.0, 0)
            return self.cache[rho]
        if self.exact:
            null = log_likelihood(self.sim, self.f, self.a, rho)
            stat = np.maximum(0, 2 * (self.sim_best - null))
            weights = np.exp(null)
            if not np.isclose(weights.sum(), 1, atol=1e-8, rtol=0):
                raise ValueError("null probabilities do not sum to one")
            p, count = float(weights[stat >= observed - 1e-9].sum()), 0
        else:
            # Independent streams for each taxon and proposed rho make resume
            # and the order of test evaluation immaterial.
            rng = np.random.default_rng(np.random.SeedSequence(
                [self.seed, int(round(rho * 1e8))]))
            count = self.replicates
            sim = rng.binomial(1, self.f + self.a * rho, size=(count, len(self.y)))
            _, best = fit_recurrence(sim, self.f, self.a)
            stat = np.maximum(0, 2 * (best - log_likelihood(sim, self.f, self.a, rho)))
            exceed = int(np.count_nonzero(stat >= observed - 1e-9))
            if abs((exceed + 1) / (count + 1) - .05) < .015:
                extra = 4 * self.replicates
                sim = rng.binomial(1, self.f + self.a * rho, size=(extra, len(self.y)))
                _, best = fit_recurrence(sim, self.f, self.a)
                stat = np.maximum(0, 2 * (best - log_likelihood(sim, self.f, self.a, rho)))
                exceed += int(np.count_nonzero(stat >= observed - 1e-9))
                count += extra
            p = (exceed + 1) / (count + 1)
        self.cache[rho] = (p, count)
        return self.cache[rho]

    def upper(self):
        if self.pvalue(1.0)[0] > .05:
            return 1.0
        grid = np.unique(np.r_[np.linspace(0, 1, 41), self.fitted[0]])
        values = np.array([self.pvalue(float(r))[0] for r in grid])
        top = float(grid[values > .05].max())
        fine = np.unique(np.r_[np.arange(top, min(1, top + .025) + .0025 / 2, .0025), top])
        accepted = [float(r) for r in fine if r <= 1 and self.pvalue(float(r))[0] > .05]
        return min(1, max(accepted) + .0025)


def calculate_upper_limits(observations, universe, catalog, output, *, seed=20261105,
                           replicates=4000):
    """Calculate all represented taxa; resume only with unchanged inputs/settings."""
    paths = {k: Path(v) for k, v in dict(
        observations=observations, universe=universe, catalog=catalog).items()}
    frame = pd.read_csv(paths["observations"], sep="\t")
    if "publication" not in frame and "source" in frame:
        frame = frame.rename(columns={"source": "publication"})
    required = {"taxid", "name", "rank", "publication", "observed_failure",
                "background_apparent_failure_probability", "cross_study_detection_probability",
                "mean_sampled_members"}
    if not required <= set(frame):
        raise ValueError(f"missing observation columns: {sorted(required - set(frame))}")
    if frame.duplicated(["taxid", "publication"]).any():
        raise ValueError("each taxon must have at most one observation per publication")
    if not frame.observed_failure.isin([0, 1]).all():
        raise ValueError("outcomes must be zero or one")
    for column in ["background_apparent_failure_probability", "cross_study_detection_probability"]:
        if not np.isfinite(frame[column]).all() or not frame[column].between(0, 1).all():
            raise ValueError(f"invalid probabilities in {column}")
    if "study_weight" in frame:
        weights = frame.study_weight
        if not np.isfinite(weights).all() or (weights <= 0).any():
            raise ValueError("study weights must be finite and positive")
        if frame.groupby("taxid").study_weight.nunique().gt(1).any():
            raise ValueError("upper limits for nonuniform publication weights are not implemented")
    if seed < 0 or replicates < 1:
        raise ValueError("seed must be nonnegative and replicates must be positive")
    taxa = pd.read_csv(paths["universe"], sep="\t")[["taxid", "name", "rank"]]
    flagged = set(pd.read_csv(paths["catalog"], sep="\t").taxid)
    if taxa.taxid.duplicated().any() or not set(frame.taxid) <= set(taxa.taxid):
        raise ValueError("universe must have unique taxids and include all observed taxa")
    if not flagged <= set(taxa.taxid):
        raise ValueError("flagged catalog contains taxa absent from the universe")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    config = {"input_sha256": {k: hashlib.sha256(p.read_bytes()).hexdigest()
                               for k, p in paths.items()},
              "seed": seed, "replicates": replicates, "method": "conditional_lr_upper95_v1"}
    manifest = output / "upper_limits_settings.json"
    checkpoint = output / "upper_limits_checkpoint.jsonl"
    if manifest.exists():
        if json.loads(manifest.read_text()) != config:
            raise ValueError("upper-limit inputs/settings changed; choose a new output directory")
    elif checkpoint.exists():
        raise ValueError("checkpoint has no settings file; choose a new output directory")
    else:
        manifest.write_text(json.dumps(config, indent=2) + "\n")
    done = {}
    if checkpoint.exists():
        for line in checkpoint.read_text().splitlines():
            row = json.loads(line)
            done[row["taxid"]] = row
    with checkpoint.open("a", buffering=1) as handle:
        for taxid, group in frame.groupby("taxid", sort=True):
            if int(taxid) in done:
                continue
            y = group.observed_failure.to_numpy(float)
            f = group.background_apparent_failure_probability.to_numpy(float)
            d = group.cross_study_detection_probability.to_numpy(float)
            fit = RecurrenceCalibrator(y, f, d, seed + int(taxid), replicates)
            upper = fit.upper()
            row = dict(taxid=int(taxid), informative_publications=len(group),
                       contradicted_publications=int(y.sum()),
                       monophyletic_publications=int(len(group) - y.sum()),
                       rho_mle=float(fit.fitted[0]), upper_recurrence_limit_95=upper,
                       mean_background_probability=float(f.mean()),
                       mean_detection_probability=float(d.mean()),
                       max_sampled_members=float(group.mean_sampled_members.max()),
                       calibration="exact_enumeration" if fit.exact else "parametric_Monte_Carlo",
                       largest_MC_sample=max(v[1] for v in fit.cache.values()),
                       evaluated_rho_values=len(fit.cache),
                       endpoint_grid_padding=.0025 if upper < 1 else 0.0,
                       evidence_status="informative_observations")
            handle.write(json.dumps(row) + "\n")
            done[int(taxid)] = row
    columns = ["taxid", "informative_publications", "contradicted_publications",
               "monophyletic_publications", "rho_mle", "upper_recurrence_limit_95",
               "mean_background_probability", "mean_detection_probability", "max_sampled_members",
               "calibration", "largest_MC_sample", "evaluated_rho_values",
               "endpoint_grid_padding", "evidence_status"]
    resource = taxa.merge(pd.DataFrame(done.values(), columns=columns), on="taxid", how="left",
                          validate="one_to_one")
    missing = resource.evidence_status.isna()
    resource.loc[missing, "evidence_status"] = "no_informative_publications"
    counts = ["informative_publications", "contradicted_publications", "monophyletic_publications"]
    resource.loc[missing, counts] = 0
    for column in counts + ["largest_MC_sample", "evaluated_rho_values"]:
        resource[column] = resource[column].astype("Int64")
    resource["currently_flagged_recurrent_nonmonophyly"] = resource.taxid.isin(flagged)
    resource.to_csv(output / "taxaaudit_recurrence_upper_limits.tsv", sep="\t", index=False,
                    float_format="%.8g")
    summary = dict(total_taxa=len(resource), taxa_with_bounds=len(done),
                   taxa_without_informative_publications=int(missing.sum()),
                   flagged_taxa=len(flagged), confidence_level=.95,
                   conditional_on="fixed fitted background/detection and independent publications")
    (output / "upper_limits_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return resource
