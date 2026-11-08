from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import chi2
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold


DETECTION_FEATURES = [
    "sampled_members",
    "sampled_parent_nonmembers",
    "represented_sibling_lineages",
]


def stable_seed(base: int, value: object) -> int:
    digest = hashlib.sha256(f"{base}:{value}".encode()).digest()
    return (base + int.from_bytes(digest[:8], "big")) % (2**32)


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="stable")
    ordered = values[order]
    adjusted = np.minimum.accumulate(
        (ordered * len(ordered) / np.arange(1, len(ordered) + 1))[::-1]
    )[::-1]
    result = np.empty(len(values), dtype=float)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def _template_weights(frame: pd.DataFrame) -> np.ndarray:
    sizes = frame.groupby("template_id").template_id.transform("size").to_numpy(float)
    weights = 1.0 / sizes
    return weights * (len(weights) / weights.sum())


def _detection_model(seed: int) -> HistGradientBoostingClassifier:
    return HistGradientBoostingClassifier(
        loss="log_loss", learning_rate=0.05, max_iter=150,
        max_leaf_nodes=15, min_samples_leaf=40, l2_regularization=1.0,
        monotonic_cst=[1, 1, 1], random_state=seed,
    )


def fit_detection_probabilities(
    reductions: pd.DataFrame,
    tree_observations: pd.DataFrame,
    publication_observations: pd.DataFrame,
    folds: int,
    seed: int,
    output: str | Path,
) -> pd.DataFrame:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    x = np.log1p(reductions[DETECTION_FEATURES].to_numpy(float))
    y = reductions.detected_contradiction.to_numpy(int)
    weights = _template_weights(reductions)
    regions = reductions.region_taxid.astype(str).to_numpy()
    splitter = GroupKFold(n_splits=min(folds, len(np.unique(regions))))
    region_fold: dict[str, int] = {}
    models: dict[int, HistGradientBoostingClassifier] = {}
    validation_prediction = np.full(len(reductions), np.nan)
    for fold, (train, test) in enumerate(splitter.split(x, y, regions)):
        model = _detection_model(seed + fold + 1)
        model.fit(x[train], y[train], sample_weight=weights[train])
        models[fold] = model
        validation_prediction[test] = model.predict_proba(x[test])[:, 1]
        for region in np.unique(regions[test]):
            region_fold[str(region)] = fold
    final = _detection_model(seed)
    final.fit(x, y, sample_weight=weights)
    feature_x = np.log1p(tree_observations[DETECTION_FEATURES].to_numpy(float))
    prediction = final.predict_proba(feature_x)[:, 1]
    sources = np.full(len(tree_observations), "final_model_unseen_region", dtype=object)
    feature_regions = tree_observations.region_taxid.astype(str).to_numpy()
    for fold, model in models.items():
        selected = np.array([region_fold.get(region) == fold for region in feature_regions])
        if selected.any():
            prediction[selected] = model.predict_proba(feature_x[selected])[:, 1]
            sources[selected] = "region_cross_fitted"
    tree_predictions = tree_observations.copy()
    tree_predictions["cross_study_detection_probability"] = np.clip(prediction, 0, 1)
    tree_predictions["detection_prediction_source"] = sources
    collapsed = (
        tree_predictions.groupby(["taxid", "publication"], as_index=False)
        .agg(
            cross_study_detection_probability=(
                "cross_study_detection_probability", "mean"
            ),
            region_cross_fitted_fraction=(
                "detection_prediction_source",
                lambda values: float(np.mean(values == "region_cross_fitted")),
            ),
        )
    )
    publication_observations = publication_observations.drop(
        columns=["cross_study_detection_probability"], errors="ignore"
    )
    result = publication_observations.merge(
        collapsed, on=["taxid", "publication"], how="left", validate="one_to_one"
    )
    if result.cross_study_detection_probability.isna().any():
        raise RuntimeError("some publication observations lack detection predictions")
    result.to_csv(output / "publication_group_detection_predictions.tsv", sep="\t", index=False)
    validation = reductions[["template_id", "publication", "region_taxid", "detected_contradiction"]].copy()
    validation["cross_fitted_prediction"] = validation_prediction
    validation.to_csv(output / "detection_cross_fitted_validation.tsv", sep="\t", index=False)
    summary = {
        "training_templates": int(reductions.template_id.nunique()),
        "training_reductions": int(len(reductions)),
        "folds": int(folds),
        "features": DETECTION_FEATURES,
        "cross_fitted_brier": float(np.mean((y - validation_prediction) ** 2)),
        "production_observations": int(len(result)),
        "production_prediction_mean": float(result.cross_study_detection_probability.mean()),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return result


def add_background(frame: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    cells = (
        frame.groupby(["publication", "rank"], sort=False).observed_failure
        .agg(cell_failures="sum", cell_observations="count").reset_index()
    )
    ranks = (
        frame.groupby("rank", sort=False).observed_failure
        .agg(rank_failures="sum", rank_observations="count").reset_index()
    )
    focal = (
        frame.groupby("taxid", sort=False).observed_failure
        .agg(focal_failures="sum", focal_observations="count").reset_index()
    )
    out = frame.merge(cells, on=["publication", "rank"], validate="many_to_one")
    out = out.merge(ranks, on="rank", validate="many_to_one")
    out = out.merge(focal, on="taxid", validate="many_to_one")
    y = out.observed_failure.to_numpy(float)
    rank_mean = (
        out.rank_failures.to_numpy(float) - out.focal_failures.to_numpy(float) + 0.5
    ) / (
        out.rank_observations.to_numpy(float) - out.focal_observations.to_numpy(float) + 1
    )
    cell_failures = out.cell_failures.to_numpy(float) - y
    cell_n = out.cell_observations.to_numpy(float) - 1

    def probabilities(strength: float) -> np.ndarray:
        return np.clip(
            (cell_failures + strength * rank_mean) / (cell_n + strength),
            1e-8, 1 - 1e-8,
        )

    def objective(log_strength: float) -> float:
        p = probabilities(math.exp(log_strength))
        return float(-np.sum(y * np.log(p) + (1 - y) * np.log1p(-p)))

    fit = minimize_scalar(
        objective, bounds=(math.log(0.05), math.log(500)), method="bounded"
    )
    strength = float(math.exp(fit.x))
    out["background_apparent_failure_probability"] = probabilities(strength)
    return out, strength


def adjust_background_for_sampling(
    frame: pd.DataFrame, folds: int, seed: int
) -> tuple[pd.DataFrame, dict[str, float]]:
    out = frame.copy()
    out["empirical_bayes_background_apparent_failure_probability"] = out[
        "background_apparent_failure_probability"
    ]
    base = np.clip(out.background_apparent_failure_probability.to_numpy(float), 1e-8, 1-1e-8)
    x = np.column_stack([
        np.log(base / (1 - base)),
        np.log1p(out.mean_sampled_members.to_numpy(float)),
    ])
    y = out.observed_failure.to_numpy(int)
    groups = out.taxid.astype(str).to_numpy()
    splitter = GroupKFold(n_splits=min(folds, len(np.unique(groups))))
    predictions = np.full(len(out), np.nan)
    for fold, (train, test) in enumerate(splitter.split(x, y, groups)):
        model = HistGradientBoostingClassifier(
            loss="log_loss", learning_rate=0.05, max_iter=200,
            max_leaf_nodes=15, min_samples_leaf=100, l2_regularization=2.0,
            monotonic_cst=[1, 1], random_state=stable_seed(seed, fold),
        )
        model.fit(x[train], y[train])
        predictions[test] = model.predict_proba(x[test])[:, 1]
    predictions = np.clip(predictions, 1e-8, 1-1e-8)
    out["background_apparent_failure_probability"] = predictions
    return out, {
        "cross_fitted_background_brier": float(np.mean((y - predictions) ** 2)),
        "cross_fitted_background_log_loss": float(
            -np.mean(y*np.log(predictions) + (1-y)*np.log1p(-predictions))
        ),
    }


def fit_group(frame: pd.DataFrame, rho_grid: np.ndarray):
    y = frame.observed_failure.to_numpy(float)
    f = frame.background_apparent_failure_probability.to_numpy(float)
    d = frame.cross_study_detection_probability.to_numpy(float)
    probability = np.clip(
        f[:, None] + (1-f[:, None]) * d[:, None] * rho_grid[None, :],
        1e-12, 1-1e-12,
    )
    null = np.clip(f[:, None], 1e-12, 1-1e-12)
    terms = y[:, None]*(np.log(probability)-np.log(null)) + (1-y[:, None])*(
        np.log1p(-probability)-np.log1p(-null)
    )
    weights = study_weights(frame)
    log_lr = (terms * weights[:, None]).sum(axis=0)
    best = int(np.argmax(log_lr))
    statistic = max(0.0, 2*float(log_lr[best]))
    first = frame.iloc[0]
    return {
        "taxid": int(first.taxid), "name": first["name"], "rank": first["rank"],
        "informative_publications": len(frame),
        "violating_publications": int(y.sum()),
        "clean_publications": int((1-y).sum()),
        "fitted_reproducible_nonmonophyly_tendency": float(rho_grid[best]),
        "likelihood_ratio_statistic": statistic,
        "one_sided_asymptotic_p": float(0.5*chi2.sf(statistic, 1)) if statistic > 1e-12 else 1.0,
        "mean_background_apparent_failure_probability": float(f.mean()),
        "mean_detection_power_violating": float(d[y == 1].mean()) if y.sum() else math.nan,
        "mean_detection_power_clean": float(d[y == 0].mean()) if (y == 0).any() else math.nan,
        "effective_clean_counterevidence": float(d[y == 0].sum()),
        "maximum_sampled_members": float(frame.mean_sampled_members.max()),
        "median_sampled_members": float(frame.mean_sampled_members.median()),
    }, float(log_lr[best])


def study_weights(frame: pd.DataFrame) -> np.ndarray:
    """Fixed publication weights, normalized within a focal group.

    Weights alter the evidence statistic, not the Bernoulli probability model.
    Unit or constant weights reproduce the unweighted statistic exactly.
    """
    w = frame.get("study_weight", pd.Series(1.0, index=frame.index)).to_numpy(float)
    if not np.all(np.isfinite(w)) or np.any(w <= 0):
        raise ValueError("study weights must be finite and strictly positive")
    return w / w.mean()


def bootstrap_p(
    frame: pd.DataFrame, rho_grid: np.ndarray, observed: float,
    replicates: int, seed: int,
) -> float:
    f = frame.background_apparent_failure_probability.to_numpy(float)
    d = frame.cross_study_detection_probability.to_numpy(float)
    probability = np.clip(
        f[:, None] + (1-f[:, None])*d[:, None]*rho_grid[None, :], 1e-12, 1-1e-12
    )
    null = np.clip(f[:, None], 1e-12, 1-1e-12)
    zero = np.log1p(-probability) - np.log1p(-null)
    one = np.log(probability) - np.log(null)
    weights = study_weights(frame)
    base = (zero * weights[:, None]).sum(axis=0)
    slope = (one-zero) * weights[:, None]
    rng = np.random.default_rng(seed)
    exceed = 0
    remaining = replicates
    while remaining:
        size = min(5000, remaining)
        simulated = rng.binomial(1, f, size=(size, len(frame)))
        best = np.max(base[None, :] + simulated @ slope, axis=1)
        exceed += int(np.sum(best >= observed - 1e-12))
        remaining -= size
    return (exceed + 1) / (replicates + 1)


def test_recurrence(
    observations: pd.DataFrame,
    universe: pd.DataFrame,
    background_folds: int,
    rho_grid_size: int,
    bootstrap_replicates: int,
    bootstrap_screen_q: float,
    filter_min: int,
    filter_max: int,
    fdr_q: float,
    background_seed: int,
    bootstrap_seed: int,
    output: str | Path,
) -> pd.DataFrame:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    fitted, strength = add_background(observations)
    fitted, background_metrics = adjust_background_for_sampling(
        fitted, background_folds, background_seed
    )
    fitted.to_csv(output / "publication_evidence.tsv", sep="\t", index=False)
    rho_grid = np.linspace(0, 1, rho_grid_size)
    rows, cache = [], {}
    for taxid, frame in fitted.groupby("taxid", sort=False):
        result, best = fit_group(frame, rho_grid)
        rows.append(result)
        cache[int(taxid)] = (best, frame)
    evidence = pd.DataFrame(rows)
    missing = universe[~universe.taxid.astype(int).isin(evidence.taxid.astype(int))]
    if not missing.empty:
        evidence = pd.concat([evidence, pd.DataFrame({
            "taxid": missing.taxid.astype(int), "name": missing.name, "rank": missing["rank"],
            "informative_publications": 0, "violating_publications": 0,
            "clean_publications": 0, "fitted_reproducible_nonmonophyly_tendency": 0.0,
            "likelihood_ratio_statistic": 0.0, "one_sided_asymptotic_p": 1.0,
        })], ignore_index=True)
    evidence["asymptotic_bh_q_full_taxon_universe"] = bh_adjust(
        evidence.one_sided_asymptotic_p.to_numpy(float)
    )
    evidence["parametric_bootstrap_p"] = np.nan
    weighted = "study_weight" in observations and observations.study_weight.nunique() > 1
    # A weighted likelihood statistic does not have the standard chi-square
    # approximation. Calibrate every observed group rather than use that screen.
    screen = (evidence.informative_publications.gt(0) if weighted else
              evidence.asymptotic_bh_q_full_taxon_universe.le(bootstrap_screen_q))
    if weighted:
        evidence["one_sided_asymptotic_p"] = np.nan
        evidence["asymptotic_bh_q_full_taxon_universe"] = np.nan
    for row in evidence[screen].itertuples(index=False):
        if int(row.taxid) not in cache:
            continue
        best, frame = cache[int(row.taxid)]
        evidence.loc[evidence.taxid.eq(row.taxid), "parametric_bootstrap_p"] = bootstrap_p(
            frame, rho_grid, best, bootstrap_replicates,
            stable_seed(bootstrap_seed, int(row.taxid))
        )
    evidence["test_p"] = evidence.parametric_bootstrap_p.fillna(1.0)
    profiles, evaluated = [], {}
    for cutoff in range(filter_min, filter_max + 1):
        mask = evidence.informative_publications.ge(cutoff).to_numpy()
        qvalues = bh_adjust(evidence.loc[mask, "test_p"].to_numpy(float))
        profiles.append({
            "minimum_informative_publications": cutoff,
            "tested_groups": int(mask.sum()),
            "discoveries": int((qvalues <= fdr_q).sum()),
        })
        evaluated[cutoff] = (mask, qvalues)
    profile = pd.DataFrame(profiles)
    maximum = int(profile.discoveries.max())
    cutoff = int(profile.loc[profile.discoveries.eq(maximum), "minimum_informative_publications"].min())
    mask, qvalues = evaluated[cutoff]
    evidence["passes_independent_publication_filter"] = mask
    evidence["filtered_bh_q"] = 1.0
    evidence.loc[mask, "filtered_bh_q"] = qvalues
    evidence["recommended_constraint_action"] = "retain_original_ncbi_clade_requirement"
    selected_mask = evidence.filtered_bh_q.le(fdr_q)
    evidence.loc[selected_mask, "recommended_constraint_action"] = "remove_exact_ncbi_clade_requirement"
    selected = evidence[selected_mask].sort_values(["filtered_bh_q", "taxid"])
    evidence.to_csv(output / "all_group_evidence.tsv", sep="\t", index=False)
    selected.to_csv(output / "supported_original_constraint_deletions.tsv", sep="\t", index=False)
    profile.to_csv(output / "independent_filter_profile.tsv", sep="\t", index=False)
    actions = {
        "remove_original_monophyly_taxids": selected.taxid.astype(int).tolist(),
        "add_synthetic_monophyly_groups": [],
        "default_policy": "deletion_only",
        "application_rule": "remove only the exact selected NCBI clade assertion",
    }
    (output / "constraint_changes.json").write_text(json.dumps(actions, indent=2) + "\n")
    summary = {
        "represented_ncbi_groups": len(evidence),
        "selected_minimum_informative_publications": cutoff,
        "groups_in_filtered_fdr_family": int(mask.sum()),
        "selected_constraint_exceptions": len(selected),
        "background_prior_strength": strength,
        "metadata_study_weights_used": bool(weighted),
        "weighting_scope": "recurrence statistic and matching null simulations only",
        **background_metrics,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return selected
