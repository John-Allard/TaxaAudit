from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import minimize_scalar

from taxaaudit.cli import main
from taxaaudit.upper_limits import (
    RecurrenceCalibrator, binary_records, calculate_upper_limits,
    fit_recurrence, log_likelihood,
)


def test_fit_matches_numerical_optimization():
    rng = np.random.default_rng(42)
    f = rng.uniform(.05, .5, 12)
    d = rng.uniform(.2, .95, 12)
    y = rng.binomial(1, .4, size=(20, 12))
    _, fitted = fit_recurrence(y, f, (1-f)*d)
    for record, best in zip(y, fitted):
        optimum = minimize_scalar(lambda r: -log_likelihood(record, f, (1-f)*d, r),
                                  bounds=(0, 1), method="bounded")
        assert best >= -optimum.fun - 1e-7


def test_exact_pvalues_and_coverage():
    f = np.array([.05, .2, .3])
    d = np.array([.9, .5, .8])
    records = binary_records(3)
    for rho in [0, .05, .3, .75, 1]:
        _, best = fit_recurrence(records, f, (1-f)*d)
        weights = np.exp(log_likelihood(records, f, (1-f)*d, rho))
        statistics = 2*(best - log_likelihood(records, f, (1-f)*d, rho))
        pvalues = []
        for row, stat in zip(records, statistics):
            actual = RecurrenceCalibrator(row, f, d, 2).pvalue(rho)[0]
            expected = weights[statistics >= stat-1e-9].sum()
            assert actual == pytest.approx(expected, abs=1e-10)
            pvalues.append(actual)
        assert weights[np.array(pvalues) > .05].sum() >= .95 - 1e-8


def test_mle_is_accepted_and_zero_detection_unrestricted():
    for y in [np.zeros(6), np.ones(6), np.array([0, 0, 1, 0, 1, 0])]:
        fit = RecurrenceCalibrator(y, np.full(6, .2), np.full(6, .7), 9)
        assert fit.pvalue(float(fit.fitted[0]))[0] == 1
        assert fit.upper() >= fit.fitted[0]
    assert RecurrenceCalibrator(np.zeros(5), np.full(5, .2), np.zeros(5), 0).upper() == 1


def test_mc_order_invariance():
    values = (np.zeros(13), np.full(13, .2), np.full(13, .7))
    a = RecurrenceCalibrator(*values, 9, replicates=100)
    first = a.pvalue(.2)
    a.pvalue(.8)
    b = RecurrenceCalibrator(*values, 9, replicates=100)
    b.pvalue(.8)
    assert b.pvalue(.2) == first


@pytest.mark.parametrize("y,f,d", [([2], [.2], [.8]), ([0], [np.nan], [.8]),
                                  ([0], [.2], [2]), ([], [], []),
                                  ([0, 1], [.2], [.8])])
def test_invalid_inputs_rejected(y, f, d):
    with pytest.raises(ValueError):
        RecurrenceCalibrator(y, f, d, 1)


def inputs(tmp_path):
    frame = pd.DataFrame(dict(taxid=[1]*3, name=["A"]*3, rank=["genus"]*3,
        publication=["p1", "p2", "p3"], observed_failure=[0, 1, 0],
        background_apparent_failure_probability=[.1]*3,
        cross_study_detection_probability=[.8]*3, mean_sampled_members=[3]*3))
    universe = pd.DataFrame(dict(taxid=[1, 2], name=["A", "B"], rank=["genus"]*2))
    flagged = pd.DataFrame(dict(taxid=[1]))
    paths = [tmp_path / name for name in ["observations.tsv", "universe.tsv", "flagged.tsv"]]
    for table, path in zip([frame, universe, flagged], paths):
        table.to_csv(path, sep="\t", index=False)
    return paths


def test_resource_resume_and_missing_evidence(tmp_path):
    paths = inputs(tmp_path)
    result = calculate_upper_limits(*paths, tmp_path / "output")
    missing = result.set_index("taxid").loc[2]
    assert pd.isna(missing.upper_recurrence_limit_95)
    assert missing.informative_publications == 0
    assert missing.evidence_status == "no_informative_publications"
    assert result.currently_flagged_recurrent_nonmonophyly.sum() == 1
    checkpoint = tmp_path / "output" / "upper_limits_checkpoint.jsonl"
    before = checkpoint.read_bytes()
    pd.testing.assert_frame_equal(result, calculate_upper_limits(*paths, tmp_path / "output"))
    assert checkpoint.read_bytes() == before
    with pytest.raises(ValueError, match="settings changed"):
        calculate_upper_limits(*paths, tmp_path / "output", seed=3)


def test_duplicate_and_nonuniform_weights_rejected(tmp_path):
    paths = inputs(tmp_path)
    frame = pd.read_csv(paths[0], sep="\t")
    frame["study_weight"] = [2, 1, 1]
    frame.to_csv(paths[0], sep="\t", index=False)
    with pytest.raises(ValueError, match="nonuniform"):
        calculate_upper_limits(*paths, tmp_path / "out")
    frame = frame.drop(columns="study_weight")
    frame.loc[1, "publication"] = "p1"
    frame.to_csv(paths[0], sep="\t", index=False)
    with pytest.raises(ValueError, match="at most one"):
        calculate_upper_limits(*paths, tmp_path / "out")


def test_cli_explicit_and_run_directory(tmp_path):
    paths = inputs(tmp_path)
    main(["upper-limits", "--observations", str(paths[0]), "--universe", str(paths[1]),
          "--catalog", str(paths[2]), "--out", str(tmp_path / "explicit")])
    catalog = tmp_path / "run" / "catalog"
    catalog.mkdir(parents=True)
    names = ["publication_evidence.tsv", "all_group_evidence.tsv",
             "supported_original_constraint_deletions.tsv"]
    for path, name in zip(paths, names):
        (catalog / name).write_bytes(path.read_bytes())
    main(["upper-limits", "--run", str(tmp_path / "run"), "--out", str(tmp_path / "via-run")])
    assert ((tmp_path / "explicit" / "taxaaudit_recurrence_upper_limits.tsv").read_bytes() ==
            (tmp_path / "via-run" / "taxaaudit_recurrence_upper_limits.tsv").read_bytes())


def test_released_resource_consistent():
    catalog = Path(__file__).parents[1] / "catalog"
    frame = pd.read_csv(catalog / "recurrence_upper_limits.tsv", sep="\t")
    flagged = pd.read_csv(catalog / "supported_original_constraint_deletions.tsv", sep="\t")
    assert len(frame) == 26291
    assert frame.taxid.is_unique
    assert frame.upper_recurrence_limit_95.notna().sum() == 18109
    assert frame.informative_publications.sum() == 73181
    assert set(frame.loc[frame.currently_flagged_recurrent_nonmonophyly, "taxid"]) == set(flagged.taxid)
    missing = frame.evidence_status.eq("no_informative_publications")
    assert missing.sum() == 8182
    assert frame.loc[missing, "upper_recurrence_limit_95"].isna().all()
    assert (frame.loc[~missing, "upper_recurrence_limit_95"] >=
            frame.loc[~missing, "rho_mle"] - 1e-8).all()


def test_scan_saves_probabilities_for_upper_limits(tmp_path, monkeypatch):
    from taxaaudit import model
    paths = inputs(tmp_path)
    frame = pd.read_csv(paths[0], sep="\t")
    universe = pd.read_csv(paths[1], sep="\t")
    monkeypatch.setattr(model, "add_background", lambda data: (data.copy(), 1.0))
    monkeypatch.setattr(model, "adjust_background_for_sampling",
                        lambda data, folds, seed: (data.copy(), {}))
    out = tmp_path / "run" / "catalog"
    model.test_recurrence(frame, universe, 2, 11, 10, .30, 2, 2, .05, 1, 1, out)
    saved = pd.read_csv(out / "publication_evidence.tsv", sep="\t")
    pd.testing.assert_frame_equal(saved, frame)
    main(["upper-limits", "--run", str(out.parent), "--out", str(tmp_path / "limits")])
    result = pd.read_csv(tmp_path / "limits" / "taxaaudit_recurrence_upper_limits.tsv", sep="\t")
    assert result.upper_recurrence_limit_95.notna().sum() == 1
