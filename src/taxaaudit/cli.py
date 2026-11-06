from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd
import numpy as np

from .model import fit_detection_probabilities, test_recurrence
from .observations import (
    generate_detection_training,
    load_manifest,
    scan_trees,
)
from .taxonomy import Taxonomy, download_taxdump


def add_run_arguments(parser: argparse.ArgumentParser) -> None:
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--trees", help="directory recursively containing Newick trees")
    source.add_argument("--manifest", help="TSV with path and publication columns")
    parser.add_argument("--taxonomy", required=True, help="unpacked NCBI taxdump directory")
    parser.add_argument("--out", required=True)
    parser.add_argument("--ranks", help="comma-separated explicit NCBI ranks; default is family and below")
    parser.add_argument("--study-weights", help="optional TSV: publication, weight. Unlisted publications receive 1; weights must be positive and fixed independently of focal outcomes")
    parser.add_argument("--max-templates", type=int, default=2000)
    parser.add_argument("--draws-per-template", type=int, default=40)
    parser.add_argument("--detection-folds", type=int, default=5)
    parser.add_argument("--background-folds", type=int, default=5)
    parser.add_argument("--rho-grid-size", type=int, default=201)
    parser.add_argument("--bootstrap-replicates", type=int, default=99999)
    parser.add_argument("--bootstrap-screen-q", type=float, default=0.30)
    parser.add_argument("--filter-min", type=int, default=2)
    parser.add_argument("--filter-max", type=int, default=20)
    parser.add_argument("--fdr-q", type=float, default=0.05)
    parser.add_argument(
        "--detection-seed", type=int, default=20260821,
        help="seed for template selection, downsampling, and detection fitting",
    )
    parser.add_argument(
        "--background-seed", type=int, default=20260728,
        help="seed for group-held-out background-model fitting",
    )
    parser.add_argument(
        "--bootstrap-seed", type=int, default=20260728,
        help="base seed for group-specific parametric bootstraps",
    )


def run(args: argparse.Namespace) -> None:
    output = Path(args.out).resolve()
    output.mkdir(parents=True, exist_ok=True)
    taxonomy = Taxonomy.from_taxdump(args.taxonomy)
    if args.ranks:
        ranks = frozenset(s.strip() for s in args.ranks.split(',') if s.strip())
        unknown = ranks - set(taxonomy.rank.values())
        if not ranks or unknown:
            raise ValueError(f"invalid or unknown NCBI ranks: {sorted(unknown)}")
        taxonomy.target_ranks = ranks
    entries = load_manifest(args.trees, args.manifest)
    if not entries:
        raise SystemExit("no Newick trees found")
    observations_dir = output / "observations"
    tree_frame, publication_frame, universe = scan_trees(
        entries, taxonomy, observations_dir
    )
    if publication_frame.empty:
        raise SystemExit("no informative group-publication observations were found")
    detection_dir = output / "detection_model"
    reductions = generate_detection_training(
        tree_frame, taxonomy, args.max_templates, args.draws_per_template,
        args.detection_seed, detection_dir,
    )
    if reductions.empty:
        raise SystemExit("no contradicted tree-group observations were available for detection training")
    predictions = fit_detection_probabilities(
        reductions, tree_frame, publication_frame, args.detection_folds,
        args.detection_seed, detection_dir,
    )
    if args.study_weights:
        weights = pd.read_csv(args.study_weights, sep="\t")
        if not {"publication", "weight"} <= set(weights):
            raise ValueError("study weights require publication and weight columns")
        if weights.publication.duplicated().any():
            raise ValueError("duplicate publication in weights file")
        weights["weight"] = pd.to_numeric(weights.weight, errors="raise")
        if not np.isfinite(weights.weight).all() or (weights.weight <= 0).any():
            raise ValueError("weights must be finite and positive")
        known = {publication for _, publication in entries}
        if not set(weights.publication) <= known:
            raise ValueError("weights file contains publications absent from the manifest")
        lookup = weights.set_index("publication").weight
        predictions["study_weight"] = predictions.publication.map(lookup).fillna(1.0)
        predictions[["publication", "study_weight"]].drop_duplicates().to_csv(
            output / "effective_study_weights.tsv", sep="\t", index=False)
    selected = test_recurrence(
        predictions, universe, args.background_folds, args.rho_grid_size,
        args.bootstrap_replicates, args.bootstrap_screen_q,
        args.filter_min, args.filter_max, args.fdr_q,
        args.background_seed, args.bootstrap_seed,
        output / "catalog",
    )
    settings = {
        key: value for key, value in vars(args).items()
        if key not in {"func"}
    }
    settings.update({
        "tree_files": len(entries),
        "publications": len({publication for _path, publication in entries}),
        "selected_groups": len(selected),
    })
    (output / "run_settings.json").write_text(
        json.dumps(settings, indent=2, default=str) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "output": str(output),
        "tree_files": len(entries),
        "selected_groups": len(selected),
    }, indent=2))


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="taxaaudit")
    subparsers = parser.add_subparsers(dest="command", required=True)
    download = subparsers.add_parser("download-taxonomy")
    download.add_argument("--out", required=True)
    download.set_defaults(func=lambda args: print(download_taxdump(args.out)))
    run_parser = subparsers.add_parser("run")
    add_run_arguments(run_parser)
    run_parser.set_defaults(func=run)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
