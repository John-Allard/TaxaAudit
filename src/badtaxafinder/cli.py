from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

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
    parser = argparse.ArgumentParser(prog="badtaxafinder")
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
