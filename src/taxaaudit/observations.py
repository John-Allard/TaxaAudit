from __future__ import annotations

from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import statistics

import numpy as np
import pandas as pd

from .taxonomy import Taxonomy
from .tree import classify_group, descendant_counts, read_mapped_tree


TREE_SUFFIXES = {".nwk", ".newick", ".tree", ".tre"}


def load_manifest(
    trees: str | None = None, manifest: str | None = None
) -> list[tuple[str, str]]:
    if bool(trees) == bool(manifest):
        raise ValueError("supply exactly one of --trees and --manifest")
    if trees:
        root = Path(trees).resolve()
        files = sorted(p for p in root.rglob("*") if p.suffix.lower() in TREE_SUFFIXES)
        return [(str(path), path.stem) for path in files]
    manifest_path = Path(str(manifest)).resolve()
    result: list[tuple[str, str]] = []
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"path", "publication"} <= set(reader.fieldnames):
            raise ValueError("manifest must contain path and publication columns")
        for row in reader:
            path = Path(row["path"])
            if not path.is_absolute():
                path = manifest_path.parent / path
            result.append((str(path.resolve()), row["publication"]))
    return result


def scan_trees(
    entries: list[tuple[str, str]], taxonomy: Taxonomy, output: str | Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    tree_rows: list[dict[str, object]] = []
    map_rows: list[dict[str, object]] = []
    universe: dict[int, tuple[str, str]] = {}
    for tree_id, (path, publication) in enumerate(entries):
        mapped = read_mapped_tree(path, publication, taxonomy)
        for label in mapped.unmapped_labels:
            map_rows.append({
                "tree_id": tree_id, "path": path, "publication": publication,
                "tip_label": label, "status": "unmapped_or_ambiguous",
            })
        if len(mapped.species_nodes) < 3:
            continue
        members_by_group: dict[int, list[int]] = defaultdict(list)
        for species in mapped.species_nodes:
            for group in taxonomy.testable_ancestors(species):
                members_by_group[group].append(species)
                universe[group] = (taxonomy.name.get(group, ""), taxonomy.rank.get(group, ""))
        counts = descendant_counts(mapped.indexed, mapped.species_nodes.values())
        for group, members in members_by_group.items():
            if len(members) < 2 or len(members) == len(mapped.species_nodes):
                continue
            nodes = [mapped.species_nodes[x] for x in members]
            status = classify_group(mapped.indexed, nodes, counts)
            parent = taxonomy.parent.get(group, group)
            member_set = set(members)
            local_nonmembers = [
                species for species in mapped.species_nodes
                if species not in member_set and parent in taxonomy.lineage(species)
            ]
            sibling_lineages = {
                taxonomy.child_under(species, parent) for species in local_nonmembers
            }
            sibling_lineages.discard(None)
            tree_rows.append({
                "taxid": group,
                "name": taxonomy.name.get(group, ""),
                "rank": taxonomy.rank.get(group, ""),
                "tree_id": tree_id,
                "publication": publication,
                "path": path,
                "mapped_species": len(mapped.species_nodes),
                "sampled_members": len(members),
                "sampled_nonmembers": len(mapped.species_nodes) - len(members),
                "sampled_parent_nonmembers": len(local_nonmembers),
                "represented_sibling_lineages": len(sibling_lineages),
                "region_taxid": taxonomy.broad_region(group),
                "clade_status": status,
            })
    tree_frame = pd.DataFrame(tree_rows)
    publication_rows: list[dict[str, object]] = []
    if not tree_frame.empty:
        for (_taxid, _publication), frame in tree_frame.groupby(
            ["taxid", "publication"], sort=False
        ):
            informative = frame[frame.clade_status.ne("compatible_unresolved")]
            if informative.empty:
                continue
            counts = informative.clade_status.value_counts()
            modes = counts[counts.eq(counts.max())].index.tolist()
            if len(modes) != 1:
                continue
            first = frame.iloc[0]
            publication_rows.append({
                "taxid": int(first.taxid), "name": first["name"], "rank": first["rank"],
                "publication": first.publication,
                "trees_in_publication": len(frame),
                "informative_trees_in_publication": len(informative),
                "publication_clade_status": modes[0],
                "observed_failure": int(modes[0] == "contradicted"),
                "mean_sampled_members": float(frame.sampled_members.mean()),
                "mean_sampled_nonmembers": float(frame.sampled_nonmembers.mean()),
                "mean_sampled_parent_nonmembers": float(
                    frame.sampled_parent_nonmembers.mean()
                ),
                "mean_represented_sibling_lineages": float(
                    frame.represented_sibling_lineages.mean()
                ),
                "region_taxid": int(first.region_taxid),
            })
    publication_frame = pd.DataFrame(publication_rows)
    universe_frame = pd.DataFrame(
        [
            {"taxid": taxid, "name": values[0], "rank": values[1]}
            for taxid, values in sorted(universe.items())
        ]
    )
    tree_frame.to_csv(output / "tree_group_observations.tsv", sep="\t", index=False)
    publication_frame.to_csv(
        output / "publication_group_observations.tsv", sep="\t", index=False
    )
    universe_frame.to_csv(output / "taxon_universe.tsv", sep="\t", index=False)
    pd.DataFrame(map_rows).to_csv(output / "unmapped_tips.tsv", sep="\t", index=False)
    (output / "tree_manifest.tsv").write_text(
        "path\tpublication\n" + "".join(f"{p}\t{s}\n" for p, s in entries),
        encoding="utf-8",
    )
    summary = {
        "tree_files": len(entries),
        "publications": len({publication for _path, publication in entries}),
        "tree_group_observations": len(tree_frame),
        "publication_group_observations": len(publication_frame),
        "represented_groups": len(universe_frame),
        "unmapped_or_ambiguous_tips": len(map_rows),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return tree_frame, publication_frame, universe_frame


def _stable_seed(seed: int, *parts: object) -> int:
    import hashlib
    digest = hashlib.sha256("\x1f".join(map(str, parts)).encode()).digest()
    return (seed + int.from_bytes(digest[:8], "big")) % (2**32)


def _select_count(rng: np.random.Generator, available: int, minimum: int) -> int:
    if available <= minimum:
        return available
    if rng.random() < 0.12:
        return available
    retention = math.exp(rng.uniform(math.log(max(minimum / available, 0.01)), 0.0))
    return min(available, max(minimum, int(rng.binomial(available, retention))))


def _sample(rng: np.random.Generator, values: list[int], count: int) -> list[int]:
    if count >= len(values):
        return list(values)
    return sorted(rng.choice(values, size=count, replace=False).tolist())


def _select_templates(frame: pd.DataFrame, maximum: int, seed: int) -> pd.DataFrame:
    contradicted = frame[frame.clade_status.eq("contradicted")].copy()
    if maximum <= 0 or len(contradicted) <= maximum:
        return contradicted
    contradicted["size_bin"] = np.floor(
        np.log2(contradicted.mapped_species.clip(lower=3))
    ).astype(int)
    contradicted["rank_bin"] = contradicted["rank"].fillna("").replace("", "no_rank")
    strata = list(contradicted.groupby(["rank_bin", "size_bin"], sort=False))
    allocation = max(1, maximum // max(1, len(strata)))
    pieces = [
        group.sample(
            n=min(len(group), allocation),
            random_state=_stable_seed(seed, "stratum", index),
        )
        for index, (_key, group) in enumerate(strata)
    ]
    selected = pd.concat(pieces).drop_duplicates()
    if len(selected) < maximum:
        remainder = contradicted.drop(index=selected.index)
        if not remainder.empty:
            selected = pd.concat([
                selected,
                remainder.sample(
                    n=min(maximum - len(selected), len(remainder)),
                    random_state=_stable_seed(seed, "remainder"),
                ),
            ])
    if len(selected) > maximum:
        selected = selected.sample(n=maximum, random_state=_stable_seed(seed, "trim"))
    return selected.drop(columns=["size_bin", "rank_bin"], errors="ignore")


def generate_detection_training(
    tree_observations: pd.DataFrame,
    taxonomy: Taxonomy,
    max_templates: int,
    draws_per_template: int,
    seed: int,
    output: str | Path,
) -> pd.DataFrame:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    templates = _select_templates(tree_observations, max_templates, seed)
    keys = {(int(r.tree_id), int(r.taxid)) for r in templates.itertuples()}
    reductions: list[dict[str, object]] = []
    for tree_id, frame in tree_observations.groupby("tree_id", sort=False):
        first = frame.iloc[0]
        mapped = read_mapped_tree(str(first.path), str(first.publication), taxonomy)
        species = set(mapped.species_nodes)
        lineages = {value: set(taxonomy.lineage(value)) for value in species}
        for row in frame.itertuples(index=False):
            if (int(tree_id), int(row.taxid)) not in keys:
                continue
            target = int(row.taxid)
            parent = taxonomy.parent.get(target, target)
            members = sorted(x for x in species if target in lineages[x])
            local = sorted(x for x in species - set(members) if parent in lineages[x])
            distant = sorted(species - set(members) - set(local))
            for draw in range(draws_per_template + 1):
                rng = np.random.default_rng(_stable_seed(seed, tree_id, target, draw))
                if draw == 0:
                    chosen_members, chosen_local, chosen_distant = members, local, distant
                else:
                    nm = _select_count(rng, len(members), 2)
                    nl = _select_count(rng, len(local), 0)
                    nd = _select_count(rng, len(distant), 0)
                    if nl + nd == 0:
                        if local:
                            nl = 1
                        elif distant:
                            nd = 1
                    chosen_members = _sample(rng, members, nm)
                    chosen_local = _sample(rng, local, nl)
                    chosen_distant = _sample(rng, distant, nd)
                retained = chosen_members + chosen_local + chosen_distant
                counts = descendant_counts(
                    mapped.indexed, [mapped.species_nodes[x] for x in retained]
                )
                status = classify_group(
                    mapped.indexed,
                    [mapped.species_nodes[x] for x in chosen_members],
                    counts,
                )
                siblings = {
                    taxonomy.child_under(x, parent) for x in chosen_local
                }
                siblings.discard(None)
                reductions.append({
                    "template_id": f"{tree_id}:{target}",
                    "tree_id": int(tree_id), "taxid": target,
                    "publication": row.publication,
                    "region_taxid": int(row.region_taxid), "draw": draw,
                    "is_full_template": int(draw == 0),
                    "sampled_members": len(chosen_members),
                    "sampled_parent_nonmembers": len(chosen_local),
                    "represented_sibling_lineages": len(siblings),
                    "detected_contradiction": int(status == "contradicted"),
                })
    result = pd.DataFrame(reductions)
    result.to_csv(output / "complete_group_downsampling.tsv", sep="\t", index=False)
    templates.to_csv(output / "template_manifest.tsv", sep="\t", index=False)
    return result

