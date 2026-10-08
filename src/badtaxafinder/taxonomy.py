from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import re
import tarfile
from urllib.request import urlretrieve


EXCLUDED_LOW_RANKS = {
    "species", "subspecies", "varietas", "forma", "strain", "isolate",
}
EXPLICIT_FAMILY_OR_BELOW = {
    "family", "subfamily", "supertribe", "tribe", "subtribe",
    "genus", "subgenus",
}
PLACEHOLDER_TERMS = (
    "unclassified", "environmental sample", "incertae sedis",
)
TIP_TAXID = re.compile(r"^(?:tx|taxid[:=_-]?)(\d+)(?:_|$)", re.I)


def _parts(line: str) -> list[str]:
    return [value.strip() for value in line.split("|")]


@dataclass
class Taxonomy:
    parent: dict[int, int]
    rank: dict[int, str]
    name: dict[int, str]
    aliases: dict[str, set[int]]
    merged: dict[int, int]
    deleted: set[int]

    @classmethod
    def from_taxdump(cls, directory: str | Path) -> "Taxonomy":
        root = Path(directory)
        parent: dict[int, int] = {}
        rank: dict[int, str] = {}
        with (root / "nodes.dmp").open(encoding="utf-8") as handle:
            for line in handle:
                fields = _parts(line)
                taxid = int(fields[0])
                parent[taxid] = int(fields[1])
                rank[taxid] = fields[2]
        names: dict[int, str] = {}
        aliases: dict[str, set[int]] = defaultdict(set)
        with (root / "names.dmp").open(encoding="utf-8") as handle:
            for line in handle:
                fields = _parts(line)
                taxid, value, name_class = int(fields[0]), fields[1], fields[3]
                if name_class == "scientific name":
                    names[taxid] = value
                if name_class in {
                    "scientific name", "synonym", "equivalent name", "includes",
                    "in-part", "misspelling", "genbank synonym",
                }:
                    for key in {value, value.replace(" ", "_")}:
                        aliases[key.casefold()].add(taxid)
        merged: dict[int, int] = {}
        if (root / "merged.dmp").exists():
            with (root / "merged.dmp").open(encoding="utf-8") as handle:
                for line in handle:
                    fields = _parts(line)
                    merged[int(fields[0])] = int(fields[1])
        deleted: set[int] = set()
        if (root / "delnodes.dmp").exists():
            with (root / "delnodes.dmp").open(encoding="utf-8") as handle:
                deleted = {int(_parts(line)[0]) for line in handle}
        return cls(parent, rank, names, dict(aliases), merged, deleted)

    def current_taxid(self, taxid: int) -> int | None:
        seen: set[int] = set()
        while taxid in self.merged and taxid not in seen:
            seen.add(taxid)
            taxid = self.merged[taxid]
        return None if taxid in self.deleted or taxid not in self.parent else taxid

    def resolve_tip(self, label: str) -> int | None:
        value = label.strip().strip("'\"")
        match = TIP_TAXID.match(value)
        if match:
            return self.current_taxid(int(match.group(1)))
        if value.isdigit():
            return self.current_taxid(int(value))
        matches = self.aliases.get(value.casefold(), set())
        current = {self.current_taxid(t) for t in matches}
        current.discard(None)
        return next(iter(current)) if len(current) == 1 else None

    def lineage(self, taxid: int) -> tuple[int, ...]:
        result: list[int] = []
        seen: set[int] = set()
        while taxid not in seen:
            result.append(taxid)
            seen.add(taxid)
            nxt = self.parent.get(taxid, taxid)
            if nxt == taxid:
                break
            taxid = nxt
        return tuple(result)

    def child_under(self, taxid: int, ancestor: int) -> int | None:
        lineage = self.lineage(taxid)
        try:
            position = lineage.index(ancestor)
        except ValueError:
            return None
        return lineage[position - 1] if position else taxid

    def broad_region(self, taxid: int) -> int:
        lineage = self.lineage(taxid)
        for wanted in ("order", "class", "phylum"):
            for node in lineage:
                if self.rank.get(node) == wanted:
                    return node
        return lineage[-2] if len(lineage) > 1 else lineage[-1]

    def testable_ancestors(self, species: int) -> tuple[int, ...]:
        lineage = self.lineage(species)
        result: list[int] = []
        below_family = False
        for node in reversed(lineage):
            node_rank = self.rank.get(node, "")
            if node_rank == "family":
                below_family = True
            if not below_family and node_rank not in EXPLICIT_FAMILY_OR_BELOW:
                continue
            if node_rank in EXCLUDED_LOW_RANKS:
                continue
            label = self.name.get(node, "")
            if not label or any(term in label.casefold() for term in PLACEHOLDER_TERMS):
                continue
            if node_rank in EXPLICIT_FAMILY_OR_BELOW or node_rank in {"no rank", ""}:
                result.append(node)
        return tuple(result)


def download_taxdump(output: str | Path) -> Path:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    archive = output / "taxdump.tar.gz"
    urlretrieve("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", archive)
    with tarfile.open(archive, "r:gz") as handle:
        wanted = {"nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"}
        for member in handle.getmembers():
            if Path(member.name).name in wanted:
                member.name = Path(member.name).name
                handle.extract(member, output, filter="data")
    return output

