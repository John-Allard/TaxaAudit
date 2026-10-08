from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Sequence

from Bio import Phylo


class IndexedTree:
    def __init__(self, tree) -> None:
        nodes = list(tree.find_clades(order="preorder"))
        self.node_index = {node: i for i, node in enumerate(nodes)}
        parent0 = [0] * len(nodes)
        self.depth = [0] * len(nodes)
        for node in nodes:
            node_i = self.node_index[node]
            for child in node.clades:
                child_i = self.node_index[child]
                parent0[child_i] = node_i
                self.depth[child_i] = self.depth[node_i] + 1
        self.up = [parent0]
        for _ in range(1, max(1, len(nodes).bit_length())):
            prev = self.up[-1]
            self.up.append([prev[prev[i]] for i in range(len(nodes))])
        self.tip_index = {
            node.name: self.node_index[node]
            for node in nodes if node.is_terminal() and node.name is not None
        }

    def lca(self, a: int, b: int) -> int:
        if self.depth[a] < self.depth[b]:
            a, b = b, a
        delta = self.depth[a] - self.depth[b]
        bit = 0
        while delta:
            if delta & 1:
                a = self.up[bit][a]
            delta >>= 1
            bit += 1
        if a == b:
            return a
        for level in range(len(self.up) - 1, -1, -1):
            if self.up[level][a] != self.up[level][b]:
                a, b = self.up[level][a], self.up[level][b]
        return self.up[0][a]


def lca_many(indexed: IndexedTree, nodes: Sequence[int]) -> int:
    result = nodes[0]
    for node in nodes[1:]:
        result = indexed.lca(result, node)
    return result


def child_below(indexed: IndexedTree, ancestor: int, node: int) -> int:
    target_depth = indexed.depth[ancestor] + 1
    delta = indexed.depth[node] - target_depth
    bit = 0
    while delta:
        if delta & 1:
            node = indexed.up[bit][node]
        delta >>= 1
        bit += 1
    return node


def descendant_counts(indexed: IndexedTree, nodes: Iterable[int]) -> Counter[int]:
    counts: Counter[int] = Counter()
    for node in nodes:
        current = node
        while True:
            counts[current] += 1
            parent = indexed.up[0][current]
            if parent == current:
                break
            current = parent
    return counts


def classify_group(
    indexed: IndexedTree,
    member_nodes: Sequence[int],
    mapped_descendants: Counter[int],
) -> str:
    mrca = lca_many(indexed, member_nodes)
    if mapped_descendants[mrca] == len(member_nodes):
        return "displayed"
    members_by_child = Counter(child_below(indexed, mrca, n) for n in member_nodes)
    if any(count < mapped_descendants[child] for child, count in members_by_child.items()):
        return "contradicted"
    return "compatible_unresolved"


@dataclass
class MappedTree:
    path: str
    publication: str
    indexed: IndexedTree
    species_nodes: dict[int, int]
    unmapped_labels: tuple[str, ...]


def read_mapped_tree(path: str, publication: str, taxonomy) -> MappedTree:
    tree = Phylo.read(path, "newick")
    indexed = IndexedTree(tree)
    labels_by_species: dict[int, list[str]] = {}
    unmapped: list[str] = []
    for label in indexed.tip_index:
        taxid = taxonomy.resolve_tip(label)
        if taxid is None:
            unmapped.append(label)
        else:
            labels_by_species.setdefault(taxid, []).append(label)
    species_nodes = {
        taxid: indexed.tip_index[labels[0]]
        for taxid, labels in labels_by_species.items() if len(labels) == 1
    }
    return MappedTree(path, publication, indexed, species_nodes, tuple(unmapped))

