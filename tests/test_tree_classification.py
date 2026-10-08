from io import StringIO

from Bio import Phylo

from badtaxafinder.tree import IndexedTree, classify_group, descendant_counts


def classify(newick, members):
    indexed = IndexedTree(Phylo.read(StringIO(newick), "newick"))
    nodes = [indexed.tip_index[x] for x in members]
    counts = descendant_counts(indexed, indexed.tip_index.values())
    return classify_group(indexed, nodes, counts)


def test_displayed_group():
    assert classify("((a,b),(c,d));", ["a", "b"]) == "displayed"


def test_contradicted_group():
    assert classify("((a,c),(b,d));", ["a", "b"]) == "contradicted"


def test_polytomy_is_compatible_unresolved():
    assert classify("(a,b,c,d);", ["a", "b"]) == "compatible_unresolved"

