# Traversing and iterating over trees

The code also provides iterators, and filtered iterators over the
branches, nodes, branchnames and nodenames of a tree (using the random tree from
[Creating and writing phylogenies](io.md))

```@example random_trees
using Phylo
nu = Nonultrametric(5);
tree = rand(nu)

collect(nodeiter(tree))
```

```@example random_trees
collect(nodenamefilter(isroot, tree))
```

TreeSets are iterators themselves

```@example random_trees
trees = rand(nu, ["Tree 1", "Tree 2"])
collect(trees)
```

```@docs
nodeiter
nodefilter
nodenameiter
nodenamefilter
branchiter
branchfilter
branchnameiter
branchnamefilter
traversal
```
