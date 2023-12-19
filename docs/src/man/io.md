# Creating phylogenies

## Reading phylogenies from disk

Phylo can read newick trees either from strings,

```@example io
using Phylo
simpletree = parsenewick("((,Tip:1.0)Internal,)Root;")
```

which will result in the following tree:

```@example io
getbranches(simpletree)
```

or from files

```@example io
tree = open(parsenewick, Phylo.path("H1N1.newick"))
```

It can read nexus trees from files too:

```@example io
ts = open(parsenexus, Phylo.path("H1N1.trees"))
```

Reading multiple trees from a nexus file returns a `TreeSet` - index to get
the individual trees

```@example reading
gettreeinfo(ts)
```

```@example reading
ts["TREE1"]
```

## Writing phylogenies to disk

Phylo can write newick trees either to strings,

```@example io
out = Phylo.outputtree(simpletree, Newick())
```

or to files

```@example io
Phylo.write("test.newick", simpletree)
```

It can write nexus trees to files too:

```@example io
Phylo.write("h1.trees", ts)
```

It will use newick as the default format for OneTree trees (e.g. a RecursiveTree),
and nexus for ManyTrees trees (e.g. a TreeSet). However, you can tell it to use nexus
for a OneTree:

```@example io
Phylo.write("test.trees", simpletree, Nexus())
```

## Creating random phylogenies

The package can be used to generate random trees using the framework from
 `Distributions`. For instance, to construct a sampler for 5 tip non-ultrametric
 trees and generate a random tree of that type

```@example random_trees
using Phylo
nu = Nonultrametric(5);
tree = rand(nu)
```

Or two trees

```@example random_trees
trees = rand(nu, ["Tree 1", "Tree 2"])
```

## Importing phylogenies from R

Phylo allows transferring trees from R's `ape` package directly via RCall.
This allows any existing R library functions to be carried out on julia trees.
Naturally the medium-term plan is to make this package feature-complete
with existing functionality in R, and as a result the R interface is not built
into the package, avoiding having RCall (and R) a dependency. Instead, if you
want to use the R interface you need to do it manually, as below:

```julia
julia> using RCall

julia> include(Phylo.path("rcall.jl", dir = "src"));

R> library(ape)
```

You can then translate back and forth using `rcopy` on
R `phylo` objects, and `RObject` constructors on julia `NamedTree`
types to keep them in Julia or `@rput` to move the object into R:

```julia
julia> rt = rcall(:rtree, 10)
RCall.RObject{RCall.VecSxp}

Phylogenetic tree with 10 tips and 9 internal nodes.

Tip labels:
  t10, t8, t1, t2, t6, t5, ...

Rooted; includes branch lengths.

julia> jt = rcopy(NamedTree, rt)
NamedTree with 10 tips, 19 nodes and 18 branches.
Leaf names are t8, t3, t7, t9, t6, ... [4 omitted] ... and t1

julia> rjt = RObject(jt); # manually translate it back to R

R> if (all.equal($rjt, $rt)) "no damage in translation"
[1] "no damage in translation"

julia> @rput rt; # Or use macros to pass R object back to R

julia> @rput jt; # And automatically translate jt back to R

R> jt

Phylogenetic tree with 10 tips and 9 internal nodes.

Tip labels:
  t10, t8, t1, t2, t6, t5, ...

Rooted; includes branch lengths.

R> if (all.equal(rt, jt)) "no damage in translation"
[1] "no damage in translation"
```

```@docs
parsenewick
parsenexus
Nonultrametric
Ultrametric
```
