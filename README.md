# Phylo

*Package for creating and manipulating phylogenies*

| **Documentation**                               | **PackageEvaluator**            | **Build Status of master**                                                    |
|:-----------------------------------------------:|:------------------------:|:-------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![][pkg-0.6-img]][pkg-0.6-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url]     |
| [![][docs-latest-img]][docs-latest-url]         | [![][pkg-0.7-img]][pkg-0.7-url] | [![][codecov-img]][codecov-url] [![][coveralls-img]][coveralls-url] |

## Installation

The package is registered in `METADATA.jl` and so can be installed with `Pkg.add`.

```julia
julia> Pkg.add("Phylo")
```

## Project Status

The package is tested against the current Julia `0.6` release, `0.7` and `1.0` on Linux, OS X, and Windows. some dependencies are
currently broken on 0.7 and 1.0, but `Phylo` works without them for now.
Hopefully it will fully functional soon (currently missing an R interface)...

## Contributing and Questions

Contributions are very welcome, as are feature requests and suggestions.
Please open an [issue][issues-url] if you encounter any problems or would
just like to ask a question.

## Summary

**Phylo** is a [Julia](http://www.julialang.org) package that provides
 functionality for generating phylogenetic trees to feed into our
 [Diversity][diversity-url] package to calculate phylogenetic
 diversity. `Phylo` is currently in *alpha*, and is missing much
 functionality that people may desire, so please
 [raise an issue][issues-url] if/when you find problems or missing
 functionality - don't assume that I know! Currently the package can
 be used to make trees manually, to generate random trees using the
 framework from `Distributions`, and to read newick and nexus format trees from
 files. For instance, to construct a sampler for 5 tip non-ultrametric
 trees, and then generate one or two random tree of that type:

```julia
julia> using Phylo

julia> nu = Nonultrametric(5);

julia> tree = rand(nu)
Phylo.BinaryTree{DataFrames.DataFrame,Dict{String,Any}} phylogenetic tree with 9 nodes and 8 branches
Leaf names:
String["tip 1", "tip 2", "tip 3", "tip 4", "tip 5"]

julia> trees = rand(nu, ["Tree 1", "Tree 2"])
TreeSet with 2 trees
```

The code also provides iterators, and filtered iterators over the
branches, nodes, branchnames and nodenames of a tree:

```julia
julia> collect(nodeiter(tree))
9-element Array{Phylo.BinaryNode{Int64},1}:
 [branch 8]-->[leaf node]
 [branch 4]-->[leaf node]
 [branch 3]-->[leaf node]
 [branch 1]-->[leaf node]
 [branch 2]-->[leaf node]
 [branch 6]-->[internal node]-->[branches 1 and 2]
 [branch 5]-->[internal node]-->[branches 3 and 4]
 [branch 7]-->[internal node]-->[branches 5 and 6]
 [root node]-->[branches 7 and 8]

julia> collect(nodenamefilter(isroot, tree))
1-element Array{String,1}:
 "Node 4"
```

TreeSets are iterators themselves

```julia
julia> collect(trees)
2-element Array{Phylo.BinaryTree{DataFrames.DataFrame,Dict{String,Any}},1}:
 Phylo.BinaryTree{DataFrames.DataFrame,Dict{String,Any}} phylogenetic tree with 9 nodes (5 leaves) and 8 branches
 Phylo.BinaryTree{DataFrames.DataFrame,Dict{String,Any}} phylogenetic tree with 9 nodes (5 leaves) and 8 branches
```

The current main purpose of this package is to provide a framework for
phylogenetics to use in our [Diversity][diversity-url] package, and
they will both be adapted as appropriate until both are functioning as
required (though they are currently working together reasonably successfully).

It can also read newick trees either from
strings or files:

```julia
julia> using Phylo

julia> simpletree = parsenewick("((,Tip:1.0)Internal,)Root;")
NamedTree phylogenetic tree with 5 nodes and 4 branches
Leaf names:
String["Node 2", "Tip", "Node 1"]

julia> getbranches(simpletree)
Dict{Int64,Phylo.Branch{String}} with 4 entries:
  4 => [node "Root"]-->[NaN length branch]-->[node "Node 2"]
  2 => [node "Internal"]-->[1.0 length branch]-->[node "Tip"]
  3 => [node "Root"]-->[NaN length branch]-->[node "Internal"]
  1 => [node "Internal"]-->[NaN length branch]-->[node "Node 1"]

julia> tree = open(parsenewick, Pkg.dir("Phylo", "test", "H1N1.newick"))
Phylo.BinaryTree{Phylo.LeafInfo,Dict{String,Any}} phylogenetic tree with 1013 nodes and 1012 branches
Leaf names:
String["44", "429", "294", "295", "227", "14", "106", "104", "174", "331"  …  "384", "173", "300", "442", "215", "480", "477", "478", "30", "418"]
```
And it can read nexus trees from files too:

```julia
julia> ts = open(parsenexus, Pkg.dir("Phylo", "test", "H1N1.trees"))
Info: Created a tree called 'TREE1'
Info: Created a tree called 'TREE2'
TreeSet with 2 trees
Tree names:
String["TREE2", "TREE1"]


julia> ts["TREE1"]
Phylo.BinaryTree{DataFrames.DataFrame,Dict{String,Any}} phylogenetic tree with 1013 nodes and 1012 branches
Leaf names:
String["H1N1_A_MIYAGI_3_2000", "H1N1_A_PARMA_6_2008", "H1N1_A_AKITA_86_2002", "H1N1_A_DAKAR_14_1997", "H1N1_A_EGYPT_84_2001", "H1N1_A_NORWAY_69_2004", "H1N1_A_STPETERSBURG_11_2001", "H1N1_A_CAPETOWN_106_2007", "H1N1_A_SWITZERLAND_5773_2001", "H1N1_A_DENMARK_11_2005"  …  "H1N1_A_HONGKONG_1441_2006", "H1N1_A_BUCURESTI_288_2008", "H1N1_A_EGYPT_186_2006", "H1N1_A_HONGKONG_1134_1998", "H1N1_A_NEWCALEDONIA_20_1999", "H1N1_A_POLAND_W1_2001", "H1N1_A_MADRID_G793_1998", "H1N1_A_NORWAY_1730_2007", "H1N1_A_KALININGRAD_7_2008","H1N1_A_HONGKONG_2070_1999"]

julia> collect(treeinfoiter(ts))
2-element Array{Dict{String,Any},1}:
 Dict{String,Any}(Pair{String,Any}("lnP", -1.0))
 Dict{String,Any}(Pair{String,Any}("lnP", 1.0))

```

And while we wait for me (or kind [contributors][pr-url]!) to fill out
the other extensive functionality that many phylogenetics packages
have in other languages, the other important feature that it offers is
a fully(?)-functional interface to R, allowing any existing R library
functions to be carried out on julia trees, and trees to be read from
disk and written using R helper functions. Naturally the medium-term
plan is to fill in as many of these gaps as possible in Julia, so the R interface does not make RCall a dependency of the package (we use the
`Requires` package to avoid dependencies). Instead, if you want to use
the R interface you just need to use both packages:

```julia
julia> using Phylo

julia> using RCall
Creating Phylo RCall interface...

R> library(ape)
```

You can then translate back and forth using `NamedTree` contructors on
R `phylo` objects, and `RObject` constructors on julia `NamedTree`
types to keep them in Julia or `@rput` to move the object into R:

```julia
julia> rt = rcall(:rtree, 10)
RCall.RObject{RCall.VecSxp}

Phylogenetic tree with 10 tips and 9 internal nodes.

Tip labels:
	t10, t8, t1, t2, t6, t5, ...

Rooted; includes branch lengths.

julia> jt = NamedTree(rt)
Phylo.BinaryTree{DataFrames.DataFrame,Dict{String,Any}} phylogenetic tree with 19 nodes and 18 branches
Leaf names:
String["t2", "t1", "t5", "t9", "t8", "t3", "t4", "t10", "t7", "t6"]

julia> @rput rt;

julia> @rput jt; # Automatically translates jt back to R

R> jt

Phylogenetic tree with 10 tips and 9 internal nodes.

Tip labels:
	t10, t8, t1, t2, t6, t5, ...

Rooted; includes branch lengths.

R> all.equal(rt, jt) # check no damage in translations
[1] TRUE
```

For the time being the code will only work with rooted trees
with named tips and branch lengths. If there's [demand][issues-url]
for other types of trees, I'll look into it.

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://richardreeve.github.io/Phylo.jl/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://richardreeve.github.io/Phylo.jl/stable

[travis-img]: https://travis-ci.org/richardreeve/Phylo.jl.svg?branch=master
[travis-url]: https://travis-ci.org/richardreeve/Phylo.jl?branch=master

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/richardreeve/Phylo.jl?svg=true&branch=master
[appveyor-url]: https://ci.appveyor.com/project/richardreeve/phylo-jl/branch/master

[coveralls-img]: https://img.shields.io/coveralls/richardreeve/Phylo.jl.svg
[coveralls-url]: https://coveralls.io/r/richardreeve/Phylo.jl?branch=master

[codecov-img]: https://codecov.io/gh/richardreeve/Phylo.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/richardreeve/Phylo.jl

[pkg-0.6-img]: http://pkg.julialang.org/badges/Phylo_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=Phylo&ver=0.6

[pkg-0.7-img]: http://pkg.julialang.org/badges/Phylo_0.7.svg
[pkg-0.7-url]: http://pkg.julialang.org/?pkg=Phylo&ver=0.7

[issues-url]: https://github.com/richardreeve/Phylo.jl/issues
[pr-url]: https://github.com/richardreeve/Phylo.jl/pulls
[diversity-url]: https://github.com/richardreeve/Diversity.jl/
