# Phylo

*Package for creating and manipulating phylogenies*

| **Documentation** | **PackageEvaluator** | **Build Status of master** |
|:-----------------:|:--------------------:|:--------------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![][pkg-0.6-img]][pkg-0.6-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] |
| [![][docs-latest-img]][docs-latest-url] | [![Works with 1.0!][pkg-1.0-img]][pkg-1.0-url] | [![][codecov-img]][codecov-url] [![][coveralls-img]][coveralls-url] |

## Installation

The package is registered in `METADATA` on Julia v0.6 and the `General`
registry on v0.7 and v1.0 and so can be installed with `add`. For example
on Julia v1.0:

```julia
(v1.0) pkg> add Phylo
 Resolving package versions...
  Updating `~/.julia/environments/v1.0/Project.toml`
  [aea672f4] + Phylo v0.3.2
  Updating `~/.julia/environments/v1.0/Manifest.toml`
  [7d9fca2a] + Arpack v0.2.2
  [9e28174c] + BinDeps v0.8.9
  [31c24e10] + Distributions v0.16.2
  [90014a1f] + PDMats v0.9.4
  [aea672f4] + Phylo v0.3.2
  [1fd47b50] + QuadGK v2.0.0
  [79098fc4] + Rmath v0.5.0
  [276daf66] + SpecialFunctions v0.7.0
  [4c63d2b9] + StatsFuns v0.7.0
  [0796e94c] + Tokenize v0.5.2
  [30578b45] + URIParser v0.4.0
  [4607b0f0] + SuiteSparse

(v1.0) pkg>
```

## Project Status

The package is tested against the current Julia v1.0 release, but also
the previous v0.6 and v0.7 versions on Linux, macOS, and Windows.

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
framework from `Distributions`, and to read newick and nexus format
trees from files. For instance, to construct a sampler for 5 tip
non-ultrametric trees, and then generate one or two random tree of
that type:

```julia
julia> using Phylo

julia> nu = Nonultrametric(5);

julia> tree = rand(nu)
Phylogenetic tree with 5 tips, 9 nodes and 8 branches.
Leaf names are tip 1, tip 2, tip 3, tip 4 and tip 5

julia> trees = rand(nu, ["Tree 1", "Tree 2"])
TreeSet with 2 trees, each with 5 tips.
Tree names are Tree 2 and Tree 1

Tree 2: Phylogenetic tree with 5 tips,9 nodes and 8 branches.

Tree 1: Phylogenetic tree with 5 tips,9 nodes and 8 branches.
```

The code also provides iterators, and filtered iterators over the
branches, nodes, branchnames and nodenames of a tree:

```julia
julia> collect(nodeiter(tree))
9-element Array{BinaryNode{Int64},1}:
 [branch 6]-->[leaf node]
 [branch 1]-->[leaf node]
 [branch 4]-->[leaf node]
 [branch 3]-->[leaf node]
 [branch 2]-->[leaf node]
 [branch 5]-->[internal node]-->[branch 1]
                             -->[branch 2]
 [branch 7]-->[internal node]-->[branch 3]
                             -->[branch 4]
 [branch 8]-->[internal node]-->[branch 5]
                             -->[branch 6]
 [root node]-->[branch 7]
            -->[branch 8]

julia> collect(nodenamefilter(isroot, tree))
1-element Array{String,1}:
 "Node 4"
```

TreeSets are iterators themselves

```julia
julia> collect(trees)
2-element Array{BinaryTree{DataFrames.DataFrame,Dict{String,Any}},1}:
 Phylogenetic tree with 5 tips, 9 nodes and 8 branches.
Leaf names are tip 1, tip 2, tip 3, tip 4 and tip 5

 Phylogenetic tree with 5 tips,9 nodes and 8 branches.
Leaf names are tip 1, tip 2, tip 3, tip 4 and tip 5
...
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
BinaryTree{DataFrames.DataFrame,Dict{String,Any}} with 3 tips, 5 nodes and 4 branches.
Leaf names are Node 1, Tip and Node 2


julia> getbranches(simpletree)
Dict{Int64,Branch{String}} with 4 entries:
  4 => [node "Root"]-->[NaN length branch]-->[node "Node 2"]…
  2 => [node "Internal"]-->[NaN length branch]-->[node "Node 1"]…
  3 => [node "Root"]-->[NaN length branch]-->[node "Internal"]…
  1 => [node "Internal"]-->[1.0 length branch]-->[node "Tip"]…

julia> tree = open(parsenewick, Phylo.path("H1N1.newick"))
BinaryTree{DataFrames.DataFrame,Dict{String,Any}} with 507 tips, 1013 nodes and 1012 branches.
Leaf names are 44, 429, 294, 295, 227, ... [501 omitted] ... and 418
```
And it can read nexus trees from files too:

```julia
julia> ts = open(parsenexus, Phylo.path("H1N1.trees"))
[ Info: Created a tree called 'TREE1'
[ Info: Created a tree called 'TREE2'
TreeSet with 2 trees, each with 507 tips.
Tree names are TREE2 and TREE1

TREE2: BinaryTree{DataFrames.DataFrame,Dict{String,Any}} with 507 tips, 1013 nodes and 1012 branches.
Leaf names are H1N1_A_MIYAGI_3_2000, H1N1_A_PARMA_6_2008, H1N1_A_AKITA_86_2002, H1N1_A_DAKAR_14_1997, H1N1_A_EGYPT_84_2001, ... [501 omitted] ... and H1N1_A_HONGKONG_2070_1999

TREE1: BinaryTree{DataFrames.DataFrame,Dict{String,Any}} with 507 tips, 1013 nodes and 1012 branches.
Leaf names are H1N1_A_MIYAGI_3_2000, H1N1_A_PARMA_6_2008, H1N1_A_AKITA_86_2002, H1N1_A_DAKAR_14_1997, H1N1_A_EGYPT_84_2001, ... [501 omitted] ... and H1N1_A_HONGKONG_2070_1999

julia> ts["TREE1"]
BinaryTree{DataFrames.DataFrame,Dict{String,Any}} with 507 tips, 1013 nodes and 1012 branches.
Leaf names are H1N1_A_MIYAGI_3_2000, H1N1_A_PARMA_6_2008, H1N1_A_AKITA_86_2002, H1N1_A_DAKAR_14_1997, H1N1_A_EGYPT_84_2001, ... [501 omitted] ... and H1N1_A_HONGKONG_2070_1999

julia> collect(treeinfoiter(ts))
2-element Array{Pair{String,Dict{String,Any}},1}:
 "TREE2" => Dict("lnP"=>-1.0)
 "TREE1" => Dict("lnP"=>1.0)
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
Phylo.BinaryTree{DataFrames.DataFrame,Dict{String,Any}} phylogenetic tree with 19 nodes and 18 branches
Leaf names:
String["t2", "t1", "t5", "t9", "t8", "t3", "t4", "t10", "t7", "t6"]

julia rjt = RObject(jt); # manually translate it back to R

R> all.equal($rjt, $rt) # check no damage in translations
[1] TRUE

julia> @rput rt; # Or use macros to pass R object back to R

julia> @rput jt; # And automatically translate jt back to R

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

[pkg-1.0-img]: http://pkg.julialang.org/badges/Phylo_1.0.svg
[pkg-1.0-url]: http://pkg.julialang.org/?pkg=Phylo&ver=1.0

[issues-url]: https://github.com/richardreeve/Phylo.jl/issues
[pr-url]: https://github.com/richardreeve/Phylo.jl/pulls
[diversity-url]: https://github.com/richardreeve/Diversity.jl/
