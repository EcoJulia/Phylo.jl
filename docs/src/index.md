# Phylo

*Package for creating and manipulating phylogenies*

**Phylo** is a [Julia](http://www.julialang.org) package that provides
 functionality for generating phylogenetic trees to feed into our
 [Diversity][diversity-url] package to calculate phylogenetic
 diversity (currently on master, accessible via `Pkg.checkout()`,
 but not released). Both are currently under development, so please
 [raise an issue][issues-url] if you find any problems. Currently the
 package can be used to make trees manually, and to generate random
 trees using the framework from `Distributions`. For instance, to
 construct a sampler for 5 tip non-ultrametric trees, and then
 generate a random tree of that type:

```julia
julia> using Phylo

julia> nu = Nonultrametric(5);

julia> tree = rand(nu)
NamedTree phylogenetic tree with 9 nodes and 8 branches
Leaf names:
String["tip 1", "tip 2", "tip 3", "tip 4", "tip 5"]
```

The code also provides iterators, and filtered iterators over the
branches, nodes, branchnames and nodenames of a tree:

```julia
julia> collect(nodeiter(tree))
9-element Array{Phylo.BinaryNode{Int64},1}:
 [branch 4]-->[leaf node]
 [branch 5]-->[leaf node]
 [branch 2]-->[leaf node]
 [branch 1]-->[leaf node]
 [branch 8]-->[leaf node]
 [branch 3]-->[internal node]-->[branches 1 and 2]
 [branch 6]-->[internal node]-->[branches 3 and 4]
 [branch 7]-->[internal node]-->[branches 5 and 6]
 [root node]-->[branches 7 and 8]

julia> collect(nodenamefilter(isroot, tree))
1-element Array{String,1}:
 "Node 4"
```

The current main purpose of this package is to provide a framework for
phylogenetics to use in our [Diversity][diversity-url] package, and
they will both be adapted as appropriate until both are functioning as
required (though they are currently working together reasonably successfully).

However, it can also read newick trees is a very hacky way:

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

julia> open(parsenewick, "h1n1.trees")
NamedTree phylogenetic tree with 1013 nodes and 1012 branches
Leaf names:
String["407", "153", "1", "54", "101", "371", "41", "464", "65", "475"  …  "336", "145", "36", "95", "414", "138", "294", "353", "232", "306"]
```

And while we wait for me (or kind [contributors][pr-url]!) to fill out
the other extensive functionality that many phylogenetics packages
have in other languages, the other important feature that it offers is
a fully(?)-functional interface to R, allowing any existing R library
functions to be carried out on julia trees, and trees to be read from
disk and written using R helper functions. Naturally the medium-term
plan is to fill in as many of these gaps as possible in Julia, and as
a result this R interface is not built into the package as it will make
RCall (and R) a dependency, which I wanted to avoid. Instead, if you
want to use the R interface you need to do it manually, as below:

```julia
julia> using RCall

julia> include(joinpath(Pkg.dir("Phylo"), "src", "rcall.jl"));

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
NamedTree phylogenetic tree with 19 nodes and 18 branches
Leaf names:
String["t10", "t8", "t1", "t2", "t6", "t5", "t3", "t4", "t7", "t9"]

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

[issues-url]: https://github.com/richardreeve/Phylo.jl/issues
[diversity-url]: https://github.com/richardreeve/Diversity.jl/

```@contents
```

```@autodocs
Modules = [Phylo]
Private = false
```

```@index
```
