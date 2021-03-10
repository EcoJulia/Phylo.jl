# Phylo

*Package for creating and manipulating phylogenies*

**Phylo** is a [Julia](http://www.julialang.org) package that provides
 functionality for generating phylogenetic trees to feed into our
 [Diversity](https://github.com/richardreeve/Diversity.jl) package to calculate phylogenetic
 diversity. Both are currently under development, so please
 [raise an issue](https://github.com/richardreeve/Phylo.jl/issues) if you find any problems. Currently the
 package can be used to make trees manually, and to generate random
 trees using the framework from `Distributions`. For instance, to construct a sampler for 5 tip non-ultrametric
 trees, and then generate one or two random tree of that type:

```julia
julia> using Phylo

julia> nu = Nonultrametric(5);

julia> tree = rand(nu)
RootedTree with 5 tips, 9 nodes and 8 branches.
Leaf names are tip 3, tip 5, tip 2, tip 4 and tip 1

julia> trees = rand(nu, ["Tree 1", "Tree 2"])
TreeSet with 2 trees, each with 5 tips.
Tree names are Tree 2 and Tree 1

Tree 1: RootedTree with 5 tips, 9 nodes and 8 branches.
Leaf names are tip 5, tip 3, tip 4, tip 2 and tip 1

Tree 2: RootedTree with 5 tips, 9 nodes and 8 branches.
Leaf names are tip 3, tip 5, tip 1, tip 2 and tip 4
```

The code also provides iterators, and filtered iterators over the
branches, nodes, branchnames and nodenames of a tree:

```julia
julia> collect(nodeiter(tree))
9-element Vector{LinkNode{OneRoot, String, Dict{String, Any}, LinkBranch{OneRoot, String, Dict{String, Any}, Float64}}}:
 LinkNode tip 1, a tip of the tree with an incoming connection (branch 16).

 LinkNode tip 2, a tip of the tree with an incoming connection (branch 10).

 LinkNode tip 3, a tip of the tree with an incoming connection (branch 11).

 LinkNode tip 4, a tip of the tree with an incoming connection (branch 14).

 LinkNode tip 5, a tip of the tree with an incoming connection (branch 9).

 LinkNode Node 6, an internal node with 1 inbound and 2 outbound connections (branches 12 and 9, 10)

 LinkNode Node 7, an internal node with 1 inbound and 2 outbound connections (branches 13 and 11, 12)

 LinkNode Node 8, an internal node with 1 inbound and 2 outbound connections (branches 15 and 13, 14)

 LinkNode Node 9, a root node with 2 outbound connections (branches 15, 16)

julia> collect(nodenamefilter(isroot, tree))
1-element Vector{String}:
 "Node 9"
```

TreeSets are iterators themselves

```julia
julia> collect(trees)

2-element Vector{RootedTree}:
 RootedTree with 5 tips, 9 nodes and 8 branches.
Leaf names are tip 3, tip 5, tip 1, tip 2 and tip 4

 RootedTree with 5 tips, 9 nodes and 8 branches.
Leaf names are tip 3, tip 5, tip 1, tip 2 and tip 4
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
RootedTree with 3 tips, 5 nodes and 4 branches.
Leaf names are Tip, Node 1 and Node 4


julia> getbranches(simpletree)
skipmissing(Union{Missing, LinkBranch{OneRoot, String, Dict{String, Any}, Float64}}[LinkBranch 1, from node Internal to node Tip (length 1.0).
, LinkBranch 2, from node Internal to node Node 1.
, LinkBranch 3, from node Root to node Internal.
, LinkBranch 4, from node Root to node Node 4.
])

julia> tree = open(parsenewick, Phylo.path("H1N1.newick"))
RootedTree with 507 tips, 1013 nodes and 1012 branches.
Leaf names are 227, 294, 295, 110, 390, ... [501 omitted] ... and 418
```
And it can read nexus trees from files too:

```julia
jjulia> ts = open(parsenexus, Phylo.path("H1N1.trees"))
[ Info: Created a tree called "TREE1"
[ Info: Created a tree called "TREE2"
TreeSet with 2 trees, each with 507 tips.
Tree names are TREE2 and TREE1

TREE1: RootedTree with 507 tips, 1013 nodes and 1012 branches.
Leaf names are H1N1_A_BRAZIL_11_1978, H1N1_A_TAHITI_8_1998, H1N1_A_TAIWAN_1_1986, H1N1_A_BAYERN_7_1995, H1N1_A_ENGLAND_45_1998, ... [501 omitted] ... and H1N1_A_PUERTORICO_8_1934

TREE2: RootedTree with 507 tips, 1013 nodes and 1012 branches.
Leaf names are H1N1_A_BRAZIL_11_1978, H1N1_A_TAHITI_8_1998, H1N1_A_TAIWAN_1_1986, H1N1_A_BAYERN_7_1995, H1N1_A_ENGLAND_45_1998, ... [501 omitted] ... and H1N1_A_PUERTORICO_8_1934

julia> ts["TREE1"]
RootedTree with 507 tips, 1013 nodes and 1012 branches.
Leaf names are H1N1_A_BRAZIL_11_1978, H1N1_A_TAHITI_8_1998, H1N1_A_TAIWAN_1_1986, H1N1_A_BAYERN_7_1995, H1N1_A_ENGLAND_45_1998, ... [501 omitted] ... and H1N1_A_PUERTORICO_8_1934

julia> gettreeinfo(ts)
Dict{String, Dict{String, Any}} with 2 entries:
  "TREE2" => Dict("lnP"=>-1.0)
  "TREE1" => Dict("lnP"=>1.0)
```

And while we wait for me (or kind [contributors](https://github.com/richardreeve/Phylo.jl/pulls)!) to fill out
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

```@contents
```

```@autodocs
Modules = [Phylo]
Private = false
```

```@index
```
