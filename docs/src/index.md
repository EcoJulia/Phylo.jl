# Phylo

*Package for creating and manipulating phylogenies*

**Phylo** is a [Julia](http://www.julialang.org) package that provides
 functionality for generating phylogenetic trees to feed into our
 [Diversity](https://github.com/richardreeve/Diversity.jl) package to calculate phylogenetic
 diversity (currently on master, accessible via `Pkg.checkout()`,
 but not released). Both are currently under development, so please
 [raise an issue](https://github.com/richardreeve/Phylo.jl/issues) if you find any problems. Currently the
 package can be used to make trees manually, and to generate random
 trees using the framework from `Distributions`. For instance, to construct a sampler for 5 tip non-ultrametric
 trees, and then generate one or two random tree of that type:

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

julia> include(joinpath(Pkg.dir("Phylo"), "src", "rcall.jl"));

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

```@contents
```

```@autodocs
Modules = [Phylo]
Private = false
```

```@index
```
