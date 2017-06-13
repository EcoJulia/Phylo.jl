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

julia> rand(nu)
NamedTree phylogenetic tree with 9 nodes and 8 branches
Leaf names:
String["tip 1", "tip 2", "tip 3", "tip 4", "tip 5"]
```

The main purpose of this package is to provide a framework for
phylogenetics to use in our [Diversity][diversity-url] package, and
they will [both](https://github.com/richardreeve/Diversity.jl/pull/18)
be adapted as appropriate until both are functioning as required.

However, the other important feature that it holds is to allow an
interface to R, allowing any existing R functionality to be carried
out on julia trees, and trees to be read from disk and written using R
helper functions. This is not built into the package as it will make
RCall (and R) a dependency, which I wanted to avoid. Instead, for now,
if you want to use the R interface you need to do it manually, as
below:

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
