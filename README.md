# Phylo

*Package for creating and manipulating phylogenies*

| **PackageEvaluator**            | **Build Status of master**                                                    |
|:------------------------:|:-------------------------------------------------------------------:|
| [![][pkg-0.5-img]][pkg-0.5-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url]     |
| [![][pkg-0.6-img]][pkg-0.6-url] | [![][codecov-img]][codecov-url] [![][coveralls-img]][coveralls-url] |

**Phylo** is a [Julia](http://www.julialang.org) package that provides
 functionality for generating phylogenetic trees to feed into our
 [Diversity][diversity-url] package to calculate phylogenetic
 diversity (currently in the [phylogenetics][phylogenetics-url]
 branch). Both are currently under development, so please
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
String["tip 4","tip 1","tip 2","tip 3","tip 5"]

```

The main purpose of this package is to provide a framework for
phylogenetics to use in our [Diversity][diversity-url] package, and
they will both be adapted as appropriate until both are functioning as
required - see richardreeve/Diversity.jl#18.

However, the other important feature that it holds is to allow an
interface to R, allowing any existing R functionality to be carried
out on julia trees, and trees to be read from disk and written using R
helper functions. This is not built into the package as it will make
RCall (and R) a dependency, which I wanted to avoid. Instead, for now,
if you want to use the R interface you need to do it manually, as
below:

```julia
julia> using RCall

julia> cd(Pkg.dir("Phylo", "src"))

julia> include("rcall.jl")

R> library(ape)
```

You can then translate back and forth using `NamedTree` contructors on
R `phylo` objects, and `RObject` constructors on julia `NamedTree` types:

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

julia> rt2 = RObject(jt)
RCall.RObject{RCall.VecSxp}

Phylogenetic tree with 10 tips and 9 internal nodes.

Tip labels:
	t10, t8, t1, t2, t6, t5, ...

Rooted; includes branch lengths.

julia> @rput rt;

julia> @rput rt2;

R> all.equal(rt, rt2)
[1] TRUE
```

For the time being the code will only work with rooted binary trees
with named tips and branch lengths. If there's [demand][issues-url]
for other types of trees, I'll look into it, but this is all I
currently need.

## Install

*Phylo* is in `METADATA` so can be installed via `Pkg.add("Phylo")`.



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

[issues-url]: https://github.com/richardreeve/Phylo.jl/issues

[pkg-0.5-img]: http://pkg.julialang.org/badges/Phylo_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/?pkg=Phylo&ver=0.5
[pkg-0.6-img]: http://pkg.julialang.org/badges/Phylo_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=Phylo&ver=0.6

[issues-url]: https://github.com/richardreeve/Phylo.jl/issues
[diversity-url]: https://github.com/richardreeve/Diversity.jl/
[phylogenetics-url]: https://github.com/richardreeve/Diversity.jl/tree/phylogenetics
