# Phylo.API

The **Phylo.API** submodule provides the API that must be extended
for new `AbstractTree`, `AbstractNode` and `AbstractBranch` subtypes.

## Usage

Providing additional code to extend the functionality of the system is simple:

```julia
using Phylo

struct SimplestTree <: AbstractTree{Int, Int}
    nodes::OrderedDict{Int, BinaryNode{Int}}
    branches::Dict{Int, Branch{Int}}
end

import Phylo.API: _addnode!
function _addnode!(tree::SimplestTree, num)
    _setnode!(tree, num, BinaryNode{Int}())
    return num
end
```

creates a new `SimplestTree` type (a subtype of `AbstractTree`) and
extends `Phylo.API._addnode!()` (and therefore the directly accessible
`addnode!()` interface) to handle the `SimplestTree` subtype of
`AbstractTree`. See docs here to see which `Phylo.API` functions have
to be extended for any new subtype, and which have default
implementations.

```@autodocs
Modules = [Phylo.API]
```

```@index
```
