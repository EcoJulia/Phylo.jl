# Getting tree attributes

## Phylo module and type docs

```@docs
Phylo
BinaryTree
NamedBinaryTree
BinaryNode
Node
NamedTree
PolytomousTree
NamedPolytomousTree
LinkTree
LinkBranch
LinkNode
RootedTree
ManyRootTree
UnrootedTree
TreeSet
```

## Methods on TreeSets

```@docs
ntrees
gettree
gettrees
nroots
getroots
gettreenames
```

## Methods on Trees

```@docs
mrca
nodeheights
getleafnames
getleaves
nleaves
nnodes
ninternal
nbranches
distance
distances
heighttoroot
heightstoroot
getroot
treenametype
gettreename
roottype
nodetype
nodedatatype
nodenametype
branchtype
branchdatatype
branchnametype
getnodenames
getnodename
hasnode
getnode
getnodes
getinternalnodes
getbranchnames
getbranchname
hasbranch
getbranch
getbranches
validate!
branchdims
```

## Methods on Nodes

```@docs
isleaf
isroot
isinternal
isunattached
degree
indegree
outdegree
hasinbound
getinbound
getoutbounds
getconnections
conn
conns
hasoutboundspace
hasinboundspace
getleafinfo
setleafinfo!
leafinfotype
getnodedata
setnodedata!
hasheight
getheight
setheight!
getparent
getancestors
getchildren
getdescendants
```

## Methods on Branches

```@docs
src
dst
getlength
hasrootheight
getrootheight
setrootheight!
getbranchdata
setbranchdata!
```
