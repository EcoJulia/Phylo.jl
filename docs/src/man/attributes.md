# Getting tree attributes

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
gettreeinfo
validate!
invalidate!
branchdims
treetype
treesettype
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
getsiblings
```

## Methods on Branches

```@docs
src
dst
conn
conns
haslength
getlength
hasrootheight
getrootheight
setrootheight!
getbranchdata
setbranchdata!
```
