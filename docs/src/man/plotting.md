# Plotting phylogenetic trees

Phylo defines recipes for all `AbstractTree`s, allowing them to be plotted with
 [Plots.jl](https://docs.juliaplots.org/latest).

### Keywords

It adds these keywords to the ones initiallysupported by Plots.jl:
- treetype: choosing `:fan` or `:dendrogram` determines the shape of the tree
- marker_group: applies the `group` keyword to node markers
- line_group: applies the `group` keyword to branch lines
- showtips: `true` (the default) shows the leaf names
- tipfont: a tuple defining the font to use for leaf names (default is `(7,))

### Example plots
For this example, we will use the phylogeny of all extant hummingbird species.

Read the tree and plot it using two different `treetype`s. The 
default is `:dendrogram`
```@example plotting
using Phylo, Plots
ENV["GKSwstype"]="nul" # hide
default(linecolor = :black) # looks nicer with black lines
hummers = open(parsenewick, "../../../test/hummingbirds.tree")
plot(hummers, size = (400, 600), showtips = false)
```

For larger trees, the `:fan` treetype may work better
```@example plotting
plot(hummers, treetype = :fan)
```

### Sorting trees for plotting
Many phylogenies look more aesthetically pleasing if the descendants from each
node are sorted in order of their size. This is called `ladderize` in some other
packages.
```@example plotting
sort!(hummmers, rev = true)
plot(hummers, treetype = :fan)
```

### Coloring branches or nodes by a variable
It is common in evolutionary studies to color the branches or node markers with
the value of some variable. Plots already offers the keyword attributes 
`marker_z` and `line_z` for these uses, and they also work on Phylo objects.

We can pass either 
- a `Vector` with the same number of elements as there are 
branches / internal nodes, where the values follow a depthfirst order
(because the tree is plotted in depthfirst order);
- a `Dict` of `node => value`, with the value to be plotted for each node 
(skipping nodes not in the Dict).

To demonstrate, let' start by defining a custom function for evolving a trait
on the phylogeny accoding to Brownian motion, using the utility funciton 
`map_depthfirst`
```@example plotting
evolve(tree) = map_depthfirst((val, node) -> val + randn(), 0., tree, Float64)
trait = evolve(hummers)
plot(hummers, treetype = :fan, line_z = trait, linecolor = :RdYlBu, linewidth = 5, showtips = false)
```

The inbuilt facilities for sampling traits on trees on Phylo returns a 
`Node => value` Dict, which can also be passed to `marker_z`
```@example plotting
brownsampler = BrownianTrait(hummers, "Trait")
plot(hummers, 
     showtips = false, marker_z = rand(brownsampler), 
     linewidth = 2, markercolor = :RdYlBu)
```


We can also use the `map_depthfirst` utility function to highlight the clade
descending from, e.g., Node 248. Here, the recursive function creates a vector 
of colors which will all be orange after encountering Node 248 but black before
```@example plotting
clade_color = map_depthfirst((val, node) -> node == "Node 248" ? :orange : val, :black, hummers)
plot(hummers, linecolor = clade_color, showtips = false, linewidth = 2)
```

The Plots attributes related to markers (`markersize`, `markershape` etc.) will
put markers on the internal nodes (or on both internal and tip nodes if a longer)
vector is passed). The `series_attributes` keyword is also supported and behaves
the same way
```@example plotting
plot(hummers, 
    size = (400, 800), 
    linecolor = :orange, linewidth = 5, 
    markersize = 10, markercolor = :steelblue, markerstrokecolor = :white,
    series_annotations = text.(1:nnodes(hummers), 6, :center, :center, :white)
)
```

The `markergroup` and `linegroup` keywords allow plotting discrete values onto
nodes or branches within the phylogeny
```@example plotting
groups = "Group" .* rand(3, ninternal(hummers))
plot(hummers, treetype = :fan, markergroup = groups)
```

``` @example plotting
# @enum TemperatureTrait lowTempPref midTempPref highTempPref # hide
# tempsampler = SymmetricDiscreteTrait(hummers, TemperatureTrait, 0.4) # hide
# plot(hummers, showtips = false,   # hide
#    marker_group = rand(tempsampler),  # hide
#    legend = :topleft, msc = :white, treetype = :fan,  # hide
#    c = [:red :blue :green]) # hide

```@docs
map_depthfirst
sort
sort!
```

