using Phylo
using Tokenize
using Tokenize.Lexers
const T = Tokenize.Tokens

function parsenewick(io::IOStream)
    tree = NamedTree()
    children = Dict{Int,Vector{String}}()
    children[0] = []
    lengths = Dict{String,Float64}()
    depth = 0
    makenode = true
    addlength = false
    currentname = Nullable{String}()
    readyforparent = false
    eot = false
    insquarebrackets = false
    for token in tokenize(io) 
        processed = false

        if token.kind == T.LSQUARE
            insquarebrackets = true
            processed = true
        elseif insquarebrackets
            if token.kind == T.RSQUARE
                insquarebrackets = false
            end
            processed = true
        else
            if token.kind == T.LPAREN # Going to make some children
                depth += 1
                children[depth] = []
                makenode = true
                currentname = Nullable{String}()
                processed = true
            end

            if (makenode || readyforparent) &&
                token.kind ∈ [T.COMMA, T.RPAREN] # We got an anonymous node
                if readyforparent
                    parent = addnode!(tree)
                    makenode = false
                    for child in children[depth]
                        haskey(lengths, child) ?
                            addbranch!(tree, parent, child, lengths[child]) :
                            addbranch!(tree, parent, child)
                    end
                    children[depth] = []
                    depth -= 1
                    push!(children[depth], parent)
                    currentname = Nullable(parent)
                    readyforparent = false
                else
                    nodename = addnode!(tree)
                    push!(children[depth], nodename)
                    currentname = Nullable{String}()
                    makenode = false
                end
            end

            if token.kind == T.COMMA # Ready for a new node
                makenode = true
                addlength = false
                processed = true
            end
            
            if token.kind == T.RPAREN # Add a parent node and connect to children
                readyforparent = true
                processed = true
            end
            
            if token.kind == T.COLON # Length coming
                addlength = true
                if readyforparent
                    parent = addnode!(tree)
                    makenode = false
                    for child in children[depth]
                        haskey(lengths, child) ?
                            addbranch!(tree, parent, child, lengths[child]) :
                            addbranch!(tree, parent, child)
                    end
                    children[depth] = []
                    depth -= 1
                    push!(children[depth], parent)
                    currentname = Nullable(parent)
                    readyforparent = false
                elseif makenode # We have an anonymous leaf
                    nodename = addnode!(tree)
                    push!(children[depth], nodename)
                    currentname = Nullable(nodename)
                    makenode = false
                end
                processed = true
            end
            
            if token.kind ∈ [T.STRING, T.IDENTIFIER, # Very likely a leaf
                             T.INTEGER, T.FLOAT] # Maybe a leaf
                if token.kind ∈ [T.INTEGER, T.FLOAT] && addlength # It's a length!
                    !isnull(currentname) ||
                        error("Found length $(untokenize(token)) " *
                              "with no associated node")
                    lengths[get(currentname)] = parse(Float64, untokenize(token))
                    addlength = false
                else # It's a leaf
                    if readyforparent
                        parent = addnode!(tree, untokenize(token))
                        makenode = false
                        for child in children[depth]
                            haskey(lengths, child) ?
                                addbranch!(tree, parent, child, lengths[child]) :
                                addbranch!(tree, parent, child)
                        end
                        children[depth] = []
                        depth -= 1
                        push!(children[depth], parent)
                        currentname = Nullable(parent)
                        readyforparent = false
                    else
                        makenode ||
                            error("Found nodename '$(untokenize(token))' " *
                                  "while not expecting one")
                        nodename = addnode!(tree, untokenize(token))
                        push!(children[depth], nodename)
                        currentname = Nullable(nodename)
                        makenode = false
                    end
                end
                processed = true
            end

            if token.kind == T.SEMICOLON # End of tree
                if readyforparent
                    parent = addnode!(tree)
                    makenode = false
                    for child in children[depth]
                        haskey(lengths, child) ?
                            addbranch!(tree, parent, child, lengths[child]) :
                            addbranch!(tree, parent, child)
                    end
                    children[depth] = []
                    depth -= 1
                    push!(children[depth], parent)
                    currentname = Nullable(parent)
                    readyforparent = false
                elseif makenode
                    nodename = addnode!(tree)
                    push!(children[depth], nodename)
                    currentname = Nullable(nodename)
                    makenode = false
                end
                eot = true
                processed = true
            end
end

if token.kind == T.ENDMARKER # End of tokens
    eot ||
        error("Reached end of tokenbs without reaching end of tree")
    processed = true
end

if !processed && token.kind != T.WHITESPACE
    info("Discarding unused token '$(untokenize(token))'")
end
end

tree = resetleaves(tree)
return tree
end
