using Phylo
using Tokenize
using Tokenize.Lexers
const T = Tokenize.Tokens

noname(name::String) = name == ""

function parsenewick(io::IOBuffer)
    tree = NamedTree()
    children = Dict{Int,Vector{String}}([0 => [], -1 => []])
    lengths = Dict{String,Float64}()
    depth = 0
    makenode = true
    addlength = false
    currentname = ""
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
                currentname = ""
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
                    currentname = parent
                    readyforparent = false
                else
                    nodename = addnode!(tree)
                    push!(children[depth], nodename)
                    currentname = ""
                    makenode = false
                end
            end

            if token.kind == T.COMMA # Ready for a new node
                makenode = true
                addlength = false
                processed = true
            elseif token.kind == T.RPAREN # Add a parent node and connect to children
                readyforparent = true
                processed = true
            elseif token.kind == T.COLON # Length coming
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
                    currentname = parent
                    readyforparent = false
                elseif makenode # We have an anonymous leaf
                    nodename = addnode!(tree)
                    push!(children[depth], nodename)
                    currentname = nodename
                    makenode = false
                end
                processed = true
            elseif token.kind ∈ [T.STRING, T.CHAR, T.IDENTIFIER, # A leaf
                                 T.INTEGER, T.FLOAT] || # A leaf or a length
                                     T.iskeyword(token.kind) ||
                                     T.isliteral(token.kind) ||
                                     T.isoperator(token.kind)
                if token.kind ∈ [T.INTEGER, T.FLOAT] && addlength # A length!
                    noname(currentname) ?
                        error("Found length $(untokenize(token)) " *
                              "with no associated node") :
                              lengths[currentname] = parse(Float64, untokenize(token))
                    addlength = false
                else # It's a leaf
                    name = token.kind == T.STRING ?
                        Meta.parse(untokenize(token)) : token.kind == T.CHAR ?
                        replace(untokenize(token), "'", "") :
                        untokenize(token)
                    if readyforparent
                        parent = addnode!(tree, name)
                        makenode = false
                        for child in children[depth]
                            haskey(lengths, child) ?
                                addbranch!(tree, parent, child, lengths[child]) :
                                addbranch!(tree, parent, child)
                        end
                        children[depth] = []
                        depth -= 1
                        push!(children[depth], parent)
                        currentname = parent
                        readyforparent = false
                    else
                        makenode ||
                            error("Found nodename '$name' " *
                                  "while not expecting one")
                        nodename = addnode!(tree, name)
                        push!(children[depth], nodename)
                        currentname = nodename
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
                    currentname = parent
                    readyforparent = false
                elseif makenode
                    nodename = addnode!(tree)
                    push!(children[depth], nodename)
                    currentname = nodename
                    makenode = false
                end
                eot = true
                processed = true
            end
end

if token.kind == T.ENDMARKER # End of tokens
    eot ||
        error("Reached end of file without reaching end of tree")
    depth == 0 ||
        error("Malformed tree ended without closing brackets")
    processed = true
end

depth >= 0 ||
    error("Malformed tree had too many closing brackets")

if !processed && token.kind != T.WHITESPACE
    warn("Discarding unused token '$(untokenize(token))'")
end
end

tree = resetleaves!(tree)
return tree
end

parsenewick(s::String) = parsenewick(IOBuffer(s))

function parsenewick(ios::IOStream)
    buf = IOBuffer()
    print(buf, read(ios, String))
    seek(buf, 0)
    return parsenewick(buf)
end
