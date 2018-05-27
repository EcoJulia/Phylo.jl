using Compat.@warn
using Compat.@info
using Tokenize
using Tokenize.Lexers
const T = Tokenize.Tokens

noname(name::String) = name == ""

tokenerror(token, str) = error("Unexpected $(token.kind) token " *
                               "'$(untokenize(token))' instead of '$str'")

isWHITESPACE(token) = token.kind == T.WHITESPACE
isEQ(token) = token.kind == T.EQ
isEQorRSQUARE(token) = token.kind ∈ [T.EQ, T.RSQUARE]
isIDENTIFIER(token, text) =
    (token.kind == T.IDENTIFIER) && (lowercase(untokenize(token)) == text)
isBEGIN(token) = isIDENTIFIER(token, "begin")
isTAXA(token) = isIDENTIFIER(token, "taxa")
isTREE(token) = isIDENTIFIER(token, "tree")
isTREES(token) = isIDENTIFIER(token, "trees")
isDIMENSIONS(token) = isIDENTIFIER(token, "dimensions")
isTAXLABELS(token) = isIDENTIFIER(token, "taxlabels")
isTRANSLATE(token) = isIDENTIFIER(token, "translate")
isEND(token) = isIDENTIFIER(token, "end")

function nextskip(tokens, state)
    token, state = next(tokens, state)
    while isWHITESPACE(token)
        token, state = next(tokens, state)
    end
    #@info untokenize(token)
    return token, state
end

function tokensgetkey(token, state, tokens, finished::Function = isEQ)
    sofar = String[]
    while !finished(token) && token.kind != T.ENDMARKER
        if token.kind ∈ [T.STRING, T.CHAR]
            push!(sofar, untokenize(token)[2:end-1])
        else
            push!(sofar, untokenize(token))
        end
        token, state = nextskip(tokens, state)
    end
    return token, state, join(sofar)
end

function checktosemi(test::Function, token, state, tokens)
    if !test(token)
        return false
    end
    token, state = nextskip(tokens, state)
    if token.kind != T.SEMICOLON
        tokenerror(token, ";")
    end
    return true
end

function parsevector(token, state, tokens, ::Type{TY}) where TY <: Real
    vec = TY[]
    while token.kind != T.RBRACE && token.kind != T.ENDMARKER
        push!(vec, parse(TY, untokenize(token)))
        token, state = nextskip(tokens, state)
        if token.kind != T.COMMA && token.kind != T.RBRACE
            tokenerror(token, ",")
        elseif token.kind == T.COMMA
            token, state = nextskip(tokens, state)
        end
    end
    return token, state, vec
end

function parsevector(token, state, tokens)
    token, state = nextskip(tokens, state)

    if token.kind == T.INTEGER
        return parsevector(token, state, tokens, Int)
    elseif token.kind == T.FLOAT
        return parsevector(token, state, tokens, Float64)
    end

    vec = String[]
    while token.kind != T.RBRACE && token.kind != T.ENDMARKER
        if token.kind ∈ [T.STRING, T.CHAR]
            push!(vec, untokenize(token)[2:end-1])
        else
            push!(vec, untokenize(token))
        end
        token, state = nextskip(tokens, state)
        if token.kind != T.COMMA && token.kind != T.RBRACE
            tokenerror(token, ",")
        elseif token.kind == T.COMMA
            token, state = nextskip(tokens, state)
        end
    end
    return token, state, vec
end

function parsedict(token, state, tokens)
    dict = Dict{String, Any}()
    token, state = next(tokens, state)
    if token.kind != T.AND
        tokenerror(token, "&")
    end
    token, state = nextskip(tokens, state)
    while token.kind != T.RSQUARE && token.kind != T.ENDMARKER
        token, state, key = tokensgetkey(token, state, tokens, isEQorRSQUARE)
        if token.kind != T.RSQUARE # Allow [&R] as a valid (empty) dict
            token, state = nextskip(tokens, state)
            if token.kind == T.LBRACE
                token, state, value = parsevector(token, state, tokens)
            else
                if token.kind == T.INTEGER
                    value = parse(Int, untokenize(token))
                elseif token.kind == T.FLOAT
                    value = parse(Float64, untokenize(token))
                elseif token.kind ∈ [T.STRING, T.CHAR]
                    value = untokenize(token)[2:end-1]
                else
                    value = untokenize(token)
                end
            end
            dict[key]=value
            token, state = nextskip(tokens, state)
            if token.kind != T.COMMA && token.kind != T.RSQUARE
                tokenerror(token, ",")
            elseif token.kind == T.COMMA
                token, state = nextskip(tokens, state)
            end
        end
    end

    return token, state, dict
end

function parsenode(token, state, tokens, tree::TREE,
                   lookup, siblings, istip) where {NL, BL, TREE <: AbstractBranchTree{NL, BL}}
    myname = ""
    endkinds = [T.SEMICOLON, T.COLON, T.COMMA, T.RPAREN, T.LSQUARE]
    if token.kind ∉ endkinds
        # We got a nodename
        token, state, name = tokensgetkey(token, state, tokens,
                                          t -> t.kind ∈ endkinds)
        if isempty(lookup)
            myname = addnode!(tree, name)
            @info "Created new named $(istip?"tip":"node") '$myname'"
        else
            if istip
                if haskey(lookup, name)
                    myname = lookup[name]
                else
                    @warn "Wrongly named tip '$name' found in tree " *
                        "with named tips, adding"
                    myname = addnode!(tree, name)
                end
                @info "Created lookup named tip '$myname'"
            else
                if haskey(lookup, name)
                    @info "Recognized internal node '$name' -> " *
                        "'$(lookup[name])' found in tree with named tips"
                    myname = lookup[name]
                else
                    @info "Unrecognized internal node '$name' " *
                        "found in tree with named tips"
                    myname = addnode!(tree, name)
                end
                @info "Created lookup named internal node '$myname'"
            end
        end
    else
        if istip && !isempty(lookup)
            @warn "Anonymous tip found in tree with named tips"
        end
        myname = addnode!(tree)
        @info "Created anonymous $(istip?"tip":"node") '$myname'"
    end
    siblings[myname] = Dict{String, Any}()
    if token.kind == T.LSQUARE
        token, state, dict = parsedict(token, state, tokens)
        siblings[myname]["dict"] = dict
        setnoderecord!(tree, myname, dict)
        token, state = nextskip(tokens, state)
    end
    if token.kind == T.COLON
        token, state = nextskip(tokens, state)
        if token.kind ∈ [T.INTEGER, T.FLOAT]
            siblings[myname]["length"] = parse(Float64, untokenize(token))
        else
            tokenerror(token, "a length")
        end
        token, state = nextskip(tokens, state)
    end
    return token, state, myname
end

function parsenewick(token, state, tokens, tree::TREE,
                     lookup = Dict(), depth = 0,
                     children = Dict{NL, Dict{String, Any}}()) where
    {NL, BL, TREE <: AbstractBranchTree{NL, BL}}
    @warn "Entering depth $depth, $(token.kind) = '$(untokenize(token))'"
    token, state = nextskip(tokens, state)
    mychildren = Dict{NL, Dict{String, Any}}()
    while token.kind != T.RPAREN && token.kind != T.ENDMARKER
        if token.kind == T.LPAREN
            token, state = parsenewick(token, state, tokens, tree,
                                       lookup, depth + 1, mychildren)
        elseif token.kind == T.COMMA
            token, state = nextskip(tokens, state)
        else
            token, state, nodename = parsenode(token, state, tokens, tree,
                                               lookup, mychildren, true)
            @warn "Added a tip called $nodename"
        end
    end
    token, state = nextskip(tokens, state)
    token, state, nodename = parsenode(token, state, tokens, tree,
                                       lookup, children, false)
    @warn "Added a node called $nodename"
    for child in keys(mychildren)
        dict = mychildren[child]
        haskey(dict, "length") ?
            addbranch!(tree, nodename, child, dict["length"]) :
            addbranch!(tree, nodename, child)
        @warn "Added branch from $nodename to $child"
    end
    @warn "Exiting depth $depth, $(token.kind) = '$(untokenize(token))'"
    if depth == 0
        # Should be at end of tree
        if token.kind == T.SEMICOLON
            # Am at end of tree
            token, state = nextskip(tokens, state)
        else
            error("At end of tree, but not ';'")
        end
        tree = resetleaves!(tree)
    end

    return token, state
end

function parsenewick(tokens::Tokenize.Lexers.Lexer,
                     ::Type{TREE}) where TREE <: AbstractBranchTree{String, Int}
    token, state = nextskip(tokens, start(tokens))
    if token.kind == T.LPAREN
        tree = TREE()
        parsenewick(token, state, tokens, tree)
        return tree
    else
        error("Unexpected $(token.kind) token '$(untokenize(token))' " *
              "at start of newick file")
    end
end

parsenewick(io::IOBuffer, ::Type{TREE}) where TREE <: AbstractBranchTree =
    parsenewick(tokenize(io), TREE)

parsenewick(s::String, ::Type{TREE}) where TREE <: AbstractBranchTree =
    parsenewick(IOBuffer(s), TREE)


function parsenewick(ios::IOStream, ::Type{TREE}) where TREE <: AbstractBranchTree
    buf = IOBuffer()
    print(buf, read(ios, String))
    seek(buf, 0)
    return parsenewick(buf, TREE)
end

parsenewick(inp) = parsenewick(inp, NamedTree)
