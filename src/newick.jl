using Tokenize
using Tokenize.Lexers
using Missings
const T = Tokenize.Tokens

noname(name::String) = name == ""

tokenerror(token, str) = error("Unexpected $(token.kind) token " *
                               "'$(untokenize(token))' instead of '$str'")

isWHITESPACE(token) = token.kind == T.WHITESPACE
isEQ(token) = token.kind == T.EQ
isEQorRSQUARE(token) = token.kind ∈ [T.EQ, T.RSQUARE]
isIDENTIFIER(token, text) =
    (token.kind == T.IDENTIFIER) & (lowercase(untokenize(token)) == text)
    isBEGIN(token) = (token.kind == T.BEGIN) | isIDENTIFIER(token, "begin")
isTAXA(token) = isIDENTIFIER(token, "taxa")
isTREE(token) = isIDENTIFIER(token, "tree")
isTREES(token) = isIDENTIFIER(token, "trees")
isDIMENSIONS(token) = isIDENTIFIER(token, "dimensions")
isTAXLABELS(token) = isIDENTIFIER(token, "taxlabels")
isTRANSLATE(token) = isIDENTIFIER(token, "translate")
isEND(token) = (token.kind == T.END) | isIDENTIFIER(token, "end") | isIDENTIFIER(token, "endblock")

function iterateskip(tokens, state = nothing)
    result = (state === nothing) ? iterate(tokens) : iterate(tokens, state)
    result === nothing && return nothing
    token, state = result
    while isWHITESPACE(token)
        result = iterate(tokens, state)
        result === nothing && return nothing
        token, state = result
    end
    ut = untokenize(token)
    if length(ut) > 1
        if ut[length(ut)] == '\r'
            token = T.Token(token.kind, token.startpos, token.endpos,
                            token.startbyte, token.endbyte,
                            token.val[1:(length(token.val)-1)],
                            token.token_error, token.dotop, token.suffix)
        end
    end
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
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
    end
    return token, state, join(sofar)
end

function checktosemi(test::Function, token, state, tokens)
    if !test(token)
        return false
    end
    result = iterateskip(tokens, state)
    result === nothing && return nothing
    token, state = result
    if token.kind != T.SEMICOLON
        tokenerror(token, ";")
    end
    return true
end

function parsevector(token, state, tokens, ::Type{TY}, sgn) where TY <: Real
    vec = TY[]
    while token.kind != T.RBRACE && token.kind != T.ENDMARKER
        if token.kind == T.PLUS
            sgn = +;
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
        elseif token.kind == T.MINUS
            sgn = -;
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
        end
        push!(vec, sgn(parse(TY, untokenize(token))))
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
        if token.kind != T.COMMA && token.kind != T.RBRACE
            tokenerror(token, ",")
        elseif token.kind == T.COMMA
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
        end
    end

    return token, state, vec
end

function parsevector(token, state, tokens)
    result = iterateskip(tokens, state)
    result === nothing && return nothing
    token, state = result

    sgn = +;
    if token.kind == T.MINUS
        sgn = -;
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
    elseif token.kind == T.PLUS
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
    end

    if token.kind == T.INTEGER
        return parsevector(token, state, tokens, Int, sgn)
    elseif token.kind == T.FLOAT
        return parsevector(token, state, tokens, Float64, sgn)
    end

    vec = String[]
    while token.kind != T.RBRACE && token.kind != T.ENDMARKER
        if token.kind ∈ [T.STRING, T.CHAR]
            push!(vec, untokenize(token)[2:end-1])
        else
            push!(vec, untokenize(token))
        end
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
        if token.kind != T.COMMA && token.kind != T.RBRACE
            tokenerror(token, ",")
        elseif token.kind == T.COMMA
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
        end
    end

    return token, state, vec
end

function parsedict(token, state, tokens)
    dict = Dict{String, Any}()
    result = iterateskip(tokens, state)
    result === nothing && return nothing
    token, state = result
    if token.kind != T.AND
        tokenerror(token, "&")
    else
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
    end

    while token.kind != T.RSQUARE && token.kind != T.ENDMARKER
        token, state, key = tokensgetkey(token, state, tokens, isEQorRSQUARE)
        if token.kind != T.RSQUARE # Allow [&R] as a valid (empty) dict
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
            if token.kind == T.LBRACE
                token, state, value = parsevector(token, state, tokens)
            else
                sgn = +;
                if token.kind == T.PLUS
                    sgn = +;
                    result = iterateskip(tokens, state)
                    result === nothing && return nothing
                    token, state = result
                elseif token.kind == T.MINUS
                    sgn = -;
                    result = iterateskip(tokens, state)
                    result === nothing && return nothing
                    token, state = result
                end

                if token.kind == T.INTEGER
                    value = sgn(parse(Int, untokenize(token)))
                elseif token.kind == T.FLOAT
                    value = sgn(parse(Float64, untokenize(token)))
                elseif token.kind ∈ [T.STRING, T.CHAR]
                    value = untokenize(token)[2:end-1]
                else
                    value = untokenize(token)
                end
            end
            dict[key] = value

            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
            if token.kind != T.COMMA && token.kind != T.RSQUARE
                tokenerror(token, ", or ]")
            elseif token.kind == T.COMMA
                result = iterateskip(tokens, state)
                result === nothing && return nothing
                token, state = result
            end
        end
    end

    if token.kind == T.RSQUARE
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
    end

    return token, state, dict
end

function parsenode(token, state, tokens, tree::TREE,
                   lookup, siblings, istip) where
    {RT, NL, N, B, TREE <: AbstractTree{OneTree, RT, NL, N, B}}
    myname = missing
    endkinds = [T.SEMICOLON, T.COLON, T.COMMA, T.RPAREN, T.LSQUARE]
    if token.kind ∉ endkinds # We got a nodename
        token, state, name = tokensgetkey(token, state, tokens,
                                          t -> t.kind ∈ endkinds)
        if isempty(lookup)
            myname = createnode!(tree, name)
        else
            if istip
                if haskey(lookup, name)
                    myname = lookup[name]
                else
                    @warn "Wrongly named tip '$name' found in tree " *
                        "with named tips, adding"
                    myname = createnode!(tree, name)
                end
            else
                if haskey(lookup, name)
                    myname = lookup[name]
                else
                    myname = createnode!(tree, name)
                end
            end
        end
    else
        if istip && !isempty(lookup)
            @warn "Anonymous tip found in tree with named tips"
        end
        myname = createnode!(tree)
    end

    foundcolon = false
    if token.kind == T.COLON
        foundcolon = true
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
    end

    siblings[getnodename(tree, myname)] = Dict{String, Any}()
    if token.kind == T.LSQUARE
        token, state, dict = parsedict(token, state, tokens)
        siblings[getnodename(tree, myname)]["dict"] = dict
        setnodedata!(tree, myname, dict)
    end

    if token.kind == T.COLON || foundcolon
        if token.kind == T.COLON
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
        end
        sgn = +;
        if token.kind == T.PLUS
            sgn = +;
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
        elseif token.kind == T.MINUS
            sgn = -;
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
        end
        if token.kind ∈ [T.INTEGER, T.FLOAT]
            siblings[getnodename(tree, myname)]["length"] =
                sgn(parse(Float64, untokenize(token)))
        else
            tokenerror(token, "a length")
        end
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
    end

    return token, state, myname
end

function parsenewick!(token, state, tokens, tree::TREE,
                     lookup = Dict(), depth = 0,
                     children = Dict{NL, Dict{String, Any}}()) where
    {RT, NL, N, B, TREE <: AbstractTree{OneTree, RT, NL, N, B}}
    mychildren = Dict{NL, Dict{String, Any}}()
    while (token.kind != T.RPAREN) & (token.kind != T.ENDMARKER)
        result = iterateskip(tokens, state)
        result === nothing &&
            error("Tree ended at depth $depth before right bracket")
        token, state = result
        if token.kind == T.LPAREN
            token, state = parsenewick!(token, state, tokens, tree,
                                        lookup, depth + 1, mychildren)
        else
            token, state, nodename = parsenode(token, state, tokens, tree,
                                               lookup, mychildren, true)
        end
    end
    result = iterateskip(tokens, state)
    result === nothing &&
        error("Tree ended at depth $depth before" *
              (depth > 0 ? "right bracket" : "semicolon"))
    token, state = result
    token, state, nodename = parsenode(token, state, tokens, tree,
                                       lookup, children, false)
    for child in keys(mychildren)
        dict = mychildren[child]
        haskey(dict, "length") ?
            createbranch!(tree, nodename, child, dict["length"]) :
            createbranch!(tree, nodename, child)
    end

    if depth == 0
        # Should be at end of tree
        if token.kind == T.SEMICOLON
            # Am at end of tree
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
        else
            error("At end of tree, but not ';'")
        end
        if !validate!(tree)
            error("Tree $(gettreename(tree)) does not validate!")
        end
    end

    return token, state
end

function parsenewick(tokens::Tokenize.Lexers.Lexer, ::Type{TREE}) where
    {RT, N, B, TREE <: AbstractTree{OneTree, RT, String, N, B}}
    result = iterateskip(tokens)
    if result === nothing
        error("Unexpected end of file at start of newick file")
    end
    token, state = result
    if token.kind == T.LPAREN
        tree = TREE()
        result = parsenewick!(token, state, tokens, tree)
        (result !== nothing) && (result[1].kind != T.ENDMARKER) &&
            @warn "Tree ended but not finished - " *
                  "got $(result[1].kind) $(untokenize(result[1]))"
        return tree
    else
        error("Unexpected $(token.kind) token '$(untokenize(token))' " *
              "at start of newick file")
    end
end

"""
    parsenewick(io::IOBuffer, ::Type{TREE}) where TREE <: AbstractTree

Parse an IOBuffer containing a newick tree and convert into a phylogenetic
tree of type TREE <: AbstractTree
"""
parsenewick(io::IOBuffer, ::Type{TREE}) where TREE <: AbstractTree =
    parsenewick(tokenize(io), TREE)

"""
    parsenewick(io::String, ::Type{TREE}) where TREE <: AbstractTree

Parse a String containing a newick tree and convert into a phylogenetic
tree of type TREE <: AbstractTree
"""
parsenewick(s::String, ::Type{TREE}) where TREE <: AbstractTree =
    parsenewick(IOBuffer(s), TREE)

"""
    parsenewick(io::IOStream, ::Type{TREE}) where TREE <: AbstractTree

Parse an IOStream containing a newick tree and convert into a phylogenetic
tree of type TREE <: AbstractTree
"""
function parsenewick(ios::IOStream, ::Type{TREE}) where TREE <: AbstractTree
    buf = IOBuffer()
    print(buf, read(ios, String))
    seek(buf, 0)
    return parsenewick(buf, TREE)
end

"""
    parsenewick(inp)

Parse some input containing a newick tree and convert into a phylogenetic
tree of type RootedTree
"""
parsenewick(inp) = parsenewick(inp, RootedTree)

function parsetaxa(token, state, tokens, taxa)
    if !isDIMENSIONS(token)
        tokenerror(token, "Dimensions")
    end
    result = iterateskip(tokens, state)
    result === nothing && return nothing
    token, state = result
    token, state, ntaxstr = tokensgetkey(token, state, tokens)
    if lowercase(ntaxstr) != "ntax"
        error("Unexpected label '$ntaxstr=' not 'ntax='")
    end
    result = iterateskip(tokens, state)
    result === nothing && return nothing
    token, state = result
    ntax = parse(Int64, untokenize(token))
    result = iterateskip(tokens, state)
    result === nothing && return nothing
    token, state = result
    if token.kind != T.SEMICOLON
        tokenerror(token, ";")
    end
    result = iterateskip(tokens, state)
    result === nothing && return nothing
    token, state = result
    if !isTAXLABELS(token)
        tokenerror(token, "Taxlabels")
    end
    result = iterateskip(tokens, state)
    result === nothing && return nothing
    token, state = result
    while token.kind != T.SEMICOLON && token.kind != T.ENDMARKER
        name = untokenize(token)
        if token.kind ∈ [T.STRING, T.CHAR]
            name = name[2:end-1]
        end
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
        taxa[name] = name
        if token.kind == T.LSQUARE
            while token.kind != T.RSQUARE
                result = iterateskip(tokens, state)
                result === nothing && return nothing
                token, state = result
            end
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
        end
    end
    if length(taxa) != ntax
        @warn "Taxa list length ($(length(taxa))) and ntax ($ntax) do not match"
    end

    result = iterateskip(tokens, state)
    result === nothing && return nothing
    token, state = result
    if !checktosemi(isEND, token, state, tokens)
        tokenerror(token, "End;")
    end
    return iterateskip(tokens, state)
end

function parsetrees(token, state, tokens, ::Type{TREE}, taxa) where
    {RT, N, B, TREE <: AbstractTree{OneTree, RT, String, N, B}}
    notaxa = isempty(taxa)
    if isTRANSLATE(token)
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
        while token.kind != T.SEMICOLON && token.kind != T.ENDMARKER
            short = untokenize(token)
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
            proper = untokenize(token)
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
            if haskey(taxa, proper)
                delete!(taxa, proper)
                taxa[short] = proper
            elseif notaxa
                @warn "Found a '$short' => '$proper' link without a taxa entry first"
                taxa[short] = proper
            else
                @warn "Missing '$proper' in taxa block, but in 'trees -> translate' block"
            end
            if token.kind == T.COMMA
                result = iterateskip(tokens, state)
                result === nothing && return nothing
                token, state = result
            end
        end
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
    end

    trees = Dict{String, TREE}()
    treedata = Dict{String, Dict{String, Any}}()
    numTrees = 0
    while isTREE(token)
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
        token, state, treename = tokensgetkey(token, state, tokens,
                                              t -> t.kind ∈ [T.LSQUARE, T.EQ])
        numTrees += 1
        if numTrees < 10
            @info "Created a tree called \"$treename\""
        elseif numTrees == 10
            @info "[Stopping reporting on tree creation]"
        end
        trees[treename] = TREE()
        createnodes!(trees[treename], collect(values(taxa)))
        if token.kind == T.LSQUARE
            token, state, treedata[treename] = parsedict(token, state, tokens)
        end
        if !isEQ(token)
            tokenerror(token, "=")
        else
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
            if token.kind == T.LSQUARE
                token, state, _ = parsedict(token, state, tokens)
            end
        end
        if token.kind != T.LPAREN
            tokenerror(token, "(")
        else
            result = parsenewick!(token, state, tokens, trees[treename], taxa)
            result === nothing && return nothing
            token, state = result
        end
    end
    if !checktosemi(isEND, token, state, tokens)
        tokenerror(token, "End;")
    end

    result = iterateskip(tokens, state)
    result === nothing && return nothing
    token, state = result
    return token, state, trees, treedata
end

function parsenexus(token, state, tokens, ::Type{TREE}) where
    {RT, NL, N, B, TREE <: AbstractTree{OneTree, RT, NL, N, B}}
    trees = missing
    treedata = missing
    taxa = Dict{NL, NL}()
    while isBEGIN(token) || token.kind == T.LSQUARE
        if token.kind == T.LSQUARE
            while token.kind != T.RSQUARE
                result = iterateskip(tokens, state)
                result === nothing && return nothing
                token, state = result
            end
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
        else
            result = iterateskip(tokens, state)
            result === nothing && return nothing
            token, state = result
            if checktosemi(isTAXA, token, state, tokens)
                result = iterateskip(tokens, state)
                result === nothing && return nothing
                token, state = result
                token, state = parsetaxa(token, state, tokens, taxa)
            elseif checktosemi(isTREES, token, state, tokens)
                result = iterateskip(tokens, state)
                result === nothing && return nothing
                token, state = result
                token, state, trees, treedata = parsetrees(token, state, tokens, TREE, taxa)
            else
                @warn "Unexpected nexus block '$(untokenize(token))', skipping..."
                result = iterateskip(tokens, state)
                result === nothing && return nothing
                token, state = result
                while !checktosemi(isEND, token, state, tokens) && token.kind != T.ENDMARKER
                    result = iterateskip(tokens, state)
                    result === nothing && return nothing
                    token, state = result
                end
                result = iterateskip(tokens, state)
                result === nothing && return nothing
                token, state = result
            end
        end
    end

    if token.kind != T.ENDMARKER
        tokenerror(token, "end of file")
    end
    return TreeSet(trees, treedata)
end

function parsenexus(tokens::Tokenize.Lexers.Lexer,
                    ::Type{TREE}) where {RT, NL, N, B,
                                         TREE <: AbstractTree{OneTree, RT, NL, N, B}}
    result = iterateskip(tokens)
    result === nothing && return nothing
    token, state = result
    if token.kind == T.COMMENT && lowercase(untokenize(token)) == "#nexus"
        result = iterateskip(tokens, state)
        result === nothing && return nothing
        token, state = result
        return parsenexus(token, state, tokens, TREE)
    else
        error("Unexpected $(token.kind) token '$(untokenize(token))' " *
              "at start of nexus file")
    end
end

"""
    parsenexus(io::IOBuffer, ::Type{TREE}) where TREE <: AbstractTree

Parse an IOBuffer containing a nexus tree and convert into a phylogenetic
tree of type TREE <: AbstractTree
"""
parsenexus(io::IOBuffer, ::Type{TREE}) where TREE <: AbstractTree =
    parsenexus(tokenize(io), TREE)

"""
    parsenexus(io::IOStream, ::Type{TREE}) where TREE <: AbstractTree

Parse an IOStream containing a nexus tree and convert into a phylogenetic
tree of type TREE <: AbstractTree
"""
function parsenexus(ios::IOStream, ::Type{TREE}) where TREE <: AbstractTree
    buf = IOBuffer()
    print(buf, read(ios, String))
    seek(buf, 0)
    return parsenexus(buf, TREE)
end

"""
    parsenexus(inp)

Parse some input containing a nexus tree and convert into a phylogenetic
tree of type RootedTree
"""
parsenexus(inp) = parsenexus(inp, RootedTree)
