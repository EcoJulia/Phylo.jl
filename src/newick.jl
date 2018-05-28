using Compat.@warn
using Compat.@info
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

function parsevector(token, state, tokens, ::Type{TY}, sgn) where TY <: Real
    vec = TY[]
    while token.kind != T.RBRACE && token.kind != T.ENDMARKER
        if token.kind == T.PLUS
            sgn = +;
            token, state = nextskip(tokens, state)
        elseif token.kind == T.MINUS
            sgn = -;
            token, state = nextskip(tokens, state)
        end
        push!(vec, sgn(parse(TY, untokenize(token))))
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

    sgn = +;
    if token.kind == T.MINUS
        sgn = -;
        token, state = nextskip(tokens, state)
    elseif token.kind == T.PLUS
        token, state = nextskip(tokens, state)
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
    token, state = nextskip(tokens, state)
    if token.kind != T.AND
        tokenerror(token, "&")
    else
        token, state = nextskip(tokens, state)
    end
    
    while token.kind != T.RSQUARE && token.kind != T.ENDMARKER
        token, state, key = tokensgetkey(token, state, tokens, isEQorRSQUARE)
        if token.kind != T.RSQUARE # Allow [&R] as a valid (empty) dict
            token, state = nextskip(tokens, state)
            if token.kind == T.LBRACE
                token, state, value = parsevector(token, state, tokens)
            else
                sgn = +;
                if token.kind == T.PLUS
                    sgn = +;
                    token, state = nextskip(tokens, state)
                elseif token.kind == T.MINUS
                    sgn = -;
                    token, state = nextskip(tokens, state)
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
            
            token, state = nextskip(tokens, state)
            if token.kind != T.COMMA && token.kind != T.RSQUARE
                tokenerror(token, ", or ]")
            elseif token.kind == T.COMMA
                token, state = nextskip(tokens, state)
            end
        end
    end
    
    if token.kind == T.RSQUARE
        token, state = nextskip(tokens, state)
    end

    return token, state, dict
end

function parsenode(token, state, tokens, tree::TREE,
                   lookup, siblings, istip) where {NL, BL,
                                                   TREE <: AbstractBranchTree{NL, BL}}
    myname = ""
    endkinds = [T.SEMICOLON, T.COLON, T.COMMA, T.RPAREN, T.LSQUARE]
    if token.kind ∉ endkinds # We got a nodename
        token, state, name = tokensgetkey(token, state, tokens,
                                          t -> t.kind ∈ endkinds)
        if isempty(lookup)
            myname = addnode!(tree, name)
        else
            if istip
                if haskey(lookup, name)
                    myname = lookup[name]
                else
                    @warn "Wrongly named tip '$name' found in tree " *
                        "with named tips, adding"
                    myname = addnode!(tree, name)
                end
            else
                if haskey(lookup, name)
                    myname = lookup[name]
                else
                    myname = addnode!(tree, name)
                end
            end
        end
    else
        if istip && !isempty(lookup)
            @warn "Anonymous tip found in tree with named tips"
        end
        myname = addnode!(tree)
    end

    foundcolon = false
    if token.kind == T.COLON
        foundcolon = true
        token, state = nextskip(tokens, state)
    end

    siblings[myname] = Dict{String, Any}()
    if token.kind == T.LSQUARE
        token, state, dict = parsedict(token, state, tokens)
        siblings[myname]["dict"] = dict
        setnoderecord!(tree, myname, dict)
    end

    if token.kind == T.COLON || foundcolon
        if token.kind == T.COLON
            token, state = nextskip(tokens, state)
        end
        sgn = +;
        if token.kind == T.PLUS
            sgn = +;
            token, state = nextskip(tokens, state)
        elseif token.kind == T.MINUS
            sgn = -;
            token, state = nextskip(tokens, state)
        end
        if token.kind ∈ [T.INTEGER, T.FLOAT]
            siblings[myname]["length"] = sgn(parse(Float64, untokenize(token)))
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
    mychildren = Dict{NL, Dict{String, Any}}()
    while token.kind != T.RPAREN && token.kind != T.ENDMARKER
        token, state = nextskip(tokens, state)
        if token.kind == T.LPAREN
            token, state = parsenewick(token, state, tokens, tree,
                                       lookup, depth + 1, mychildren)
        else
            token, state, nodename = parsenode(token, state, tokens, tree,
                                               lookup, mychildren, true)
        end
    end
    token, state = nextskip(tokens, state)
    token, state, nodename = parsenode(token, state, tokens, tree,
                                       lookup, children, false)
    for child in keys(mychildren)
        dict = mychildren[child]
        haskey(dict, "length") ?
            addbranch!(tree, nodename, child, dict["length"]) :
            addbranch!(tree, nodename, child)
    end
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

function parsetaxa(token, state, tokens, taxa)
    if !isDIMENSIONS(token)
        tokenerror(token, "Dimensions")
    end
    token, state = nextskip(tokens, state)
    token, state, ntaxstr = tokensgetkey(token, state, tokens)
    if lowercase(ntaxstr) != "ntax"
        error("Unexpected label '$ntaxstr=' not 'ntax='")
    end
    token, state = nextskip(tokens, state)
    ntax = parse(Int64, untokenize(token))
    token, state = nextskip(tokens, state)
    if token.kind != T.SEMICOLON
        tokenerror(token, ";")
    end
    token, state = nextskip(tokens, state)
    if !isTAXLABELS(token)
        tokenerror(token, "Taxlabels")
    end
    token, state = nextskip(tokens, state)
    while token.kind != T.SEMICOLON && token.kind != T.ENDMARKER
        name = untokenize(token)
        token, state = nextskip(tokens, state)
        taxa[name] = name
    end
    if length(taxa) != ntax
        @warn "Taxa list length ($(length(taxa))) and ntax ($ntax) do not match"
    end

    token, state = nextskip(tokens, state)
    if !checktosemi(isEND, token, state, tokens)
        tokenerror(token, "End;")
    end
    return nextskip(tokens, state)
end

function parsetrees(token, state, tokens,
                    ::Type{TREE}, taxa) where TREE <: AbstractBranchTree{String, Int}
    notaxa = isempty(taxa)
    if isTRANSLATE(token)
        token, state = nextskip(tokens, state)
        while token.kind != T.SEMICOLON && token.kind != T.ENDMARKER
            short = untokenize(token)
            token, state = nextskip(tokens, state)
            proper = untokenize(token)
            token, state = nextskip(tokens, state)
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
                token, state = nextskip(tokens, state)
            end
        end
    end

    token, state = nextskip(tokens, state)
    trees = Dict{String, TREE}()
    treedata = Dict{String, Dict{String, Any}}()
    while isTREE(token)
        token, state = nextskip(tokens, state)
        token, state, treename = tokensgetkey(token, state, tokens,
                                              t -> t.kind ∈ [T.LSQUARE, T.EQ])
        @info "Created a tree called '$treename'"
        trees[treename] = TREE()
        addnodes!(trees[treename], collect(values(taxa)))
        if token.kind == T.LSQUARE
            token, state, treedata[treename] = parsedict(token, state, tokens)
        end
        if !isEQ(token)
            tokenerror(token, "=")
        else
            token, state = nextskip(tokens, state)
            if token.kind == T.LSQUARE
                token, state, _ = parsedict(token, state, tokens)
            end
        end
        if token.kind != T.LPAREN
            tokenerror(token, "(")
        else
            token, state = parsenewick(token, state, tokens, trees[treename], taxa)
        end
    end
    if !checktosemi(isEND, token, state, tokens)
        tokenerror(token, "End;")
    end

    token, state = nextskip(tokens, state)
    return token, state, trees, treedata
end

function parsenexus(token, state, tokens,
                    ::Type{TREE}) where {NL, BL,
                                         TREE <: AbstractBranchTree{NL, BL}}
    trees = missing
    treedata = missing
    taxa = Dict{NL, NL}()
    while isBEGIN(token)
        token, state = nextskip(tokens, state)
        if checktosemi(isTAXA, token, state, tokens)
            token, state = nextskip(tokens, state)
            token, state = parsetaxa(token, state, tokens, taxa)
        elseif checktosemi(isTREES, token, state, tokens)
            token, state = nextskip(tokens, state)
            token, state, trees, treedata = parsetrees(token, state, tokens, TREE, taxa)
        else
            @warn "Unexpected nexus block '$(untokenize(token))', skipping..."
            token, state = nextskip(tokens, state)
            while !checktosemi(isEND, token, state, tokens) && token.kind != T.ENDMARKER
                token, state = nextskip(tokens, state)
            end
        end
    end
    
    if token.kind != T.ENDMARKER
        tokenerror(token, "end of file")
    end
    return trees, treedata
end

function parsenexus(tokens::Tokenize.Lexers.Lexer,
                    ::Type{TREE}) where TREE <: AbstractBranchTree{String, Int}
    token, state = nextskip(tokens, start(tokens))
    if token.kind == T.COMMENT && lowercase(untokenize(token)) == "#nexus"
        token, state = nextskip(tokens, token)
        return parsenexus(token, state, tokens, TREE)
    else
        error("Unexpected $(token.kind) token '$(untokenize(token))' " *
              "at start of nexus file")
    end
end

parsenexus(io::IOBuffer, ::Type{TREE}) where TREE <: AbstractBranchTree =
    parsenexus(tokenize(io), TREE)

function parsenexus(ios::IOStream, ::Type{TREE}) where TREE <: AbstractBranchTree
    buf = IOBuffer()
    print(buf, read(ios, String))
    seek(buf, 0)
    return parsenexus(buf, TREE)
end

parsenexus(inp) = parsenexus(inp, NamedTree)
