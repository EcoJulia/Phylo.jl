"""
    Nexus

Type to specify nexus format for input or output.
"""
struct Nexus <: NewickLike end

treeOutputType(::Type{<: AbstractTree{ManyTrees}}) = Nexus

function outputtree!(io::IO, tree::AbstractTree{TT, OneRoot}, ::Nexus) where TT
    println(io, "#NEXUS")
    println(io)
    leaves = getleafnames(tree)
    d = nothing
    if !any(isnothing.(tryparse.(Int, leaves)))
        leaves = sort(parse.(Int, leaves))
    else
        println(io, "BEGIN TAXA;")
        println(io, "    DIMENSIONS NTAX=$(nleaves(tree));")
        println(io, "    TAXLABELS")
        for leaf in leaves
            println(io, "        $(replace(leaf, r"[ \t\n]" => "_"))")
        end
        println(io, "    ;")
        println(io, "END;")

        d = Dict{eltype(leaves), Int}()
        println(io)
    end
    println(io, "BEGIN TREES;")
    if eltype(leaves) ≠ Int
        println(io, "    TRANSLATE")
        for (i, leaf) in enumerate(leaves)
            print(io, "        $i $(replace(leaf, r"[ \t\n]" => "_"))")
            d[leaf] = i
            if i < nleaves(tree)
                println(io, ",")
            else
                println(io, ";")
            end
        end
    end

    for tn in gettreenames(tree)
        println(io)
        print(io, "TREE $tn")
        t1 = gettree(tree, tn)
        outputmetacomment!(io, gettreeinfo(tree)[tn], Nexus())
        print(io, " = [&R] ")
        outputtree!(io, t1, getroot(t1), Newick(d))
        println(io)
    end
    println(io, "END;")
end

function parsetaxa(token, state, tokens, taxa)
    if !isDIMENSIONS(token)
        tokenerror(token, "Dimensions")
    end
    result = iterateskip(tokens, state)
    isnothing(result) && return nothing
    token, state = result
    token, state, ntaxstr = tokensgetkey(token, state, tokens)
    if lowercase(ntaxstr) != "ntax"
        error("Unexpected label '$ntaxstr=' not 'ntax='")
    end
    result = iterateskip(tokens, state)
    isnothing(result) && return nothing
    token, state = result
    ntax = parse(Int64, untokenize(token))
    result = iterateskip(tokens, state)
    isnothing(result) && return nothing
    token, state = result
    if token.kind != T.SEMICOLON
        tokenerror(token, ";")
    end
    result = iterateskip(tokens, state)
    isnothing(result) && return nothing
    token, state = result
    if !isTAXLABELS(token)
        tokenerror(token, "Taxlabels")
    end
    result = iterateskip(tokens, state)
    isnothing(result) && return nothing
    token, state = result
    while token.kind != T.SEMICOLON && token.kind != T.ENDMARKER
        name = untokenize(token)
        if token.kind ∈ [T.STRING, T.CHAR]
            name = name[2:end-1]
        end
        result = iterateskip(tokens, state)
        isnothing(result) && return nothing
        token, state = result
        taxa[name] = name
        if token.kind == T.LSQUARE
            while token.kind != T.RSQUARE
                result = iterateskip(tokens, state)
                isnothing(result) && return nothing
                token, state = result
            end
            result = iterateskip(tokens, state)
            isnothing(result) && return nothing
            token, state = result
        end
    end
    if length(taxa) != ntax
        @warn "Taxa list length ($(length(taxa))) and ntax ($ntax) do not match"
    end

    result = iterateskip(tokens, state)
    isnothing(result) && return nothing
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
        isnothing(result) && return nothing
        token, state = result
        while token.kind != T.SEMICOLON && token.kind != T.ENDMARKER
            short = untokenize(token)
            result = iterateskip(tokens, state)
            isnothing(result) && return nothing
            token, state = result
            proper = untokenize(token)
            result = iterateskip(tokens, state)
            isnothing(result) && return nothing
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
                isnothing(result) && return nothing
                token, state = result
            end
        end
        result = iterateskip(tokens, state)
        isnothing(result) && return nothing
        token, state = result
    end

    trees = Dict{String, TREE}()
    treedata = Dict{String, Dict{String, Any}}()
    numTrees = 0
    while isTREE(token)
        result = iterateskip(tokens, state)
        isnothing(result) && return nothing
        token, state = result
        token, state, treename = tokensgetkey(token, state, tokens,
                                              t -> t.kind ∈ [T.LSQUARE, T.EQ])
        numTrees += 1
        if numTrees < 10
            @info "Created a tree called '$treename'"
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
            isnothing(result) && return nothing
            token, state = result
            if token.kind == T.LSQUARE
                token, state, _ = parsedict(token, state, tokens)
            end
        end
        if token.kind != T.LPAREN
            tokenerror(token, "(")
        else
            result = parsenewick!(token, state, tokens, trees[treename], taxa)
            isnothing(result) && return nothing
            token, state = result
        end
    end
    if !checktosemi(isEND, token, state, tokens)
        tokenerror(token, "End;")
    end

    result = iterateskip(tokens, state)
    isnothing(result) && return nothing
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
                isnothing(result) && return nothing
                token, state = result
            end
            result = iterateskip(tokens, state)
            isnothing(result) && return nothing
            token, state = result
        else
            result = iterateskip(tokens, state)
            isnothing(result) && return nothing
            token, state = result
            if checktosemi(isTAXA, token, state, tokens)
                result = iterateskip(tokens, state)
                isnothing(result) && return nothing
                token, state = result
                token, state = parsetaxa(token, state, tokens, taxa)
            elseif checktosemi(isTREES, token, state, tokens)
                result = iterateskip(tokens, state)
                isnothing(result) && return nothing
                token, state = result
                token, state, trees, treedata = parsetrees(token, state, tokens, TREE, taxa)
            else
                @warn "Unexpected nexus block '$(untokenize(token))', skipping..."
                result = iterateskip(tokens, state)
                isnothing(result) && return nothing
                token, state = result
                while !checktosemi(isEND, token, state, tokens) && token.kind != T.ENDMARKER
                    result = iterateskip(tokens, state)
                    isnothing(result) && return nothing
                    token, state = result
                end
                result = iterateskip(tokens, state)
                isnothing(result) && return nothing
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
    isnothing(result) && return nothing
    token, state = result
    if token.kind == T.COMMENT && lowercase(untokenize(token)) == "#nexus"
        result = iterateskip(tokens, state)
        isnothing(result) && return nothing
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

parsenewicklike(::Type{T}, str, format::Nexus) where T <: AbstractTree =
    parsenexus(str, eltype(T))
