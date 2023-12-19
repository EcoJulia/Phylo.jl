struct Nexus <: NewickLike end

treeOutputType(::Type{<: AbstractTree{ManyTrees}}) = Nexus

function outputtree!(io::IO, tree::AbstractTree{TT, OneRoot}, ::Nexus) where TT
    println(io, "#NEXUS")
    println(io, "BEGIN TAXA;")
    println(io, "    DIMENSIONS NTAX=$(nleaves(tree));")
    println(io, "    TAXLABELS")
    leaves = getleafnames(tree)
    for leaf in leaves
        println(io, "        $(replace(leaf, r"[ \t\n]" => "_"))")
    end
    println(io, "    ;")
    println(io, "END;")

    d = Dict{eltype(leaves), Int}()
    println(io, "BEGIN TREES;")
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

    for tn in gettreenames(tree)
        println(io)
        print(io, "TREE $tn")
        outputmetacomment!(io, gettreeinfo(tree)[tn], Nexus())
        print(io, " = [&R] ")
        outputtree!(io, tree[tn], getroot(tree[tn]), Newick(d))
        println(io)
    end
    println(io, "END;")
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
