#=
# Tree → Subtree ";"
# Subtree → Leaf | Internal
# Leaf → Name
# Internal → "(" BranchSet ")" Name_or_Support
# BranchSet → Branch | Branch "," BranchSet
# Branch → Subtree Length
# Name_or_Support → empty | string | number
# Length → empty | ":" FPNum
# FPNum number

# open_meta + ("[^"]*"+|[^,=\s]+) + spc + (=\s*(\{[^=}]*\}|"[^"]*"+|[^,]+))? + close_meta

using ParserCombinator

spc = Drop(Star(Space()))
blank = E""
quotes = E"\""
open_meta = E"[&" # |> "__start_meta"
close_meta = E"]" # |> "__end_meta"
open_branchset = E"(" # |> "__start_subtree"
close_branchset = E")" # |> "__end_subtree"
next_one = E"," # |> "__next_one"
equals = E"="

@with_names begin
@with_pre spc begin
    str =  p"[a-zA-Z_][a-zA-Z\d_]*" | (quotes + p"[\"]+" + quotes)
    name = str | blank
    length_sep = E":" # |> "__length"
    length_val = Parse(p"-?(\d*\.?\d+|\d+\.\d*)([eE]-?\d+)?", Float64)
    datum = (p"\"[^\"]*\"" | p"[^,=\s\"]+") + spc + equals + spc + (p"{[^=}]*}" | p"\"[^\"]*\"" | p"[^,]+")
    data = Delayed()
    data.matcher = datum + ((next_one + data) | blank)
    metacomments = (open_meta + data + close_meta) | blank
    node_metacomments = metacomments
    branch_metacomments = metacomments
    blength = (length_sep + spc + branch_metacomments + length_val) | blank
    leaf = (str | blank) + node_metacomments # |> "__leaf"
    subtree = Delayed()
    branch = subtree + node_metacomments + blength
    branchset = Delayed()
    branchset.matcher = branch + ((next_one + branchset) | blank)
    fpnum = Parse(p"(\d*\.?\d+|\d+\.\d*)", Float64)
    internal = open_branchset + branchset + close_branchset + (fpnum | name) + node_metacomments # |> "__internal"
    subtree.matcher = leaf | internal
    endnewick = E";" + spc + Eos() # |> "__end_tree"
    newick = subtree + endnewick
end
end

@time open("test/Qian2016.tree", "r") do io
    parse_try(io, newick)
end
io = open("test/Qian2016.tree", "r")
out = read(io, String)
parse_dbg(out, Trace(newick))
@time parse_one(out, newick)

tree = "(A9CIT6:0.59280764,P51981:0.55926221,(Q8A861:0.99105703,((Q81IL5:0.76431643,((((A1B198:0.94287572,(Q2U1E8:0.71410953,Q5LT20:0.55480466)0.975579:0.21734008)1.000000:0.55075633,(Q92YR6:1.11236546,(((Q13PB7:1.33807955,(Q161M1:1.17720944,A1AYL4:0.93440931)0.784619:0.18922325)0.878426:0.09769089,((((Q1ASJ4:1.28537785,(A8LS88:0.87558406,(C1DMY1:0.14671933,(P77215:0.02112667,Q8ZNF9:0.01593493)1.000000:0.35900384)1.000000:0.81055398)0.999947:0.27041496)0.403932:0.04809748,(A4YVM8:1.35286455,(Q9RKF7:0.83804265,((Q8ZL58:0.21550115,Q12GE3:0.23170031)1.000000:0.65091551,(Q7CU39:0.71681321,(Q8P3K2:0.27030998,(Q1GLV3:0.34268834,Q7L5Y1:0.42965239)0.843769:0.09133540)1.000000:1.04593860)0.978319:0.19792131)0.995516:0.18393058)0.684889:0.10148106)0.965570:0.12638685)0.999013:0.10597436,(((A0A0H3LM82:1.53892471,(O06741:1.56982104,(G0L7B8:0.68911617,A9CEQ8:0.63642012)0.999148:0.26000097)0.760390:0.07007390)0.534387:0.04760860,(A0A0H3LT39:0.95322505,(Q3HKK5:1.65509086,(A8H7M5:1.21743086,A8H9D1:0.47372214)0.711619:0.14571049)0.992889:0.20930969)0.824619:0.10359095)0.957749:0.09170993,(((Q8ZNH1:0.35062372,Q5NN22:0.45517287)1.000000:0.51175191,(Q7D1T6:0.17783006,Q63IJ7:0.15483880)1.000000:0.87953156)0.879253:0.11535090,(Q28RT0:1.01784576,(B9JNP7:1.11261669,(B2UCA8:0.76348582,(((A6M2W4:0.21565444,A4W7D6:0.17558479)1.000000:0.37879482,(D4GJ14:0.06604958,(C6CBG9:0.01850349,B5R541:0.03323447)0.999985:0.05738427)1.000000:0.29155266)1.000000:0.31191076,((C6D9S0:0.03964480,Q8FHC7:0.03209453)1.000000:0.13486764,((A4XF23:0.22295702,(Q1QT89:0.18122799,B3PDB1:0.14414146)0.999998:0.06015729)0.968272:0.04285536,(B0T0B1:0.05224760,Q9AAR4:0.07234240)1.000000:0.14812779)0.999678:0.08813897)1.000000:0.36147407)1.000000:0.47112287)0.928316:0.13341811)0.550875:0.06620266)0.961209:0.13640648)1.000000:0.27744895)0.988757:0.09058974)0.866782:0.06631759,(A9CL63:0.35266112,Q92ZS5:0.19783599)1.000000:1.18369094)0.402267:0.02752454)0.999888:0.19319607,(Q9F3A5:0.58548261,((Q3KB33:0.33553686,(A0R5B5:0.21484600,D6Y7Y6:0.27848012)1.000000:0.24287080)0.942008:0.15046146,(Q1QUN0:0.25266568,(C0WBB5:0.25833170,(P0AES2:0.09059530,A6VQF1:0.05673552)1.000000:0.13075826)0.999738:0.14696003)1.000000:0.81155149)1.000000:0.47139015)1.000000:0.50975221)0.966695:0.13031427)0.863808:0.08043994)0.973288:0.08656638,(Q5P025:0.97873157,(C5CFI0:1.08126948,((A0A0H3KH80:0.00000001,(A0A0H2WWB5:0.00000001,Q53635:0.00000001)-1.000000:0.00000001)-1.000000:1.77379709,((Q838J7:0.42816989,Q927X3:0.43775742)1.000000:0.24306201,(Q5SJX8:0.32151423,Q9RYA6:0.31040549)1.000000:0.38120655)0.761411:0.11121216)0.449066:0.07920285)1.000000:0.65145569)0.934431:0.11633312)0.809578:0.05916034,(A0QTN8:1.01771865,((Q8DJP8:2.14313928,(Q8NN12:0.75309341,(P05404:0.60462842,(Q4K9X1:0.32854393,A6T9N5:0.34558716)1.000000:0.22924230)0.999793:0.18933082)0.610131:0.11821100)0.811426:0.14174769,(A8HTB8:0.54448819,(Q5LM96:0.23356964,Q28SI7:0.15669417)1.000000:0.47534595)1.000000:0.49079583)0.999955:0.20013288)0.717358:0.06863833)0.729022:0.08518511)0.998543:0.11833475,(Q607C7:1.08575204,(((O34508:0.48939847,B0TZW0:1.11933860)0.879395:0.12702564,(Q9WXM1:0.80769788,(A9B055:0.36657533,(A5UXJ3:0.34547016,A9GEI3:0.26128229)0.974762:0.12479711)1.000000:0.52751375)1.000000:0.33586771)0.610795:0.07157613,(Q11T61:0.92177095,Q834W6:0.44934225)0.999732:0.14429535)0.996719:0.09875344)0.787853:0.04813874)0.983918:0.24536033)1.000000:0.66916649);"
tree = "(A9CIT6[&a=2]:0.59280764,P51981:0.55926221,(Q8A861:0.99105703,((Q81IL5:0.76431643,((((A1B198:0.94287572,(Q2U1E8:0.71410953,Q5LT20:0.55480466)0.975579:0.21734008)1.000000:0.55075633,(Q92YR6:1.11236546,(((Q13PB7:1.33807955,(Q161M1:1.17720944,A1AYL4:0.93440931)0.784619:0.18922325)0.878426:0.09769089,((((Q1ASJ4:1.28537785,(A8LS88:0.87558406,(C1DMY1:0.14671933,(P77215:0.02112667,Q8ZNF9:0.01593493)1.000000:0.35900384)1.000000:0.81055398)0.999947:0.27041496)0.403932:0.04809748)))))))))))));"
parse_dbg(tree, Trace(newick))

=#

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
    result = isnothing(state) ? iterate(tokens) : iterate(tokens, state)
    isnothing(result) && return nothing
    token, state = result
    while isWHITESPACE(token)
        result = iterate(tokens, state)
        isnothing(result) && return nothing
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
        isnothing(result) && return nothing
        token, state = result
    end
    return token, state, join(sofar)
end

function checktosemi(test::Function, token, state, tokens)
    if !test(token)
        return false
    end
    result = iterateskip(tokens, state)
    isnothing(result) && return nothing
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
            isnothing(result) && return nothing
            token, state = result
        elseif token.kind == T.MINUS
            sgn = -;
            result = iterateskip(tokens, state)
            isnothing(result) && return nothing
            token, state = result
        end
        push!(vec, sgn(parse(TY, untokenize(token))))
        result = iterateskip(tokens, state)
        isnothing(result) && return nothing
        token, state = result
        if token.kind != T.COMMA && token.kind != T.RBRACE
            tokenerror(token, ",")
        elseif token.kind == T.COMMA
            result = iterateskip(tokens, state)
            isnothing(result) && return nothing
            token, state = result
        end
    end

    return token, state, vec
end

function parsevector(token, state, tokens)
    result = iterateskip(tokens, state)
    isnothing(result) && return nothing
    token, state = result

    sgn = +;
    if token.kind == T.MINUS
        sgn = -;
        result = iterateskip(tokens, state)
        isnothing(result) && return nothing
        token, state = result
    elseif token.kind == T.PLUS
        result = iterateskip(tokens, state)
        isnothing(result) && return nothing
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
        isnothing(result) && return nothing
        token, state = result
        if token.kind != T.COMMA && token.kind != T.RBRACE
            tokenerror(token, ",")
        elseif token.kind == T.COMMA
            result = iterateskip(tokens, state)
            isnothing(result) && return nothing
            token, state = result
        end
    end

    return token, state, vec
end

function parsedict(token, state, tokens)
    dict = Dict{String, Any}()
    result = iterateskip(tokens, state)
    isnothing(result) && return nothing
    token, state = result
    if token.kind != T.AND
        tokenerror(token, "&")
    else
        result = iterateskip(tokens, state)
        isnothing(result) && return nothing
        token, state = result
    end

    while token.kind != T.RSQUARE && token.kind != T.ENDMARKER
        token, state, key = tokensgetkey(token, state, tokens, isEQorRSQUARE)
        if token.kind != T.RSQUARE # Allow [&R] as a valid (empty) dict
            result = iterateskip(tokens, state)
            isnothing(result) && return nothing
            token, state = result
            if token.kind == T.LBRACE
                token, state, value = parsevector(token, state, tokens)
            else
                sgn = +;
                if token.kind == T.PLUS
                    sgn = +;
                    result = iterateskip(tokens, state)
                    isnothing(result) && return nothing
                    token, state = result
                elseif token.kind == T.MINUS
                    sgn = -;
                    result = iterateskip(tokens, state)
                    isnothing(result) && return nothing
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
            isnothing(result) && return nothing
            token, state = result
            if token.kind != T.COMMA && token.kind != T.RSQUARE
                tokenerror(token, ", or ]")
            elseif token.kind == T.COMMA
                result = iterateskip(tokens, state)
                isnothing(result) && return nothing
                token, state = result
            end
        end
    end

    if token.kind == T.RSQUARE
        result = iterateskip(tokens, state)
        isnothing(result) && return nothing
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
        isnothing(result) && return nothing
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
            isnothing(result) && return nothing
            token, state = result
        end
        sgn = +;
        if token.kind == T.PLUS
            sgn = +;
            result = iterateskip(tokens, state)
            isnothing(result) && return nothing
            token, state = result
        elseif token.kind == T.MINUS
            sgn = -;
            result = iterateskip(tokens, state)
            isnothing(result) && return nothing
            token, state = result
        end
        if token.kind ∈ [T.INTEGER, T.FLOAT]
            siblings[getnodename(tree, myname)]["length"] =
                sgn(parse(Float64, untokenize(token)))
        else
            tokenerror(token, "a length")
        end
        result = iterateskip(tokens, state)
        isnothing(result) && return nothing
        token, state = result
    end

    return token, state, getnodename(tree, myname)
end

function parsenewick!(token, state, tokens, tree::TREE,
                     lookup = Dict(), depth = 0,
                     children = Dict{NL, Dict{String, Any}}()) where
    {RT, NL, N, B, TREE <: AbstractTree{OneTree, RT, NL, N, B}}
    mychildren = Dict{NL, Dict{String, Any}}()
    while (token.kind != T.RPAREN) & (token.kind != T.ENDMARKER)
        result = iterateskip(tokens, state)
        isnothing(result) &&
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
    isnothing(result) &&
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
            isnothing(result) && return nothing
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
    if isnothing(result)
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
