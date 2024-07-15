# SPDX-License-Identifier: BSD-2-Clause

using Automa

@enum NXToken begin
    nxerror = 0
    nxstartsubtree
    nxendsubtree
    nxcomma
    nxnatural
    nxnumber
    nxwhitespace
    nxtoken
    nxquotedstring
    nxcolon
    nxsemicolon
    nxstartmeta
    nxendmeta
    nxequals
    nxstartmetavector
    nxendmetavector
    nxbegin
    nxend
    nxtaxa
    nxtaxlabels
    nxtrees
    nxtranslate
    nxtree
    nxrootedtree
    nxunrootedtree
end

const NXTOKEN_LOOKUP = Dict{NXToken,
                            RE}(nxstartsubtree => re"\(",
                                nxendsubtree => re"\)",
                                nxcomma => re",",
                                nxnatural => re"[0-9]+",
                                nxnumber => re"-?(([0-9]+\.[0-9]*)|(\.?[0-9]+))",
                                nxwhitespace => re"[ \t\r\n]+",
                                nxtoken => re"[a-zA-Z0-9_\.]+",
                                nxquotedstring => re"'(('')|([\(\)\[\]{}/\\,;:=\*\"`+-<> \t\r\na-zA-Z0-9_\.]+))+'",
                                nxcolon => re":",
                                nxsemicolon => re";",
                                nxstartmeta => re"\[&",
                                nxendmeta => re"\]",
                                nxequals => re"=",
                                nxstartmetavector => re"{",
                                nxendmetavector => re"}",
                                nxbegin => re"[bB][eE][gG][iI][nN]",
                                nxend => re"[eE][nN][dD]",
                                nxtaxa => re"[tT][aA][xX][aA]",
                                nxtaxlabels => re"[tT][aA][xX][lL][aA][bB][eE][lL][sS]",
                                nxtrees => re"[tT][rR][eE][eE][sS]",
                                nxtranslate => re"[tT][rR][aA][nN][sS][lL][aA][tT][eE]",
                                nxtree => re"[tT][rR][eE][eE]",
                                nxrootedtree => re"\[&[rR]\]",
                                nxunrootedtree => re"\[&[uU]\]")

nxquoted2token(str) = replace(str, r"[ \t\n\r]" => s"_")
nxtoken2quoted(str) = "'" * replace(str, r"_" => s" ", r"'" => s"''") * "'"
nxuntoken(str) = replace(str, r"_" => s" ")
nxunquote(str) = replace(str[2:(end - 1)], r"''" => s"'", r"[_\t\n\r]" => s" ")

@enum NewickPos begin
    nwSubTree = 100
    nwBranchSet
    nwBranch
    nwMetaComment
    nwMetaVector
    nwEndTree
end

function generate_treeiterator(tokens...; version = 1)
    make_tokenizer((nxerror,
                    [token => NXTOKEN_LOOKUP[token] for token in tokens]),
                   version = Int(version)) |> eval
    return nothing
end

generate_treeiterator(nxwhitespace, nxtoken, nxquotedstring, nxstartsubtree,
                      version = nwSubTree)
generate_treeiterator(nxwhitespace, nxtoken, nxquotedstring, nxstartmeta,
                      nxcolon, nxcomma, nxendsubtree, version = nwBranchSet)
generate_treeiterator(nxwhitespace, nxnumber, nxstartmeta, nxcomma,
                      nxendsubtree, nxcolon, version = nwBranch)
generate_treeiterator(nxwhitespace, nxtoken, nxquotedstring, nxequals,
                      nxendmeta, nxcomma, nxstartmetavector,
                      version = nwMetaComment)
generate_treeiterator(nxwhitespace, nxtoken, nxquotedstring, nxequals,
                      nxcomma, nxendmetavector, version = nwMetaVector)
generate_treeiterator(nxwhitespace, nxsemicolon, version = nwEndTree)

nxstart() = (1, 0, nxerror)

@enum NexusPos begin
    nxSubTree = 200
    nxBranchSet
    nxBranch
    nxMetaComment
    nxMetaVector
    nxEndTree
end

function parsestring(text, context; location = nxstart(),
                     skipwhitespace = true)
    tk = Automa.tokenize(NXToken, text, Int(context))
    status = iterate(tk, location)
    isnothing(status) && return nothing

    safe = !skipwhitespace || status[1][3] ≠ nxwhitespace
    while !safe
        if status[1][3] == nxerror
            if length(text) - location[1] < 10
                error("Parsing failed, text is '$(text[location[1]:end])'")
            else
                error("Parsing failed, text starts with '" *
                      text[location[1]:(location[1] + 9)] * "[...]'")
            end
            return nothing
        elseif status[1][3] ≠ nxwhitespace
            safe = true
        else
            location = status[2]
            status = iterate(tk, location)
            isnothing(status) && return nothing
        end
    end
    (this, next) = status
    return (from = this[1], to = this[1] + this[2] - 1, token = this[3],
            next = next)
end
