
using ImportMacros
@import LightXML as XML
using StructArrays

function parseCA!(rr2::RR2, constrArray, node)
    for constr in node
        attr = XML.attributes_dict(constr)
        ub = parse(UInt8, attr["max"])
        pen = attr["type"] == "SOFT" ? parse(UInt8, attr["penalty"]) : 0
        modestr = XML.name(constr) == "CA1" ? "mode" : "mode1"
        if attr[modestr] == "H"
            mode = MODE_H
        elseif attr[modestr] == "A"
            mode = MODE_A
        else
            mode = MODE_HA
        end

        if XML.name(constr) == "CA1"
            slots = IdxSet(parse.(Int64, split(attr["slots"], ";")).+1)
            teams1 = IdxSet(2^parse(Int64, attr["teams"]))
            teams2 = IdxSet(2^rr2.nteams - 1)
            push!(constrArray, CAcon(teams1, teams2, slots, ub, mode, pen))
        elseif XML.name(constr) == "CA2"
            slots = IdxSet(parse.(Int64, split(attr["slots"], ";")).+1)
            teams1 = IdxSet(parse.(Int64, split(attr["teams1"], ";")).+1)
            teams2 = IdxSet(parse.(Int64, split(attr["teams2"], ";")).+1)
            push!(constrArray, CAcon(teams1, teams2, slots, ub, mode, pen))
        elseif XML.name(constr) == "CA3"
            teams1 = IdxSet(parse.(Int64, split(attr["teams1"], ";")).+1)
            teams2 = IdxSet(parse.(Int64, split(attr["teams2"], ";")).+1)
            intp = parse(Int64, attr["intp"])
            binrepr = 2^intp-1
            for start in 1:rr2.nslots - intp
                push!(constrArray, CAcon(teams1, teams2, IdxSet(binrepr), ub, mode, pen))
                binrepr <<= 1
            end
        elseif XML.name(constr) == "CA4"
            teams1 = IdxSet(parse.(Int64, split(attr["teams1"], ";")).+1)
            teams2 = IdxSet(parse.(Int64, split(attr["teams2"], ";")).+1)
            if attr["mode2"] == "GLOBAL"
                slots = IdxSet(parse.(Int64, split(attr["slots"], ";")).+1)
                push!(constrArray, CAcon(teams1, teams2, slots, ub, mode, pen))
            else
                for s in parse.(Int64, split(attr["slots"], ";"))
                    slots = IdxSet(2^s)
                    push!(constrArray, CAcon(teams1, teams2, slots, ub, mode, pen))
                end
            end
        else
            throw(ArgumentError("invalid node"))
        end
    end
end

function parseBR!(constrArray, node)
    for constr in node
        attr = XML.attributes_dict(constr)
        slots = IdxSet(parse.(Int64, split(attr["slots"], ";")).+1)
        pen = attr["type"] == "SOFT" ? parse(UInt8, attr["penalty"]) : 0
        ub = parse(UInt8, attr["intp"])
        modestr = XML.name(constr) == "BR1" ? "mode2" : "homeMode"
        if attr[modestr] == "H"
            mode = MODE_H
        elseif attr[modestr] == "A"
            mode = MODE_A
        else
            mode = MODE_HA
        end
        if XML.name(constr) == "BR1"
            for team in parse.(Int64, split(attr["teams"], ";"))
                teams = IdxSet(2^team)
                push!(constrArray, BRcon(teams,slots,ub,mode,pen))
            end
        else
            teams = IdxSet(parse.(Int64, split(attr["teams"], ";")).+1)
            push!(constrArray, BRcon(teams,slots,ub,mode,pen))
        end
    end
end

function parseGA!(constrArray, node)
    for constr in node
        attr = XML.attributes_dict(constr)
        pen = attr["type"] == "SOFT" ? parse(UInt8, attr["penalty"]) : 0
        slots = IdxSet(parse.(Int64, split(attr["slots"], ";")).+1)
        ub = parse(UInt8, attr["max"])
        lb = parse(UInt8, attr["min"])
        splittotuple = z->tuple(parse.(Int64,split(z, ",")).+1...)
        meetings = splittotuple.(split(strip(attr["meetings"], [';']), ";"))
        push!(constrArray, GAcon(slots, lb, ub, meetings, pen))
    end

end

function parseFA!(constrArray, node)
    for constr in node
        attr = XML.attributes_dict(constr)
        teams = IdxSet(parse.(Int64, split(attr["teams"], ";")).+1)
        slots = IdxSet(parse.(Int64, split(attr["slots"], ";")).+1)
        ub = parse(UInt8, attr["intp"])
        pen = attr["type"] == "SOFT" ? parse(UInt8, attr["penalty"]) : 0
        push!(constrArray, FAcon(teams, slots, ub, pen))
    end
end

function parseSE!(constrArray, node)
    for constr in node
        attr = XML.attributes_dict(constr)
        teams = IdxSet(parse.(Int64, split(attr["teams"], ";")).+1)
        pen = attr["type"] == "SOFT" ? parse(UInt8, attr["penalty"]) : 0
        lb = parse(UInt8, attr["min"])
        push!(constrArray, SEcon(teams, lb, pen))
    end
end

function parseRobinX(filename)
    xdoc = XML.parse_file(filename)
    root = XML.root(xdoc)

    formatNode = root["Structure"][1]["Format"][1]
    resourceNode = root["Resources"][1]
    constraintNode = root["Constraints"][1]

    nteams = length(resourceNode["Teams"][1]["team"])
    nslots = 2nteams-2
    isPhased = formatNode["gameMode"] == "P"

    rr2 = RR2(nteams, nslots, isPhased)

    # capacity constraints
    constrCA = Vector{CAcon}(undef, 0)
    parseCA!(rr2, constrCA, constraintNode["CapacityConstraints"][1]["CA1"])
    parseCA!(rr2, constrCA, constraintNode["CapacityConstraints"][1]["CA2"])
    parseCA!(rr2, constrCA, constraintNode["CapacityConstraints"][1]["CA3"])
    parseCA!(rr2, constrCA, constraintNode["CapacityConstraints"][1]["CA4"])

    # break constraints
    constrBR = Vector{BRcon}(undef, 0)
    parseBR!(constrBR, constraintNode["BreakConstraints"][1]["BR1"])
    parseBR!(constrBR, constraintNode["BreakConstraints"][1]["BR2"])

    # game constraints
    constrGA = Vector{GAcon}(undef, 0)
    parseGA!(constrGA, constraintNode["GameConstraints"][1]["GA1"])

    # fairness constraints
    constrFA = Vector{FAcon}(undef, 0)
    parseFA!(constrFA, constraintNode["FairnessConstraints"][1]["FA2"])

    # separation constraints
    constrSE = Vector{SEcon}(undef, 0)
    parseSE!(constrSE, constraintNode["SeparationConstraints"][1]["SE1"])

    return nteams, isPhased, constrCA, constrBR, constrGA, constrFA, constrSE
end
