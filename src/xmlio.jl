
"""
    parse_robinx(filename::String)

Parses an instance given in the RobinX format.

# Arguments
- `filename::String` : The project path to the file to be parsed

# Returns
- `Instance` : the problem instance
"""
function parse_robinx(filename::String)
    xdoc = XML.parse_file(filename)
    root = XML.root(xdoc)

    formatNode = root["Structure"][1]["Format"][1]
    resourceNode = root["Resources"][1]
    constraintNode = root["Constraints"][1]

    nteams = length(resourceNode["Teams"][1]["team"])
    nslots = 2nteams - 2
    isPhased = XML.content(only(formatNode["gameMode"])) == "P"

    rr2 = RR2(nteams, nslots, isPhased)

    # capacity constraints
    constrCA = Vector{CAcon}(undef, 0)
    parse_ca!(rr2, constrCA, constraintNode["CapacityConstraints"][1]["CA1"])
    parse_ca!(rr2, constrCA, constraintNode["CapacityConstraints"][1]["CA2"])
    parse_ca!(rr2, constrCA, constraintNode["CapacityConstraints"][1]["CA3"])
    parse_ca!(rr2, constrCA, constraintNode["CapacityConstraints"][1]["CA4"])

    # break constraints
    constrBR = Vector{BRcon}(undef, 0)
    parse_br!(constrBR, constraintNode["BreakConstraints"][1]["BR1"])
    parse_br!(constrBR, constraintNode["BreakConstraints"][1]["BR2"])

    # game constraints
    constrGA = Vector{GAcon}(undef, 0)
    parse_ga!(constrGA, constraintNode["GameConstraints"][1]["GA1"])

    # fairness constraints
    constrFA = Vector{FAcon}(undef, 0)
    parse_fa!(constrFA, constraintNode["FairnessConstraints"][1]["FA2"])

    # separation constraints
    constrSE = Vector{SEcon}(undef, 0)
    parse_se!(constrSE, constraintNode["SeparationConstraints"][1]["SE1"])

    rr2 = RR2(nteams, 2*nteams-2, isPhased)

    return Instance(rr2, constrCA, constrBR, constrGA, constrFA, constrSE)
end

"""
    parse_ca!(rr2::RR2, constrArray::Vector{CAcon}, node)

Parses all capacity constraints given by `node` and adds them to `constrArray`.
"""
function parse_ca!(rr2::RR2, constrArray::Vector{CAcon}, node)
    for constr in node
        attr = XML.attributes_dict(constr)
        ub = parse(Int64, attr["max"])
        pen = attr["type"] == "SOFT" ? parse(Int64, attr["penalty"]) : 0
        modestr = XML.name(constr) == "CA1" ? "mode" : "mode1"
        if attr[modestr] == "H"
            mode = MODE_H
        elseif attr[modestr] == "A"
            mode = MODE_A
        else
            mode = MODE_HA
        end

        if XML.name(constr) == "CA1"
            slots = IdxSet(parse.(Int64, split(attr["slots"], ";")) .+ 1)
            teams1 = IdxSet(2^parse(Int64, attr["teams"]))
            teams2 = IdxSet(2^rr2.nteams - 1)
            push!(constrArray, CAcon(teams1, teams2, slots, ub, mode, pen))
        elseif XML.name(constr) == "CA2"
            slots = IdxSet(parse.(Int64, split(attr["slots"], ";")) .+ 1)
            teams1 = IdxSet(parse.(Int64, split(attr["teams1"], ";")) .+ 1)
            teams2 = IdxSet(parse.(Int64, split(attr["teams2"], ";")) .+ 1)
            push!(constrArray, CAcon(teams1, teams2, slots, ub, mode, pen))
        elseif XML.name(constr) == "CA3"
            teams1 = parse.(Int64, split(attr["teams1"], ";"))
            teams2 = IdxSet(parse.(Int64, split(attr["teams2"], ";")) .+ 1)
            intp = parse(Int64, attr["intp"])
            binrepr = 2^intp - 1
            for team in teams1
                for start = 0:rr2.nslots-intp
                    push!(
                        constrArray,
                        CAcon(
                            IdxSet(2^team), # teams1
                            teams2,
                            IdxSet(binrepr << start), # slots
                            ub,
                            mode,
                            pen,
                        ),
                    )
                end
            end
        elseif XML.name(constr) == "CA4"
            teams1 = IdxSet(parse.(Int64, split(attr["teams1"], ";")) .+ 1)
            teams2 = IdxSet(parse.(Int64, split(attr["teams2"], ";")) .+ 1)
            if attr["mode2"] == "GLOBAL"
                slots = IdxSet(parse.(Int64, split(attr["slots"], ";")) .+ 1)
                push!(constrArray, CAcon(teams1, teams2, slots, ub, mode, pen))
            else
                for s in parse.(Int64, split(attr["slots"], ";"))
                    slots = IdxSet(2^s)
                    push!(
                        constrArray,
                        CAcon(teams1, teams2, slots, ub, mode, pen),
                    )
                end
            end
        else
            throw(ArgumentError("invalid node"))
        end
    end
end

"""
    parse_br!(rr2::RR2, constrArray::Vector{BRcon}, node)

Parses all break constraints given by `node` and adds them to `constrArray`.
"""
function parse_br!(constrArray, node)
    for constr in node
        attr = XML.attributes_dict(constr)
        slots = IdxSet(parse.(Int64, split(attr["slots"], ";")) .+ 1)
        pen = attr["type"] == "SOFT" ? parse(Int64, attr["penalty"]) : 0
        ub = parse(Int64, attr["intp"])
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
                push!(constrArray, BRcon(teams, slots, ub, mode, pen))
            end
        else
            teams = IdxSet(parse.(Int64, split(attr["teams"], ";")) .+ 1)
            push!(constrArray, BRcon(teams, slots, ub, mode, pen))
        end
    end
end

"""
    parse_ga!(rr2::RR2, constrArray::Vector{GAcon}, node)

Parses all game constraints given by `node` and adds them to `constrArray`.
"""
function parse_ga!(constrArray, node)
    for constr in node
        attr = XML.attributes_dict(constr)
        pen = attr["type"] == "SOFT" ? parse(Int64, attr["penalty"]) : 0
        slots = IdxSet(parse.(Int64, split(attr["slots"], ";")) .+ 1)
        ub = parse(Int64, attr["max"])
        lb = parse(Int64, attr["min"])
        splittotuple = z -> tuple(parse.(Int64, split(z, ",")) .+ 1...)
        meetings = splittotuple.(split(strip(attr["meetings"], [';']), ";"))
        push!(constrArray, GAcon(slots, lb, ub, meetings, pen))
    end
end

"""
    parse_fa!(rr2::RR2, constrArray::Vector{FAcon}, node)

Parses all fairness constraints given by `node` and adds them to `constrArray`.
"""
function parse_fa!(constrArray, node)
    for constr in node
        attr = XML.attributes_dict(constr)
        teams = IdxSet(parse.(Int64, split(attr["teams"], ";")) .+ 1)
        slots = IdxSet(parse.(Int64, split(attr["slots"], ";")) .+ 1)
        ub = parse(Int64, attr["intp"])
        pen = attr["type"] == "SOFT" ? parse(Int64, attr["penalty"]) : 0
        push!(constrArray, FAcon(teams, slots, ub, pen))
    end
end

"""
    parse_se!(rr2::RR2, constrArray::Vector{SEcon}, node)

Parses all separation constraints given by `node` and adds them to `constrArray`.
"""
function parse_se!(constrArray, node)
    for constr in node
        attr = XML.attributes_dict(constr)
        teams = IdxSet(parse.(Int64, split(attr["teams"], ";")) .+ 1)
        pen = attr["type"] == "SOFT" ? parse(Int64, attr["penalty"]) : 0
        lb = parse(Int64, attr["min"])
        push!(constrArray, SEcon(teams, lb, pen))
    end
end

"""
    toxml(sched::Matrix{Int}, filename::String)

Writes a xml representation of the schedule `sched` to a file specified by
`filename`.
"""
function toxml(sched::Matrix{Int}, filename::String)
    open(filename, "w") do io
        write(
            io,
            "<Solution>
            <MetaData>
            <SolutionName>mySolution.xml</SolutionName>
            <InstanceName>myInstance.xml</InstanceName>
            <ObjectiveValue infeasibility=\"0\" objective=\"0\"/>
            </MetaData>
            <Games>",
        )
        for i = 1:size(sched, 1)
            for s = 1:size(sched, 2)
                if sign(sched[i, s]) == 1
                    write(
                        io,
                        "<ScheduledMatch home=\"$(i-1)\" away=\"$(sched[i,s]-1)\" slot=\"$(s-1)\"/>\n",
                    )
                end
            end
        end
        write(io, "</Games>\n</Solution>")
    end
end

function tofile(sched::Matrix{Int}, filename::String)
    open(filename, "w") do io
        for i = 1:size(sched, 1)
            write(io, join(sched[i,:], " "))
            write(io, "\n")
        end
    end
end
