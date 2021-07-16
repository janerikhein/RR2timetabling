using JuMP, Gurobi
using Combinatorics
#=
variables:
    x[i,j,s] - Team i plays team j at home in slot s.
    bh[i,s] - Team i has a home break in slot s
    ba[i,s] - Team i has an away break in slot s

constraints:.
    hgames[s] - Slot s features exactly nteams/2 home games.
    onefactor[i,s] - Team i participates exacty once in slot s.
    alldiffh[i,j] - Team i plays exactly one home game against team j.
    phased[i,j] - Teams i,j meet exactly once in the first part of the
                  tournament
    hasym[i,j,s] - games (i,j) and (j,i) cannot occur in the same slot

annonymous constraints:
    Team i does not play itself in any slot.

=#
function baseIPModel(nteams, isphased)
    nslots = 2*nteams-2
    model = direct_model(Gurobi.Optimizer())

    @variable(model, x[1:nteams, 1:nteams, 1:nslots], Bin)
    @variable(model, bh[1:nteams, 1:nslots], Bin)
    @variable(model, ba[1:nteams, 1:nslots], Bin)


    @constraint(model, [i = 1:nteams], x[i,i, : ] .== 0 )
    @constraint(model, hgames[s = 1:nteams], sum(x[:,:,s]) == div(nteams,2))
    @constraint(model, onefactor[i = 1:nteams, s = 1:nslots],
                sum(x[i,:,s]) + sum(x[:,i,s]) == 1)
    @constraint(model, alldiffh[i= 1:nteams, j=1:nteams ; i!=j],
                sum(x[i,j,:]) == 1)
    @constraint(model, [i = 1:nteams, j=1:nteams, s=1:nslots],
                x[i,j,s] + x[j,i,s] <= 1)
    if isphased
        @constraint(model, phased[i = 1:nteams, j = 1:nteams ; i!=j],
            sum(x[i,j,1:div(nslots,2)]) + sum(x[j,i,1:div(nslots,2)]) == 1)
    end

    @constraint(model, [i = 1:nteams, s = 2:nslots],
                ba[i,s] - bh[i,s] + sum(x[i,:,s]) + sum(x[i,:,s-1]) == 1)

    @constraint(model, [i = 1:nteams, s = 2:nslots], bh[i,s] + ba[i,s] <= 1)

    model
end

function addOffsetVar!(model, con)
    offset = @variable(model, lower_bound=0)
    add_to_expression!(model[:obj], con.pen, offset)
    offset
end

@inline function exprCA(model, ca)
    x = model[:x]
    offset = ca.pen == 0 ? 0 : addOffsetVar!(model, ca)
    if ca.mode == MODE_H
        return @expression(model, sum(x[i,j,s] for i in ca.teams1
                for j in ca.teams2 for s in ca.slots) - offset)
    elseif ca.mode == MODE_A
        return @expression(model, sum(x[j,i,s] for i in ca.teams1
                for j in ca.teams2 for s in ca.slots) - offset)
    elseif ca.mode == MODE_HA
        return @expression(model, sum(x[j,i,s]+x[i,j,s] for i in ca.teams1
                for j in ca.teams2 for s in ca.slots) - offset)
    end
end

@inline function exprBR(model, br)
    bh, ba = model[:bh], model[:ba]
    offset = br.pen == 0 ? 0 : addOffsetVar!(model, br)
    if br.mode == MODE_H
        return @expression(model, sum(bh[i,s] for i in br.teams
                for s in br.slots) - offset)
    elseif br.mode == MODE_A
        return @expression(model, sum(ba[i,s] for i in br.teams
                for s in br.slots) - offset)
    elseif br.mode == MODE_HA
        return @expression(model, sum(bh[i,s] + ba[i,s] for i in br.teams
                for s in br.slots) - offset)
    end
end

# upperbound = true if expr <= ub
@inline function exprGA(model, ga, upperbound)
    x = model[:x]
    offset = ga.pen == 0 ? 0 : addOffsetVar!(model, ga)
    if upperbound
        return @expression(model, sum(x[ga.meetings[i][1],ga.meetings[i][2], s]
            for i =1:length(ga.meetings) for s in ga.slots) - offset)
    else
        return @expression(model, sum(x[ga.meetings[i][1],ga.meetings[i][2], s]
            for i =1:length(ga.meetings) for s in ga.slots) + offset)
    end
end

@inline function exprSE(model, se, i, j, upperbound)
    x = model[:x]
    offset = se.pen == 0 ? 0 : addOffsetVar!(model, se)
    if upperbound
        return @expression(model, sum(s*x[i,j,s] - s*x[j,i,s] for s=1:size(x,3)) - offset)
    else
        return @expression(model, sum(s*x[i,j,s] - s*x[j,i,s] for s=1:size(x,3)) + offset)
    end
end

function buildModel(rr2, caCons, brCons, gaCons, faCons, seCons)
    @assert(length(seCons) <= 1)
    @assert(length(faCons) <= 1)

    model = baseIPModel(rr2.nteams, rr2.isphased)

    # objective function expression
    @expression(model, obj, AffExpr(0))

    length(caCons) == 0 || @constraint(model, ca[i = 1:length(caCons)],
        exprCA(model, caCons[i]) <= caCons[i].ub)
    length(brCons) == 0 || @constraint(model, br[i = 1:length(brCons)],
        exprBR(model, brCons[i]) <= brCons[i].ub)
    if length(gaCons) != 0

        @constraint(model, gaub[i = 1:length(gaCons)],
            exprGA(model, gaCons[i], true) <= gaCons[i].ub)
        @constraint(model, galb[i = 1:length(gaCons)],
            gaCons[i].lb <= exprGA(model, gaCons[i], false))
    end

    if length(faCons) > 0
        x = model[:x]
        cfa = faCons[1]
        teams = collect(i for i in cfa.teams)
        slots = collect(s for s in cfa.slots)
        max = @variable(model, [1:length(teams), 1:length(teams)])
        min = @variable(model, [1:length(teams), 1:length(teams)])
        for s in slots
            @constraint(model, [i=1:length(teams), j=1:length(teams); i!=j],
                min[i,j] <= sum(x[teams[i], :, @view slots[1:s]]) <= max[i,j])
            @constraint(model, [i=1:length(teams), j=1:length(teams); i!=j],
                min[i,j] <= sum(x[teams[j], :, @view slots[1:s]]) <= max[i,j])
        end
        offsetfa = cfa.pen == 0 ? 0 : addOffsetVar!(model, cfa)
        @constraint(model, [i=1:length(teams), j=1:length(teams); i!=j],
                    max[i,j] - min[i,j] - offsetfa <= cfa.ub)
    end

    if length(seCons) > 0
        cse = seCons[1]
        for i in cse.teams
            for j in cse.teams
                j<=i && continue
                homefirst = @variable(model, binary=true)
                @constraint(model, sehf[i=1:length(cse.teams), j=1:length(cse.teams)],
                    exprSE(model, cse, i, j, true) <= -cse.lb - 2*nslots*homefirst + 2*nslots-1)
                @constraint(model, seaf[i=1:length(c.teams), j=1:length(cse.teams)],
                    exprSE(model, cse, i, j, false) >= cse.lb - 2*nslots*homefirst)
            end
        end
    end

    @objective(model, MIN, model[:obj])
    model
end
#=
parses solution data into a schedule
=#
function parseSolution(sol)
    sched = Array{Int64}(undef, size(sol,1), size(sol,3))
    for s = 1:size(sched,2)
        for i = 1:size(sched,1)
            for j = 1:size(sched,1)
                if sol[i,j,s] == 1.0
                    sched[i,s] = j
                    sched[j,s] = -i
                end
            end
        end
    end
    sched
end
