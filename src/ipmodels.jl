"""
    base_ip_model(nteams::Int, isphased::Bool)

Builds the base IP model, adding all base constraints for a RR2-timetable with
`nteams` teams and enforces a phased structure if `isphased == true`.

# Arguments

- `nteams::Int`: the number of teams participating
- `isphased::Bool`: sets the format of the tournament to be phased, enforcing
    each pair of teams to meet once in each half of the tournament

# Returns

- `Model` - Gurobi Optimizer model

# Model Description

The model uses the following binary variables:
- `x[i,j,s]` - Team `i` plays team `j` at home in slot `s`.
- `bh[i,s]` - Team `i` has a home break in slot `s`.
- `ba[i,s]` - Team `i` has an away break in slot `s`.
"""
function base_ip_model(nteams::Int, isphased::Bool)
    nslots = 2 * nteams - 2
    model = direct_model(Gurobi.Optimizer(GRB_ENV))

    @variable(model, x[1:nteams, 1:nteams, 1:nslots], Bin)
    @variable(model, bh[1:nteams, 1:nslots], Bin)
    @variable(model, ba[1:nteams, 1:nslots], Bin)


    @constraint(model, [i = 1:nteams], x[i, i, :] .== 0)
    @constraint(model, hgames[s = 1:nteams], sum(x[:, :, s]) == div(nteams, 2))
    @constraint(
        model,
        onefactor[i = 1:nteams, s = 1:nslots],
        sum(x[i, :, s]) + sum(x[:, i, s]) == 1
    )
    @constraint(
        model,
        alldiffh[i = 1:nteams, j = 1:nteams; i != j],
        sum(x[i, j, :]) == 1
    )
    @constraint(
        model,
        [i = 1:nteams, j = 1:nteams, s = 1:nslots],
        x[i, j, s] + x[j, i, s] <= 1
    )
    if isphased
        @constraint(
            model,
            phased[i = 1:nteams, j = 1:nteams; i != j],
            sum(x[i, j, 1:div(nslots, 2)]) + sum(x[j, i, 1:div(nslots, 2)]) ==
            1
        )
    end

    @constraint(
        model,
        [i = 1:nteams, s = 2:nslots],
        ba[i, s] - bh[i, s] + sum(x[i, :, s]) + sum(x[i, :, s-1]) == 1
    )

    @constraint(model, [i = 1:nteams, s = 2:nslots], bh[i, s] + ba[i, s] <= 1)

    return model
end

"""
    ip_model(inst::Instance)

Builds the base IP-model and populates it with the constraints given by `inst`.
"""
function ip_model(inst::Instance)

    @assert(length(inst.constrSE) <= 1)
    @assert(length(inst.constrFA) <= 1)

    rr2 = inst.rr2
    model = base_ip_model(rr2.nteams, rr2.isphased)

    # expression used to build objective function iterativly
    @expression(model, obj, AffExpr(0))

    # adding capacity constraints
    length(inst.constrCA) == 0 || @constraint(
        model,
        ca[i = 1:length(inst.constrCA)],
        expr!(model, inst.constrCA[i]) <= inst.constrCA[i].ub
    )

    # adding break constraints
    length(inst.constrBR) == 0 || @constraint(
        model,
        br[i = 1:length(inst.constrBR)],
        expr!(model, inst.constrBR[i]) <= inst.constrBR[i].ub
    )

    # adding game constraints
    if length(inst.constrGA) != 0
        @constraint(
            model,
            gaub[i = 1:length(inst.constrGA)],
            expr!(model, inst.constrGA[i], true) <= inst.constrGA[i].ub
        )
        @constraint(
            model,
            galb[i = 1:length(inst.constrGA)],
            inst.constrGA[i].lb <= expr!(model, inst.constrGA[i], false)
        )
    end

    # adding fairness constraints
    if length(inst.constrFA) > 0
        cfa = only(inst.constrFA)
        for i in cfa.teams
            for j in cfa.teams
                j <= i && continue
                offset = add_offset_var(model, cfa)
                for s in cfa.slots
                    @constraint(
                        model,
                        expr!(model, cfa, i, j, s) <= cfa.ub + offset
                    )
                    @constraint(
                        model,
                        expr!(model, cfa, j, i, s) <= cfa.ub + offset
                    )
                end
            end
        end
    end

    # adding separation constraints
    if length(inst.constrSE) > 0
        cse = only(inst.constrSE)
        for i in cse.teams
            for j in cse.teams
                j <= i && continue
                homefirst = @variable(model, binary = true)
                @constraint(
                    model,
                    expr!(model, cse, i, j, true) >=
                    cse.lb + 2 * rr2.nslots * homefirst - 2 * rr2.nslots
                )
                @constraint(
                    model,
                    expr!(model, cse, i, j, false) >=
                    cse.lb - 2 * rr2.nslots * homefirst
                )
            end
        end
    end

    @objective(model, Min, model[:obj])
    model
end

"""
    game_alloc_model(
        rr2::RR2,
        patternset::Vector{Int64},
        cons::Vector{Constraint};
        feas_relax = false,
    )

Creates the game allocation model for a given `patternset`. If
`feas_relax == true` a feasibility relaxation is computed, minimizing the offset
of the generic constraints `cons`.
"""
function game_alloc_model(
    rr2::RR2,
    patternset::Vector{Int64},
    cons::Vector{Constraint};
    feas_relax = false,
)
    ps(j, s) = patternset[j] >>> (s - 1) & 1
    gam = direct_model(Gurobi.Optimizer())

    @variable(gam, x[1:rr2.nteams, 1:rr2.nteams, 1:rr2.nslots], Bin)

    # setting pattern set
    for s = 1:rr2.nslots
        for i = 1:rr2.nteams
            for j = 1:rr2.nteams
                if ps(i, s) == 1
                    fix(x[j, i, s], 0; force = true)
                else
                    fix(x[i, j, s], 0; force = true)
                end
            end
        end
    end
    # one-factor
    @constraint(
        gam,
        [i = 1:rr2.nteams, s = 1:rr2.nslots],
        sum(x[i, :, s]) + sum(x[:, i, s]) == 1
    )
    # alldiff
    @constraint(
        gam,
        [i = 1:rr2.nteams, j = 1:rr2.nteams; i != j],
        sum(x[i, j, :]) == 1
    )

    # phased cons
    if rr2.isphased
        @constraint(
            gam,
            [i = 1:rr2.nteams, j = i+1:rr2.nteams],
            sum(x[i, j, 1:rr2.nslots÷2]) + sum(x[j, i, 1:rr2.nslots÷2]) == 1
        )
        @constraint(
            gam,
            [i = 1:rr2.nteams, j = i+1:rr2.nteams],
            sum(x[i, j, rr2.nslots÷2+1:end]) +
            sum(x[j, i, rr2.nslots÷2+1:end]) == 1
        )
    end

    offsets = Vector{VariableRef}()
    for con in cons
        @assert isa(con, CAcon) || isa(con, GAcon)
        if isa(con, CAcon)
            os = feas_relax ? @variable(gam, lower_bound = 0) : 0
            feas_relax && push!(offsets, os)
            @constraint(gam, expr!(gam, con) <= con.ub + os)
        elseif isa(con, GAcon)
            os1 = feas_relax ? @variable(gam, lower_bound = 0) : 0
            os2 = feas_relax ? @variable(gam, lower_bound = 0) : 0
            feas_relax && push!(offsets, os1, os2)
            @constraint(gam, expr!(gam, con, true) <= con.ub + os1)
            @constraint(gam, con.lb <= expr!(gam, con, false) + os2)
        end
    end
    if feas_relax
        @objective(gam, Min, sum(offsets))
    end
    return gam
end

"""
    add_offset_var(model::Model, con::Constraint)

Adds an offset variable for a given soft constraint and updates the objective
expression.
"""
function add_offset_var(model::Model, con::Constraint)
    offset = @variable(model, lower_bound = 0)
    add_to_expression!(model[:obj], con.pen, offset)
    offset
end


"""
    expr!(model::Model, con::Constraint, [args...])

Adds an expression to `model` representing the left-hand-side of the
IP-constraint modelling `con`, or if additional `args` are given, one of the
expresions used to model `con`. If the constraint is `soft` offset variables are
added accordingly.

# Arguments

- `model::Model` - The direct model of the Gurobi Optimizer
- `con::Constraint` - The constraint to be modeled
"""
function expr! end

"""
    expr!(model::Model, ca::CAcon)

Adds a capacity constraint expression for `ca`.
"""
function expr!(model::Model, ca::CAcon)
    x = model[:x]
    offset = ca.pen == 0 ? 0 : add_offset_var(model, ca)
    if ca.mode == MODE_H
        return @expression(
            model,
            sum(
                x[i, j, s] for i in ca.teams1 for j in ca.teams2 for
                s in ca.slots
            ) - offset
        )
    elseif ca.mode == MODE_A
        return @expression(
            model,
            sum(
                x[j, i, s] for i in ca.teams1 for j in ca.teams2 for
                s in ca.slots
            ) - offset
        )
    elseif ca.mode == MODE_HA
        return @expression(
            model,
            sum(
                x[j, i, s] + x[i, j, s] for i in ca.teams1 for j in ca.teams2
                for s in ca.slots
            ) - offset
        )
    end
end

"""
    expr!(model::Model, br::BRcon)

Adds a break constraint expression for `br`.
"""
function expr!(model::Model, br::BRcon)
    bh, ba = model[:bh], model[:ba]
    offset = br.pen == 0 ? 0 : add_offset_var(model, br)
    if br.mode == MODE_H
        return @expression(
            model,
            sum(bh[i, s] for i in br.teams for s in br.slots) - offset
        )
    elseif br.mode == MODE_A
        return @expression(
            model,
            sum(ba[i, s] for i in br.teams for s in br.slots) - offset
        )
    elseif br.mode == MODE_HA
        return @expression(
            model,
            sum(bh[i, s] + ba[i, s] for i in br.teams for s in br.slots) -
            offset
        )
    end
end

"""
    expr!(model::Model, ga::GAcon, upperbound)

Adds upper (`upperbound == true`) or lower (`upperbound == false`) bound
expession for `ga`.
"""
function expr!(model::Model, ga::GAcon, upperbound)
    x = model[:x]
    offset = ga.pen == 0 ? 0 : add_offset_var(model, ga)
    if upperbound
        return @expression(
            model,
            sum(
                x[ga.meetings[i][1], ga.meetings[i][2], s] for
                i = 1:length(ga.meetings) for s in ga.slots
            ) - offset
        )
    else
        return @expression(
            model,
            sum(
                x[ga.meetings[i][1], ga.meetings[i][2], s] for
                i = 1:length(ga.meetings) for s in ga.slots
            ) + offset
        )
    end
end

"""
    expr!(model::Model, se::SEcon, i::Int, j::Int, homefirst::Bool)

Adds expression for the separation of games between teams `i` and `j`.
`homefirst == true` is used to model the case that `i` plays the first game at
home, whereas `homefirst == false` is used to model the case that `i` plays the
first game away.
"""
function expr!(model::Model, se::SEcon, i::Int, j::Int, homefirst::Bool)
    x = model[:x]
    offset = se.pen == 0 ? 0 : add_offset_var(model, se)
    if homefirst
        return @expression(
            model,
            -sum(s * x[i, j, s] - s * x[j, i, s] for s = 1:size(x, 3)) +
            offset - 1
        )
    else
        return @expression(
            model,
            sum(s * x[i, j, s] - s * x[j, i, s] for s = 1:size(x, 3)) + offset -
            1
        )
    end
end

"""
    expr!(model::Model, fa::FAcon, i::Int, j::Int, s::Int)

Adds expression of the difference in home games played after slot `s` between
teams `i` and `j`.
"""
function expr!(model::Model, fa::FAcon, i::Int, j::Int, s::Int)
    x = model[:x]
    return @expression(model, sum(x[i, :, 1:s]) - sum(x[j, :, 1:s]))
end

"""
    parse_solution!(sched::Matrix{Int}, sol, model)

parses solution data `sol` into a schedule `sched`.
"""
function parse_solution!(sched::Matrix{Int}, sol, model)
    for s = 1:size(sched, 2)
        for i = 1:size(sched, 1)
            for j = 1:size(sched, 1)
                if sol[i, j, s] == 1.0
                    sched[i, s] = j
                    sched[j, s] = -i
                end
            end
        end
    end
    return sched
end

"""
    parse_solution(model)

Retrieves a solution schedule from `model`.
"""
function parse_solution(model)
    sol = value.(model[:x])
    sched = Array{Int64}(undef, size(sol, 1), size(sol, 3))
    return parse_solution!(sched, sol, model)
end
