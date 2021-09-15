"""
    pgba_heuristic(inst::Instance, maxbreaks::Int)

Finds a heuristic start solution satisfying all hard contraints of a given
instance making use of a pattern generating decomposition scheme.

# Arguments
- `inst::Instance` - the instance to be solved
- `maxbreaks::Int` - limits the number of patterns generated by setting an upper
    bound on the number of breaks a pattern is allowed to have

# Returns
- `Matrix{Int}` - a feasible schedule
"""
function pgba_heuristic(
    inst::Instance,
    maxbreaks::Int,
    schedbuf::Matrix{Int};
    gamtimelimit = GRB_INFINITY,
)

    nteams, nslots = inst.rr2.nteams, inst.rr2.nslots

    # pattern data
    pid = 0
    patternidx = Dict{Int,Int}()
    team_to_patidx = [Vector{Int}() for i = 1:nteams]
    patternCons, patternSetCons, genCons = group_constraints(inst; f = ishard)
    lazycutbuf = Vector{PGBACut}()
    for nbreaks = 0:maxbreaks
        #adding new patterns
        for i = 1:nteams
            pat = generate_patterns(inst.rr2, nbreaks, patternCons[i])
            for p in pat
                idx = get!(patternidx, p) do
                    pid += 1
                end
                push!(team_to_patidx[i], patternidx[p])
            end
        end
        # model is always infeasible if no patterns with breaks are included
        nbreaks == 0 && continue
        # building model
        patterns = sort!(collect(keys(patternidx)), by = k -> patternidx[k])
        model = pgba_model(inst.rr2, patterns, team_to_patidx, patternSetCons)
        x = model[:x]
        # adding previously generated lazy cuts to model
        for cut in lazycutbuf
            @constraint(
                model,
                sum(x[a[1], a[2]] for a in cut.assignments) <= cut.ub
            )
        end
        # solve with callback
        MOI.set(model, MOI.RawParameter("LazyConstraints"), 1)
        MOI.set(
            model,
            Gurobi.CallbackFunction(),
            (cb_data, cb_where) -> pgba_callback(
                inst.rr2,
                model,
                patterns,
                genCons,
                lazycutbuf,
                schedbuf,
                gamtimelimit,
                cb_data,
                cb_where,
            ),
        )
        optimize!(model)

        # retrieve solution if feasible solution found
        if termination_status(model) == MOI.OPTIMAL
            return schedbuf
        end
    end
end

"""
    pgba_model(
        rr2::RR2,
        patterns::Vector{Int64},
        teamToPatidx::Vector{Vector{Int64}},
        cons::Vector{T},
    ) where {T<:Constraint}

Builds the master problem of the PGBA decomposition.

# Arguments
- `rr2::RR2` - tournament format data
- `patterns::Vector{Int64}` - each pattern is represented as a `Int64` where the
    `s`-th bit indicates if the pattern has a home game in slot `s`.
- `teamToPatidx::Vector{Vector{Int64}}` - `teamToPatidx[i]` is a vector of
    indices of patterns feasible for team `i` in `patterns`.
- `cons::Vector{T}` - Vector of patternset constraints

# Model Description

The model uses binary variables `x[i,p]` indicating if team `i` gets
assigned pattern `patterns[p]`.
"""
function pgba_model(
    rr2::RR2,
    patterns::Vector{Int64},
    teamToPatidx::Vector{Vector{Int64}},
    cons::Vector{T},
) where {T<:Constraint}

    model = direct_model(Gurobi.Optimizer(GRB_ENV))
    npatterns = length(patterns)
    nteams, nslots = rr2.nteams, rr2.nslots

    # 1 if pattern at index p has home game in slot s, 0 otw.
    h(p, s) = patterns[p] >>> (s - 1) & 1

    @variable(model, x[1:nteams, 1:npatterns], Bin)

    # each team gets assigned exactly one of its feasible patterns
    @constraint(model, [i = 1:nteams], sum(x[i, teamToPatidx[i]]) == 1)
    @constraint(model, [i = 1:nteams], sum(x[i, :]) == 1)

    # same number of home and away games in each slot
    @constraint(
        model,
        [s = 1:nslots],
        sum(sum(h(p, s) * x[i, p] for p in 1:npatterns) for i in 1:nteams) ==
        nteams ÷ 2
    )

    for con in cons
        add_patternset_con!(model, con, nteams, patterns)
    end

    return model
end

"""
    pgba_callback(
        rr2::RR2,
        model::Model,
        patterns::Vector{Int},
        gencons::Vector{Constraint},
        lazycuts::Vector{PGBACut},
        schedbuf::Matrix{Int},
        useiis::Bool,
        cb_data,
        cb_where::Cint,
    )

Callback function used to generate the Gurobi callback for the master problem.

# Arguments

- `rr2::RR2`- tournament format data
- `model::Model`- the direct Gurobi model
- `patterns::Vector{Int64}` - each pattern is represented as a `Int64` where the
    `s`-th bit indicates if the pattern has a home game in slot `s`.
- `gencons::Vector{Constraint}`- generic constraints of the instance
- `lazycuts::Vector{PGBACut}` - stores all cuts generated in callbacks
- `schedbuf::Matrix{Int}` - solution schedule buffer
- `cb_data`, `cb_where::Cint`- Gurobi callback arguments
"""
function pgba_callback(
    rr2::RR2,
    model::Model,
    patterns::Vector{Int},
    gencons::Vector{Constraint},
    lazycuts::Vector{PGBACut},
    schedbuf::Matrix{Int},
    gamtimelimit::Float64,
    cb_data::Gurobi.CallbackData,
    cb_where::Cint,
)
    cb_where != GRB_CB_MIPSOL && return
    x = model[:x]
    patternset, patidx = cb_load_solution(model, patterns, cb_data, cb_where)

    # check meeting ranges and add cuts accordingly. if cuts are added return.
    infeasPairs = meeting_ranges_check(model, rr2, patternset)
    if !isempty(infeasPairs)
        for ij in infeasPairs
            addcut!(model, lazycuts, cb_data, ij, patidx, true)
        end
        return
    end

    # check pattern diversity with incr. subsetsize. Once cut is added return.
    for subsetsize = 3:rr2.nteams
        infeasTeams =
            pattern_diversity_check(model, rr2, subsetsize, patternset)
        if infeasTeams !== nothing
            addcut!(model, lazycuts, cb_data, infeasTeams, patidx, false)
            @show subsetsize
            return
        end
    end

    val = assign_games!(
        schedbuf,
        rr2,
        patternset,
        gencons;
        timelimit = gamtimelimit,
        usestart = false,
    )

    if val == GRB_INFINITY
        addcut!(model, lazycuts, cb_data, 1:rr2.nteams, patidx, false)
    end
end

"""
    get_patterns_from_breaks(nslots::Int, breaks::Vector{Int})

Computes the two potential patterns of length `nslots` with breaks in slots
given by `breaks`.
"""
function get_patterns_from_breaks(nslots::Int, breaks::Vector{Int})
    pat1 = sum(2^i for i = 0:2:nslots-1)
    pat2 = sum(2^i for i = 1:2:nslots-1)
    for bridx = 1:2:length(breaks)
        idx1 = breaks[bridx]
        idx2 = bridx == length(breaks) ? nslots : breaks[bridx+1] - 1
        pat1 = flipbits(pat1, idx1, idx2)
        pat2 = flipbits(pat2, idx1, idx2)
    end
    return pat1, pat2
end

"""
    generate_patterns(rr2::RR2, numbreaks::Int, patterncons::Vector{Constraint})

Generate all feasible patterns with a given number of breaks wrt. to some
pattern constraints.

# Arguments
- `rr2::RR2` - tournament format data
- `numbreaks::Int` - number of breaks each generated pattern must have
- `patterncons::Vector{Constraint}` - pattern constraints to be satisfied

# Returns
- `Set{Int64}` - the set of patterns satisfying `patterncons` where each pattern
    is represented as a 64-bit integer
"""
function generate_patterns(
    rr2::RR2,
    numbreaks::Int,
    patterncons::Vector{T},
) where {T<:Constraint}

    patterns = Set{Int64}()
    for breaks in combinations(2:rr2.nslots, numbreaks)
        pat1, pat2 = get_patterns_from_breaks(rr2.nslots, breaks)
        # only check pat1 since pat2 is simply bitwise complement
        if isPattern(rr2.nslots, pat1)
            # validating pattern constraints for pat1
            if all(c -> validate_pattern_constraint(rr2, pat1, c), patterncons)
                push!(patterns, pat1)
            end
            # validating pattern constraints for pat2
            if all(c -> validate_pattern_constraint(rr2, pat2, c), patterncons)
                push!(patterns, pat2)
            end
        end
    end
    return patterns
end


"""
    add_patternset_con!(
        model::Model,
        con::T,
        nteams::Int,
        patterns::Vector{Int}
    ) where {T <: Union{CAcon, BRcon}

Adds a patternset constraint `con` to `model`.
"""
function add_patternset_con! end

function add_patternset_con!(
    model::Model,
    con::CAcon,
    nteams::Int,
    patterns::Vector{Int},
)
    @assert con.mode != MODE_HA
    x = model[:x]
    if con.mode == MODE_H
        @constraint(
            model,
            sum(
                count_ones(patterns[p] & con.slots.binrepr) * x[i, p] for
                i in con.teams1, p in 1:length(patterns)
            ) <= con.ub
        )
    else # mode_A
        @constraint(
            model,
            sum(
                count_ones(~patterns[p] & con.slots.binrepr) * x[i, p] for
                i in con.teams1, p in 1:length(patterns)
            ) <= con.ub
        )
    end
end

function add_patternset_con!(
    model::Model,
    con::BRcon,
    nteams::Int,
    patterns::Vector{Int},
)
    @assert con.mode == MODE_HA

    x = model[:x]
    @constraint(
        model,
        sum(
            (
                count_ones(
                    patterns[p] & (patterns[p] << 1) & con.slots.binrepr,
                ) + count_ones(
                    ~patterns[p] & (~patterns[p] << 1) & con.slots.binrepr,
                )
            ) * x[i, p] for i in con.teams, p = 1:length(patterns)
        ) <= con.ub
    )
end

function pattern_diversity_model(
    rr2::RR2,
    model::Model,
    patternset::Vector{Int64},
    subsetsize::Int,
    slotrange::UnitRange{Int64},
)

    ps(j, s) = patternset[j] >>> (s - 1) & 1
    ubm = direct_model(Gurobi.Optimizer(GRB_ENV_2))
    MOI.set(ubm, MOI.RawParameter("OutputFlag"), 0)
    @variable(ubm, b[slotrange] >= 0)
    @variable(ubm, a[1:rr2.nteams], Bin)
    @variable(ubm, d[slotrange], Bin)
    @constraint(ubm, sum(a) == subsetsize)

    @constraint(
        ubm,
        [s in slotrange],
        b[s] - sum(ps(j, s) * a[j] for j = 1:rr2.nteams) + subsetsize -
        subsetsize * d[s] >= 0
    )
    @constraint(
        ubm,
        [s in slotrange],
        b[s] - sum((1 - ps(j, s)) * a[j] for j = 1:rr2.nteams) +
        subsetsize * d[s] >= 0
    )
    @objective(ubm, Min, sum(b))
    return ubm
end

function pattern_diversity_check(model, rr2, subsetsize, patternset)
    if rr2.isphased
        for phase = 1:2
            slots = phase == 1 ? (1:rr2.nslots÷2) : (rr2.nslots÷2+1:rr2.nslots)
            ubm = pattern_diversity_model(
                rr2,
                model,
                patternset,
                subsetsize,
                1:rr2.nslots÷2,
            )
            optimize!(ubm)
        end
        if objective_value(ubm) < (subsetsize * (subsetsize - 1)) / 2
            a = ubm[:a]
            return filter(i -> value(a[i]) ≈ 1, 1:rr2.nteams)

        end
    else
        ubm = pattern_diversity_model(
            rr2,
            model,
            patternset,
            subsetsize,
            1:rr2.nslots,
        )
        optimize!(ubm)
        if objective_value(ubm) < (subsetsize * (subsetsize - 1))
            a = ubm[:a]
            return filter(i -> value(a[i]) ≈ 1, 1:rr2.nteams)
        end
    end
    return nothing
end

function meeting_ranges_check(model, rr2, patternset)
    infeasPairs = Vector{Tuple{Int,Int}}()
    for (i, j) in combinations(1:rr2.nteams, 2)
        hg = patternset[i] & ~patternset[j] # slots for game (i,j)
        ag = ~patternset[i] & patternset[j] # slots for game (j,i)
        phase1 = 2^(rr2.nslots ÷ 2) - 1
        phase2 = phase1 << (rr2.nslots ÷ 2)
        if hg == 0 || ag == 0
            push!(infeasPairs, (i, j))
        elseif rr2.isphased
            if hg & phase1 == ag & phase1 == 0 ||
               hg & phase2 == ag & phase2 == 0
                push!(infeasPairs, (i, j))
            end
        end
    end
    return infeasPairs
end

"""
    cb_load_solution(model::Model, cb_data, cb_where::Cint)

Loads solution values for the PGBA master `model` during a callback.

# Returns
`Vector{Int}`, `Vector{Int}` - patternset and pattern index vectors
"""
function cb_load_solution(model::Model, patterns::Vector{Int}, cb_data, cb_where::Cint)
    @assert cb_where == GRB_CB_MIPSOL

    Gurobi.load_callback_variable_primal(cb_data, cb_where)
    x_val = callback_value.(Ref(cb_data), model[:x])
    patternset = Vector{Int64}(undef, size(x_val, 1))
    idx = Vector{Int64}(undef, size(x_val, 1))
    for p = 1:size(x_val, 2)
        for i = 1:size(x_val, 1)
            if x_val[i, p] ≈ 1
                patternset[i] = patterns[p]
                idx[i] = round(Int, p)
            end
            continue
        end
    end
    return patternset, idx
end



function addcut!(model, lazycuts, cut::PGBACut, cb_data)
    x = model[:x]
    lzcon = @build_constraint(sum(x[a[1], a[2]] for a in cut.assignments) <= cut.ub)
    push!(lazycuts, cut)
    MOI.submit(model, MOI.LazyConstraint(cb_data), lzcon)
end


"""
Adds a cut to model requiring at least one team in teams to have a different
pattern assigned to them
"""
function addcut!(model, lazycuts, cb_data, teams, patidx, team_specific::Bool)
    x = model[:x]
    if team_specific
        cut = PGBACut(length(teams)-1, collect((i,patidx[i]) for i in teams))
    else
        cut = PGBACut(length(teams)-1, vcat([(j, patidx[i]) for i in teams, j=1:size(x,1)]...))
    end
    addcut!(model, lazycuts, cut, cb_data)
end