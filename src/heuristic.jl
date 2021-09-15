"""
    heur_solution(inst::Instance)

Computes a heuristic start solution schedule, satisfying all hard constraints
of `inst`.

# Arguments
-`inst::Instance` - the problem instance

# Description
The heuristic first generates a start solution, satisfying all base constraints.
Then a Tabu Search algorithm is used to minimize the penalty induced by the hard
constraints of `inst`.
"""
function heur_solution(inst::Instance; maxbreaks = 3, maxiter = 200, maxdiviter = 15, tabulength = 20, gamtimelimit = 60.0, maxphases=10, start=nothing)
    # compute initial solution satisfyng all pattern & patternset constraints
    if start===nothing
        schedbuf = Matrix{Int}(undef, inst.rr2.nteams, inst.rr2.nslots)
        start = pgba_heuristic(inst, maxbreaks, schedbuf; gamtimelimit = gamtimelimit)
    end
    start::Matrix{Int}
    patcons, patsetcons, gencons = group_constraints(inst; f = ishard)

    # adding phased constraints as game constraints
    if inst.rr2.isphased
        slots = IdxSet(2^(inst.rr2.nslots÷2)-1)
        for i = 1:inst.rr2.nteams
            for j = i+1:inst.rr2.nteams
                push!(gencons, GAcon(slots, 1, 1, [(i,j), (j,i)], 0))
            end
        end
        rr2 = RR2(inst.rr2.nteams, inst.rr2.nslots, false)
    else
        rr2 = inst.rr2
    end

    relinst = relaxed(Instance(rr2, vcat(vcat(patcons...), patsetcons, gencons)))
    tso, state = initialize_search(relinst, start, maxiter, tabulength)

    patconidx = filter(
        i -> is_pattern_constraint(inst.rr2, tso.constraints[i]),
        1:length(tso.constraints),
    )
    patsetconidx = filter(
        i -> is_pattern_set_constraint(inst.rr2, tso.constraints[i]),
        1:length(tso.constraints),
    )
    offsetbuf, schedbuf, difbuf = init_buffers(tso, state)
    is_feas = trues(inst.rr2.nteams)
    lastfeas = [trues(tabulength) for i in 1:inst.rr2.nteams]
    for phc = 1:maxphases
        while state.lastimprov + maxiter >= state.it
            @show state.it, state.val
            search_step!(tso, state, offsetbuf, schedbuf, difbuf)
            if state.optval == 0
                return state.optsched
            end
        end
        println("diversifying")
        state.diversify = true
        state.sched .= state.optsched
        state.val = 0
        state.lastimprov = state.it
        for it = 1:maxdiviter
            search_step!(tso, state, offsetbuf, schedbuf, difbuf)
        end
        state.val = sum(state.offset[i]*state.pen[i] for i in 1:length(state.offset))
        @show state.val, state.it, state.lastimprov, maxiter
        state.diversify = false
    end
end

"""
    relaxed(inst::Instance)

Relaxes a given instance. Returns a copy of `inst` with all soft constraints
deleted and all hard constraints markd as soft wth a penalty of 1.
"""
function relaxed(inst::Instance)
    cons = Vector{Vector{Constraint}}()
    for field in (:constrCA, :constrBR, :constrGA, :constrFA, :constrSE)
        temp = Vector{Constraint}()
        for c in getfield(inst, field)
            if c.pen == 0
                push!(temp, @set c.pen = 1)
            end
        end
        push!(cons, temp)
    end
    return Instance(inst.rr2, cons...)
end

"""
    assign_games(
        sched::Matrix{Int},
        rr2::RR2,
        patternset::Vector{Int},
        gencons::Vector{Constraint},
    )

Solves the game assignment subproblem, minimizing the offset penalties of
`gencons`. The schedule `sched` is used as a start solution if
`use_start == true` and gets overriden by the optimal solution.
"""
function assign_games!(
    sched::Matrix{Int},
    rr2::RR2,
    patternset::Vector{Int},
    gencons::Vector{Constraint};
    timelimit = GRB_INFINITY,
    usestart = true,
)
    @show timelimit
    gam = game_alloc_model(rr2, patternset, gencons; feas_relax = true)
    x = gam[:x]
    if usestart
        @assert is_schedule(rr2, sched)
        for s = 1:size(sched, 2)
            for i = 1:size(sched, 1)
                if sched[i, s] > 0
                    set_start_value(x[i, sched[i, s], s], 1)
                end
            end
        end
    end
    #MOI.set(gam, MOI.RawParameter("OutputFlag"), 0)
    MOI.set(
        gam,
        Gurobi.CallbackFunction(),
        (cb_data, cb_where) -> gam_callback(gam, cb_data, cb_where, timelimit),
    )
    optimize!(gam)
    status = termination_status(gam)
    if status == MOI.OPTIMAL || status == MOI.INTERRUPTED
        parse_solution!(sched, value.(x), gam)
        return objective_value(gam)
    else
        return GRB_INFINITY
    end
end

"""
    gam_callback(timelimit, cb_data, cb_where)

Callback used to stop GAM if the time limit is exceeded and a feasible solution
was found
"""
function gam_callback(model, cb_data, cb_where, timelimit)
    if cb_where == GRB_CB_MIP
        time, obj = Ref{Cdouble}(), Ref{Cdouble}()
        GRBcbget(cb_data, cb_where, GRB_CB_RUNTIME, time)
        GRBcbget(cb_data, cb_where, GRB_CB_MIP_OBJBST, obj)
        if time[] > timelimit && obj[] < GRB_INFINITY
            GRBterminate(backend(model))
        end
    end
end
