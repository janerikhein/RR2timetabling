"""
    heur_solution(inst::Instance)

Computes a heuristic solution schedule, satisfying all hard constraints
minimizing the panalties

# Arguments
-`inst::Instance` - the problem instance

# Description
The heuristic first generates a start solution, satisfying all base constraints.
Then a Tabu Search algorithm is used to minimize the penalty induced by the hard
constraints of `inst`.
"""
function heur_solution(
    inst::Instance;
    maxbreaks = 3,
    maxiter = 1000,
    infeaslength = 5,
    tabulength = 20,
    lastfeaslength = 50,
    gamtimelimit = 1200.0,
    maxphases = 20,
    inithardpen = 100,
    relaxphased = true,
    start = nothing,
)
    # compute initial solution satisfyng all pattern & patternset constraints
    if start === nothing
        schedbuf = Matrix{Int}(undef, inst.rr2.nteams, inst.rr2.nslots)
        start = pgba_heuristic(
            inst,
            maxbreaks,
            schedbuf;
            gamtimelimit = gamtimelimit,
        )
    end
    start::Matrix{Int}

    if inst.rr2.isphased && relaxphased
        phasedcons = Vector{Constraint}()
        slots = IdxSet(2^(inst.rr2.nslots ÷ 2) - 1)
        for i = 1:inst.rr2.nteams
            for j = i+1:inst.rr2.nteams
                push!(phasedcons, GAcon(slots, 1, 1, [(i, j), (j, i)], 0))
            end
        end
        rr2 = RR2(inst.rr2.nteams, inst.rr2.nslots, false)
        inst = Instance(rr2, vcat(constraints(inst), phasedcons))
    end

    cons = constraints(inst)
    hardCons = ishard.(cons)
    #lastfeas = [trues(lastfeaslength) for c = 1:length(cons)]
    pen = Int.([ishard(c) ? inithardpen : c.pen for c in cons])

    tso, state =
        initialize_search(inst.rr2, cons, pen, start, maxiter, tabulength, infeaslength)
    offsetbuf, schedbuf, difbuf = init_buffers(tso, state)
    starttime = time()

    tabu = deepcopy(state.tabu)
    add_tabu = Vector{Move}()
    #allow_infeas = false
    for phc = 1:maxphases
#        allow_infeas = !allow_infeas
        state.lastimprov = state.it
        @show state.lastimprov, maxiter, state.it
#        @show allow_infeas
        while state.lastimprov + maxiter >= state.it
            search_step!(
                tso,
                state,
                offsetbuf,
                schedbuf,
                difbuf,
                accept_infeas = false,
            )
            # logging
            print("it ", state.it)
            print(" hpen ", state.hpen, " spen ", state.spen)
            print(" optval ", state.optval, " lastimprov ", state.lastimprov)
            print(" time ", round(time() - starttime, digits = 3), "\n")
            if state.hpen == 0 && state.spen == state.optval
                println(
                    "new best feasible solution found with value ",
                    state.optval,
                )
                println(state.optsched)
                if state.spen == 0
                    return state.optsched, state.optval
                end
                tabu = deepcopy(state.tabu)
                add_tabu = Vector{Move}()
            end
        end
        if state.optval == typemax(Int) # no feasible solution found
            return state.optsched, typemax(Int)
        end
        map!(i -> i + state.it - state.lastimprov, values(tabu))
        state.tabu = deepcopy(tabu)
        state.sched .= state.optsched
        state.hpen, state.spen = 0, state.optval
        for mv in add_tabu
            state.tabu[mv] = state.it + tabulength
        end
        println("accepting infeasible solutions...")
        for k = 1:infeaslength
            println("spen, hpen =", state.spen, " ", state.hpen)
            push!(
                add_tabu,
                search_step!(
                    tso,
                    state,
                    offsetbuf,
                    schedbuf,
                    difbuf,
                    accept_infeas = true,
                ),
            )
        end
    end
    return state.optsched, state.optval
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
