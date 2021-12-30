"""
Initializes the `TabuSearchOptimizer` and the `SearchState`.
"""
function initialize_search(
    rr2::RR2,
    cons::Vector{Constraint},
    pen::Vector{Int},
    start::Matrix{Int},
    maxiter::Int,
    tabulength::Int,
    infeaslength::Int,
)
    entries = conentries(rr2, cons)
    tso = TabuSearchOptimizer(
        maxiter,
        tabulength,
        infeaslength,
        rr2,
        cons,
        entries,
    )
    hpen, spen = 0, 0
    offset = [eval_constraint(start, c) for c in cons]
    for (i, c) in enumerate(cons)
        if c.pen == 0
            hpen += offset[i] * pen[i]
        else
            spen += offset[i] * pen[i]
        end
    end
    moves = getmoves(rr2)

    tabu = Dict((m, 0) for m in moves)
    blocked = Dict((m, 0) for m in getmoves(rr2))
    optval = hpen == 0 ? spen : typemax(Int)
    state = SearchState(
        0, #it
        hpen,
        spen,
        copy(start), #sched
        optval,
        copy(start), #optsched
        0, #lastimprov
        offset,
        moves,
        tabu,
        blocked,
        pen,
    )
    return tso, state
end

"""
Initialize buffers used in the Tabu Search. `offsetbuf` is used for storing
constraint offsets, `schedbuf` is used as a schedule buffer and `difbuf` is used
for move evaluation.
"""
function init_buffers(tso::TabuSearchOptimizer, state::SearchState)
    offsetbuf = copy(state.offset)
    schedbuf = similar(state.sched)
    difbuf = BitMatrix(undef, tso.rr2.nteams, tso.rr2.nslots)
    return offsetbuf, schedbuf, difbuf
end

"""
    search_step!(
        tso::TabuSearchOptimizer,
        state::SearchState,
        offsetbuf::Vector{Int},
        schedbuf::Matrix{Int},
        difbuf::BitMatrix,
    )

Selects a suitable neighbor solution and moves to it.
"""
function search_step!(
    tso::TabuSearchOptimizer,
    state::SearchState,
    offsetbuf::Vector{Int},
    schedbuf::Matrix{Int},
    difbuf::BitMatrix;
    accept_infeas = false,
)
    state.it += 1
    opt_feas_mv, opt_infeas_mv = nothing, nothing
    opt_infeas_hpen, opt_infeas_spen = typemax(Int), typemax(Int)
    opt_feas_spen = typemax(Int)

    # compute constrained entries of violated constraints
    infeas = falses(size(schedbuf))
    for idx in 1:length(tso.constraints)
        if state.offset[idx] > 0
            infeas = infeas .| @view tso.conentries[:,:,idx]
        end
    end
    if !any(infeas)
        infeas .= true
    end

    # randomize move order
    shuffle!(state.moves)

    for mv in state.moves
        # identical move already evaluated this iteration
        state.blocked[mv] == state.it && continue
        # reset buffers
        offsetbuf .= state.offset
        schedbuf .= state.sched

        move!(schedbuf, mv)
        # check if move changes constrained entry
        difbuf .= (schedbuf .!= state.sched)
        !any(difbuf .& infeas) && continue

        # evaluate move and block equivalent moves
        eval_move!(
            offsetbuf,
            schedbuf,
            difbuf,
            state.sched,
            tso.constraints,
            tso.conentries,
        )
        block_equiv!(state.blocked, mv, difbuf, state.it)

        hpen, spen = 0, 0
        hpenold, spenold = 0, 0
        for i in 1:length(tso.constraints)
            if tso.constraints[i].pen == 0
                hpen += state.pen[i] * offsetbuf[i]
                hpenold += state.pen[i] * state.offset[i]
            else
                spen += state.pen[i] * offsetbuf[i]
                spenold += state.pen[i] * state.offset[i]
            end
        end

        if state.tabu[mv] >= state.it # move is marked tabu
            # aspiration
            if hpen == 0 && spen < state.optval
                opt_feas_mv, opt_feas_spen = mv, spen
                break
            end
        else
            if hpen == 0
                if spen < spenold
                    opt_feas_mv, opt_feas_spen = mv, spen
                    break
                elseif spen < opt_feas_spen
                    opt_feas_mv, opt_feas_spen = mv, spen
                end
            else
                if (spen < spenold && hpen <= hpenold) || (spen <= spenold && hpen < hpenold)
                    opt_infeas_mv = mv
                    opt_infeas_hpen, opt_infeas_spen = hpen, spen
                    break
                elseif hpenold == 0 && (spen, hpen) < (opt_infeas_spen, opt_infeas_hpen)
                    opt_infeas_mv = mv
                    opt_infeas_hpen, opt_infeas_spen = hpen, spen
                elseif hpenold > 0 && (hpen, spen) < (opt_infeas_hpen, opt_infeas_spen)
                    opt_infeas_mv = mv
                    opt_infeas_hpen, opt_infeas_spen = hpen, spen
                end
            end
        end
    end

    if opt_feas_mv === nothing # no feasible neighbor solution found
        @assert opt_infeas_mv !== nothing
        bestmv = opt_infeas_mv
        state.hpen, state.spen = opt_infeas_hpen, opt_infeas_spen
    else
        if accept_infeas && opt_infeas_spen < opt_feas_spen
            @assert opt_infeas_mv !== nothing
            bestmv = opt_infeas_mv
            state.hpen, state.spen = opt_infeas_hpen, opt_infeas_spen
        else
            bestmv = opt_feas_mv
            state.hpen, state.spen = 0, opt_feas_spen
        end
    end
    schedbuf .= state.sched
    move!(schedbuf, bestmv)
    difbuf .= (schedbuf .!= state.sched)
    eval_move!(
        state.offset,
        schedbuf,
        difbuf,
        state.sched,
        tso.constraints,
        tso.conentries,
    )
    move!(state.sched, bestmv)

    if state.hpen == 0 && state.spen < state.optval
        state.optsched .= state.sched
        state.optval = state.spen
        state.lastimprov = state.it
    end

    # mark bestmv and equivalent moves as tabu
    block_equiv!(state.tabu, bestmv, difbuf, state.it + tso.tabulength)
    return bestmv
end

"""
    eval_move!(
        offset::Vector{Int},
        schedbuf::Matrix{Int},
        difbuf::BitMatrix,
        sched::Matrix{Int},
        constraints::Vector{Constraint},
        entries::BitArray{3},
    )

evaluates a move and updates constraint offsets. Only iterates through
constraints that are constraining (team, slot) entries changed by move.

# Arguments

`schedbuf::Matrix{Int}` - schedule buffer
`difbuf::BitMatrix` - schedule difference buffer
`sched::Matrix{Int}` - current solution
`constraints::Vector{Constraint}` - constraints of the model
`entries::BitArray{3}` - used for determining relevant constraints
"""
function eval_move!(
    offsetbuf::Vector{Int},
    schedbuf::Matrix{Int},
    difbuf::BitMatrix,
    sched::Matrix{Int},
    constraints::Vector{Constraint},
    entries::BitArray{3},
)
    for (i, c) in enumerate(constraints)
        if any(difbuf .& entries[:, :, i])
            offsetbuf[i] = eval_constraint(schedbuf, c)
        end
    end
end

"""
Marks equivalent moves to mv as blocked.
"""
function block_equiv!(
    blocked::Dict{Move,Int},
    mv::Move,
    dif::BitMatrix,
    it::Int,
)
    blocked[mv] = it
    if mv[1] == :partial_swap_teams!
        t1, t2, s = mv[2:end]
        for slot = 1:size(dif, 2)
            if dif[t1, slot]
                blocked[(:partial_swap_teams!, t1, t2, slot)] = it
            end
        end
    elseif mv[1] == :partial_swap_rounds!
        t, s1, s2 = mv[2:end]
        for team = 1:size(dif, 1)
            if dif[team, s1]
                blocked[(:partial_swap_rounds!, team, s1, s2)] = it
            end
        end
    end
end
