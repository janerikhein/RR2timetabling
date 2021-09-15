"""
Initializes the `TabuSearchOptimizer` and the `SearchState`.
"""
function initialize_search(
    inst::Instance,
    start::Matrix{Int},
    maxiter::Int,
    tabulength::Int,
)
    cons = sort!(constraints(inst), by = !ishard)
    pen = [c.pen for c in cons]
    entries = conentries(inst.rr2, cons)
    tso = TabuSearchOptimizer(maxiter, tabulength, inst.rr2, cons, entries)
    val = 0
    offset = [eval_constraint(start, c) for c in cons]
    for (i, c) in enumerate(cons)
        @assert c.pen != 0 || offset[i] == 0
        val += offset[i] * c.pen
    end
    moves = getmoves(inst.rr2)
    tabu = Dict((m, 0) for m in moves)
    blocked = Dict((m, 0) for m in getmoves(inst.rr2))

    state = SearchState(
        0,
        val,
        copy(start),
        val,
        copy(start),
        0,
        offset,
        moves,
        tabu,
        blocked,
        pen,
        zeros(Int, size(start)),
        false
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
    difbuf::BitMatrix,
)
    state.it += 1
    bestmv = nothing
    bestval = typemax(Int)

    infeas = falses(size(schedbuf))
    for idx in 1:length(tso.constraints)
        if state.offset[idx] > 0
            infeas = infeas .| @view tso.conentries[:,:,idx]
        end
    end

    shuffle!(state.moves)

    for mv in state.moves
        # identical move already evaluated this iteration
        state.blocked[mv] == state.it && continue

        # reset buffers
        offsetbuf .= state.offset
        schedbuf .= state.sched

        move!(schedbuf, mv)
        difbuf .= (schedbuf .!= state.sched)
        if !state.diversify && !any(difbuf .& infeas)
            continue
        end

        # evaluate move and block equivalent moves
        isfeas = eval_move!(
            offsetbuf,
            schedbuf,
            difbuf,
            state.sched,
            tso.constraints,
            tso.conentries,
        )
        block_equiv!(state.blocked, mv, difbuf, state.it)
        !isfeas && continue
        ncons = length(tso.constraints)
        mvval = sum(state.pen[i] * offsetbuf[i] for i = 1:ncons)
        if state.diversify
            mvval += sum(state.fr[c] for c in findall(difbuf))
        end

        if state.tabu[mv] >= state.it # move is marked tabu
            if !state.diversify && mvval < state.optval # global improvement (aspiration)
                bestmv, bestval = mv, mvval
                break
            end
        else
            if mvval < state.val
                bestmv, bestval = mv, mvval # local improvement
                break
            elseif mvval < bestval
                bestmv, bestval = mv, mvval
            end
        end
    end

    if bestmv === nothing # no feasible neighbor solution found
        return false
    end
    @show bestmv
    schedbuf .= state.sched
    move!(schedbuf, bestmv)
    difbuf .= (schedbuf .!= state.sched)
    state.fr .+= difbuf
    eval_move!(
        state.offset,
        schedbuf,
        difbuf,
        state.sched,
        tso.constraints,
        tso.conentries,
    )
    move!(state.sched, bestmv)
    state.val = bestval
    if !state.diversify && state.val < state.optval
        state.optsched .= state.sched
        state.optval = state.val
        state.lastimprov = state.it
    end

    if state.val == state.optval
        state.optsched .= state.sched
    end
    # mark bestmv and equivalent moves as tabu
    block_equiv!(state.tabu, bestmv, difbuf, state.it + tso.tabulength)
    return true
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
            if c.pen == 0 && offsetbuf[i] > 0
                return false
            end
        end
    end
    return true
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
