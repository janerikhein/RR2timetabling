
"""
    swapHomes!(S::Array{Int, 2}, i::Int, j::Int)

swaps home/away roles of two teams `i` and `j` in the schedule.
"""
function swap_homes!(S::Array{Int, 2}, i::Int, j::Int)
    slot1 = 0
    slot2 = 0
    @inbounds for col = 1:size(S, 2)
        if S[i, col] == j
            slot1 = col
        end
        if S[j, col] == i
            slot2 = col
        end
        if slot1 !== 0 && slot2 !== 0
            break
        end
    end
    @assert slot1 != 0 && slot2 != 0
    S[i, slot1], S[i, slot2] = S[i, slot2], S[i, slot1]
    S[j, slot1], S[j, slot2] = S[j, slot2], S[j, slot1]
end

"""
    swapRounds!(S::Array{Int, 2}, k::Int, l::Int)

Swaps the two rounds `k` and `l` in the schedule.
"""
function swap_rounds!(S::Array{Int, 2}, k, l)
    @inbounds for row = 1:size(S, 1)
        S[row, k], S[row, l] = S[row, l], S[row, k]
    end
end

"""
    slot_swap_teams!(S::Array{Int, 2}, i::Int, j::Int, k::Int)

Swaps games of the two teams `i` and `j` in a single slot `k`. Does `not`
produce a RR2-tournament schedule.
"""
@inline function slot_swap_teams!(S::Array{Int, 2}, i::Int, j::Int, k::Int)
    oppi, oppj = S[i, k], S[j, k]
    if sign(oppi) == sign(oppj)
        S[abs(oppi), k], S[abs(oppj), k] = S[abs(oppj), k], S[abs(oppi), k]
    else
        S[abs(oppi), k], S[abs(oppj), k] = -S[abs(oppj), k], -S[abs(oppi), k]
    end
    S[i, k], S[j, k] = S[j, k], S[i, k]
end


"""
    swapTeams!(S::Array{Int, 2}, i::Int, j::Int)

Swaps games of two teams `i` and `j` in all slots except the two slots in which
they meet.
"""
function swap_teams!(S::Array{Int, 2}, i::Int, j::Int)
    for col = 1:size(S, 2)
        if S[i, col] == j || S[j, col] == i
            continue
        end
        slot_swap_teams!(S, i, j, col)
    end
end

"""
    partialSwapRounds!(S::Array{Int, 2}, i::Int, k::Int, l::Int)

Swaps games in two slots `k` and `l` for a team `i` and updates the schedule in
a deterministic way to produce a double round-robin tournament.
"""
function partial_swap_rounds!(S::Array{Int, 2}, i::Int, k::Int, l::Int)
    swapped = 2^i
    team1, team2 = abs(S[i, k]), abs(S[i, l])
    idx1, idx2 = k, l
    S[i, k], S[i, l] = S[i, l], S[i, k]
    while true
        cont = false
        if swapped & (1 << team1) === 0
            S[team1, k], S[team1, l] = S[team1, l], S[team1, k]
            swapped += 2^team1
            team1 = abs(S[team1, idx1])
            idx1 = idx1 == k ? l : k
            cont = true
        end
        if swapped & (1 << team2) === 0
            S[team2, k], S[team2, l] = S[team2, l], S[team2, k]
            swapped += 2^team2
            team2 = abs(S[team2, idx2])
            idx2 = idx2 == k ? l : k
            cont = true
        end
        cont || break
    end
end

"""
    partial_swap_teams!(S::Array{Int, Int}, i::Int, j::Int, k::Int)

Swaps games of two teams `i` and `j` in a single slot `k` and updates the
schedule in a deterministic way to produce a double round-robin tournament.
"""
function partial_swap_teams!(S::Array{Int,2}, i::Int, j::Int, k::Int)
    swapped = 2^k
    slot_swap_teams!(S, i, j, k)
    match1, match2 = S[i, k], S[j, k]
    offset = 0
    fixedi, fixedj = false, false
    for col in Iterators.cycle(vcat(k+1:size(S, 2), 1:k))
        if !fixedi && S[i, col] == match1
            if swapped & (1 << col) !== 0
                fixedi = true
            else
                slot_swap_teams!(S, i, j, col)
                swapped += 2^col
                match1 = S[i, col]
            end
        end
        if !fixedj && S[j, col] == match2
            if swapped & (1 << col) !== 0
                fixedj = true
            else
                slot_swap_teams!(S, i, j, col)
                swapped += 2^col
                match2 = S[j, col]
            end
        end

        fixedi && fixedj && break
    end
end


"""
    getmoves(rr2::RR2)

Compute all possible neighbor moves for a given tournament. If the touramenent
is phased the neighborhood is restricted, only alowing moves that retain the
phased structure.

# Parameters
-`rr2::RR2`- Tournament structure data

# Returns
-`Vector{Move}` - Vector of
    moves. A move tuple contains the neighborhood move function to
    be called (given as a `Symbol`), followed by the call arguments.
"""
function getmoves(rr2::RR2)
    sh = [(:swap_homes!, i, j) for (i, j) in combinations(1:rr2.nteams, 2)]
    st = [(:swap_teams!, i, j) for (i, j) in combinations(1:rr2.nteams, 2)]

    if rr2.isphased
        pst = Vector{Move}()
        sr = vcat(
            [
                (:swap_rounds!, k, l) for
                (k, l) in combinations(1:rr2.nslots÷2, 2)
            ],
            [
                (:swap_rounds!, k, l) for
                (k, l) in combinations(rr2.nslots÷2+1:rr2.nslots, 2)
            ],
        )
        psr = vcat(
            vec([
                (:partial_swap_rounds!, i, k, l) for i = 1:rr2.nteams,
                (k, l) in combinations(1:rr2.nslots÷2, 2)
            ]),
            vec([
                (:partial_swap_rounds!, i, k, l) for i = 1:rr2.nteams,
                (k, l) in combinations(rr2.nslots÷2+1:rr2.nslots, 2)
            ]),
        )
    else
        pst = vec([
            (:partial_swap_teams!, i, j, k) for
            (i, j) in combinations(1:rr2.nteams, 2), k = 1:rr2.nslots
        ])
        sr = [
            (:swap_rounds!, k, l) for (k, l) in combinations(1:rr2.nslots÷2, 2)
        ]
        psr = vec([
            (:partial_swap_rounds!, i, k, l) for i = 1:rr2.nteams,
            (k, l) in combinations(1:rr2.nslots, 2)
        ])
    end
    vcat(sh, sr, st, psr, pst)
end

"""
    move!(sched::Matrix{Int}, mv::Move)

Apply a move `mv` to a schedule `sched`.
"""
move!(sched::Matrix{Int}, mv::Move) = eval(mv[1])(sched, mv[2:end]...)
