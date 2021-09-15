
"""
    eval_constraint(S::Array{Int,2}, c::CAcon)
    eval_constraint(S::Array{Int,2}, c::BRcon)
    eval_constraint(S::Array{Int,2}, c::GAcon)
    eval_constraint(S::Array{Int,2}, c::FAcon)
    eval_constraint(S::Array{Int,2}, c::SEcon)

evaluate the penalty for a constraint on a schedule.

# Arguments

- `S::Array{Int,2}`: The schedule to be evaluated on. Assuming teams `i` and `j`
    meet in slot `s` we have `S[i,s] == j`, if team `i` plays at home and
    `S[i,s] == -j` if team `i` plays away.
- `c::T where T <: Constraint`: The constraint to be evaluated.

# Returns

- `Int`: the offset of the constraint on the given schedule.
"""
function eval_constraint end

function eval_constraint(S::Array{Int,2}, c::CAcon)
    count = 0
    for i in c.teams1
        for s in c.slots
            match = S[i, s]
            if c.mode == MODE_HA || sign(match) == c.mode
                if abs(match) in c.teams2
                    count += 1
                end
            end
        end
    end
    return max(count - c.ub, 0)
end

function eval_constraint(S::Array{Int,2}, c::BRcon)
    count = 0
    for i in c.teams
        for s in c.slots
            s == 1 && continue
            if c.mode == MODE_HA || sign(S[i, s]) == c.mode
                if sign(S[i, s]) == sign(S[i, s-1])
                    count += 1
                end
            end
        end
    end
    return max(count - c.ub, 0)
end

function eval_constraint(S::Array{Int,2}, c::GAcon)
    count = 0
    for s in c.slots
        for i = 1:size(S, 1)
            if sign(S[i, s]) == 1 && (i, S[i, s]) in c.meetings
                count += 1
            end
        end
    end
    return max(max(c.lb - count, 0), max(count - c.ub, 0))
end

function eval_constraint(S::Array{Int,2}, c::FAcon)
    hgames = zeros(Int64, size(S, 1))
    maxdiff = zeros(Int64, size(S, 1), size(S, 1))
    lastslot = 0
    for s in c.slots
        for sl = lastslot+1:s
            for i = 1:size(S, 1)
                if sign(S[i, sl]) == 1
                    hgames[i] += 1
                end
            end
        end
        for i in c.teams
            for j in c.teams
                j <= i && continue
                if maxdiff[i, j] < abs(hgames[i] - hgames[j])
                    maxdiff[i, j] = abs(hgames[i] - hgames[j])
                end
            end
        end
        lastslot = s
    end
    count = 0
    for i in c.teams
        for j in c.teams
            j <= i && continue
            count += max(maxdiff[i, j] - c.ub, 0)
        end
    end
    return count
end

function eval_constraint(S::Array{Int,2}, c::SEcon)
    count = 0
    games = Array{Int64}(undef, size(S, 1), size(S, 1))
    for s = 1:size(S, 2)
        for i = 1:size(S, 1)
            if sign(S[i, s]) == 1
                games[i, S[i, s]] = s
            end
        end
    end
    for i in c.teams
        for j in c.teams
            j <= i && continue
            if abs(games[i, j] - games[j, i]) < c.lb
                count += c.lb - abs(games[i, j] - games[j, i]) + 1
            end
        end
    end
    return count
end

"""
    validate_pattern_constraint(rr2::RR2, pat::Int64, con::CAcon)
    validate_pattern_constraint(rr2::RR2, pat::Int64, con::BRcon)

Checks if a pattern constraint `con` is valid for a given pattern `pat`.
"""
function validate_pattern_constraint end

function validate_pattern_constraint(rr2::RR2, pat::Int64, con::CAcon)
    @assert is_pattern_constraint(rr2, con)

    if con.mode == MODE_H
        return count_ones(pat & con.slots.binrepr) <= con.ub
    else # mode_A
        return count_ones(flipbits(pat, 1, rr2.nslots) & con.slots.binrepr) <=
               con.ub
    end
end

function validate_pattern_constraint(rr2::RR2, pat::Int64, con::BRcon)
    @assert is_pattern_constraint(rr2, con)

    breakcount = 0
    if con.mode == MODE_H || con.mode == MODE_HA
        breakcount += count_ones((pat & pat << 1) & con.slots.binrepr)
    end
    if con.mode == MODE_A || con.mode == MODE_HA
        pat = flipbits(pat, 1, rr2.nslots)
        breakcount += count_ones((pat & pat << 1) & con.slots.binrepr)
    end
    return breakcount <= con.ub
end

"""
Groups `cons` into pattern and patternset constraints and adds them to
`patterncons` or `patternSetCons` respectivly.
"""
function group_constraints!(
    patternCons::Vector{Vector{Constraint}},
    patternSetCons::Vector{Constraint},
    cons::Vector{BRcon},
)
    for brc in cons
        if is_pattern_constraint(brc)
            team = trailing_zeros(brc.teams.binrepr) + 1
            push!(patternCons[team], brc)
        else
            push!(patternSetCons, brc)
        end
    end
end

"""
Groups `cons` into pattern, patternset and general constraints and adds them to
`patterncons`, `patternSetCons` or `genCons` respectivly.
"""
function group_constraints!(
    patternCons::Vector{Vector{Constraint}},
    patternSetCons::Vector{Constraint},
    genCons::Vector{Constraint},
    rr2::RR2,
    cons::Vector{CAcon},
)
    for cac in cons
        if is_pattern_constraint(rr2, cac)
            team = trailing_zeros(cac.teams1.binrepr) + 1
            push!(patternCons[team], cac)
        elseif is_pattern_set_constraint(rr2, cac)
            push!(patternSetCons, cac)
        else
            push!(genCons, cac)
        end
    end
end

"""
Adds `cons` to `genCons` and adds a relaxation into two capacity constraints to
`patternCons` or `patternSetCons` respectively.
"""
function group_constraints!(
    patternCons::Vector{Vector{Constraint}},
    patternSetCons::Vector{Constraint},
    genCons::Vector{Constraint},
    rr2::RR2,
    cons::Vector{GAcon},
)
    for gac in cons
        # special case in which GA con can be relaxed to two CA cons
        if gac.lb > 0
            teams1, teams2 = unique.(zip(gac.meetings...))

            # teams1 must play at least lb home games in slots
            cacon1 = CAcon(
                IdxSet(teams1),
                IdxSet(2^(rr2.nteams) - 1),
                gac.slots,
                length(teams1) * length(gac.slots) - gac.lb,
                MODE_A,
                0,
            )
            # teams2 must play at least lb away games in slots
            cacon2 = CAcon(
                IdxSet(teams2),
                IdxSet(2^(rr2.nteams) - 1),
                gac.slots,
                length(teams2) * length(gac.slots) - gac.lb,
                MODE_H,
                0,
            )
            length(teams1) == 1 ? push!(patternCons[only(teams1)], cacon1) :
                push!(patternSetCons, cacon1)
            length(teams2) == 1 ? push!(patternCons[only(teams2)], cacon2) :
                push!(patternSetCons, cacon2)
        end
        push!(genCons, gac)
    end
end

"""
    group_constraints(inst::Instance)

Groups all constraints of `inst` into pattern constraints, pattern-set
constraints and general constraints.
"""
function group_constraints(inst::Instance; f = x -> true)
    nteams = inst.rr2.nteams
    patternCons = [Vector{Constraint}() for i = 1:nteams]
    patternSetCons = Vector{Constraint}()
    genCons = Vector{Constraint}()
    group_constraints!(
        patternCons,
        patternSetCons,
        filter(f, inst.constrBR),
    )
    group_constraints!(
        patternCons,
        patternSetCons,
        genCons,
        inst.rr2,
        filter(f, inst.constrCA),
    )
    group_constraints!(
        patternCons,
        patternSetCons,
        genCons,
        inst.rr2,
        filter(f, inst.constrGA),
    )
    return patternCons, patternSetCons, genCons
end

"""
    conentries(rr2::RR2, cons::Vector{Constraint})

Generates a `BitArray{3}` where the index `(c,i,s)` indicates whether `cons[c]`
constrains games played by team `i` in slot `s`.
"""
function conentries(rr2::RR2, cons::Vector{Constraint})
    A = falses(rr2.nteams, rr2.nslots, length(cons))
    for (i, c) in enumerate(cons)
        constrained_entries!(view(A, :, :, i), rr2, c)
    end
    return A
end

"""
    constrained_entries(M::AbstractArray, rr2::RR2, con::T) where T <:Constraint

Compute all (team, slot) entries that are constrained by `con` and stores them
in `M` where M is either a `BitMatrix` or a `SubArray` of a BitMatrix.
"""
function constrained_entries! end

@inline function constrained_entries!(
    M::AbstractArray,
    rr2::RR2,
    con::Union{CAcon,GAcon},
)
    for i in teams(con)
        for s in slots(con)
            M[i, s] = 1
        end
    end
end

@inline function constrained_entries!(M::AbstractArray, rr2::RR2, con::BRcon)
    for i in teams(con)
        for s in slots(con)
            M[i, s] = 1
            if s > 1
                M[i, s-1] = 1
            end
        end
    end
end

@inline function constrained_entries!(M::AbstractArray, rr2::RR2, con::FAcon)
    maxslot = max(slots(con)...)
    for i in teams(con)
        M[i, 1:maxslot] .= 1
    end
end

@inline function constrained_entries!(M::AbstractArray, rr2::RR2, con::SEcon)
    for i in teams(con)
        M[i, :] .= 1
    end
end
