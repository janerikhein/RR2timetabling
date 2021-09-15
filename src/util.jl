# flips bits between ith and jth bit in n
flipbits(n, i, j) = n ⊻ (2^(j - i + 1) - 1) << (i - 1)

# checks if given pattern has same number of home and away games
isPattern(nslots, pat) =
    count_ones(pat) == count_zeros(pat) - sizeof(pat) * 8 + nslots

ishard(c::T) where {T<:Constraint} = c.pen == 0

@inline teams(c::Constraint) = c.teams
@inline teams(c::CAcon) = IdxSet(c.teams1.binrepr | c.teams2.binrepr)
@inline teams(c::GAcon) = unique(Iterators.flatten(c.meetings))

@inline patconteam(c::CAcon) = only(c.teams1) # get team for a patetrn con
@inline patconteam(c::BRcon) = only(c.teams)

@inline patsetconteams(c::CAcon) = c.teams1
@inline patsetconteams(c::BRcon) = c.teams

@inline slots(c::Constraint) = c.slots

@inline constraints(i::Instance) =
    vcat(i.constrCA, i.constrBR, i.constrGA, i.constrFA, i.constrSE)

# base case
@inline is_pattern_constraint(rr2::RR2, con::Constraint) = false

@inline is_pattern_constraint(rr2::RR2, con::BRcon) =
    is_pattern_constraint(con::BRcon)
@inline is_pattern_constraint(con::BRcon) = length(con.teams) == 1
@inline is_pattern_constraint(rr2::RR2, con::CAcon) =
    length(con.teams1) == 1 && length(con.teams2) == rr2.nteams

@inline is_pattern_set_constraint(rr2::RR2, con::Constraint) = false
@inline is_pattern_set_constraint(rr2::RR2, con::CAcon) =
    length(con.teams1) > 1 && length(con.teams2) == rr2.nteams
@inline is_pattern_set_constraint(rr2::RR2, con::BRcon) = length(con.teams) > 1

function patternset(sched::Matrix{Int})
    return [
        sum(2^(s - 1) for s in 1:size(sched, 2) if sign(sched[i, s]) == 1) for
        i = 1:size(sched, 1)
    ]
end

"""
    is_schedule(rr2::RR2, S::Array{Int, 2})

Checks if a schedule satisfies all base constraints.

# Arguments

- `rr2::RR2`: Format data of the problem instance
- `S::Array{Int,2}`: The potential schedule to be verified. Assuming teams `i`
    and `j` meet in slot `s` we have `S[i,s] == j`, if team `i` plays at home
    and `S[i,s] == -j` if team `i` plays away.

# Notes

This function is purely intended for debugging purposes.
"""
function is_schedule(rr2::RR2, S::Array{Int,2})
    @assert size(S, 1) == rr2.nteams && size(S, 2) == rr2.nslots

    for s = 1:rr2.nslots
        # exactly rr2.nteams÷2 teams play at home in each slot
        if count(x -> sign(x) == 1, S[:, s]) != rr2.nteams ÷ 2
            return false
        end
        # each team participates in each slot exactly once
        if sort(abs.(S[:, s])) != 1:rr2.nteams
            return false
        end
    end

    # each team plays against each other team once at home and once away
    for r = 1:rr2.nteams
        if sort(S[r, :]) !=
           vcat(-rr2.nteams:-r-1, -r+1:-1, 1:r-1, r+1:rr2.nteams)
            return false
        end
    end
    # check if first half of tournament is a single round robin
    if rr2.isphased
        for i = 1:rr2.nteams
            if sort(abs.(S[i, 1:rr2.nslots÷2])) != vcat(1:i-1, i+1:rr2.nteams)
                return false
            end
        end
    end

    return true
end
