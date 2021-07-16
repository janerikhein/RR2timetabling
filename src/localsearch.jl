import DataStructures


function swapHomes!(S,i,j)
    slot1 = 0
    slot2 = 0
    @inbounds for col=1:size(S,2)
        if S[i,col] == j
            slot1 = col
        end
        if S[j,col] == i
            slot2 = col
        end
        if slot1!==0 && slot2!==0
            break
        end
    end
    @assert slot1 != 0 && slot2 != 0
    S[i,slot1], S[i,slot2] = S[i,slot2], S[i,slot1]
    S[j,slot1], S[j,slot2] = S[j,slot2], S[j,slot1]
end

function swapRounds!(S,i,j)
    @inbounds for row=1:size(S,1)
        S[row,i], S[row,j] = S[row,j], S[row,i]
    end
end

@inline function slotSwapTeams!(S,i,j,k)
    oppi,oppj = S[i,k],S[j,k]
    if sign(oppi) == sign(oppj)
        S[abs(oppi),k],S[abs(oppj),k] = S[abs(oppj),k],S[abs(oppi),k]
    else
        S[abs(oppi),k],S[abs(oppj),k] = -S[abs(oppj),k],-S[abs(oppi),k]
    end
    S[i,k], S[j,k] = S[j,k], S[i,k]
end

function swapTeams!(S,i,j)
    for col=1:size(S,2)
        if S[i,col] == j || S[j,col] == i
            continue
        end
        slotSwapTeams!(S,i,j,col)
    end
end



function partialSwapRounds!(S,i,k,l)
    swapped = 2^i
    team1, team2 = abs(S[i,k]), abs(S[i,l])
    S[i,k], S[i,l] = S[i,l], S[i,k]
    while true
        cont = false
        if swapped & (1<<team1) === 0
            S[team1,k], S[team1,l] = S[team1,l], S[team1,k]
            swapped += 2^team1
            team1 = abs(S[team1,k])
            cont = true
        end
        if swapped & (1<<team2) === 0
            S[team2,k], S[team2,l] = S[team2,l], S[team2,k]
            swapped += 2^team2
            team2 = abs(S[team2,l])
            cont = true
        end
        cont || break
    end
end

function partialSwapTeams!(S,i,j,k)
    swapped = 2^k
    slotSwapTeams!(S,i,j,k)
    match1, match2 = S[i,k], S[j,k]
    offset = 0
    fixedi,fixedj = false, false
    for col in Iterators.cycle(vcat(k+1:size(S,2), 1:k))
        if !fixedi && S[i,col] == match1
            if swapped & (1<<col) !== 0
                fixedi = true
            else
                slotSwapTeams!(S,i,j,col)
                swapped += 2^col
                match1 = S[i,col]
            end
        end
        if !fixedj && S[j,col] == match2
            if swapped & (1<<col) !== 0
                fixedj = true
            else
                slotSwapTeams!(S,i,j,col)
                swapped += 2^col
                match2 = S[j,col]
            end
        end

        fixedi && fixedj && break
    end
end


function isSchedule(S)
    nteams = size(S,1)
    nslots = 2*nteams-2
    for s=1:nslots
        if count(x->sign(x) == 1, S[:,s]) != nteams/2
            return false
        end
        if sort(abs.(S[:,s])) != 1:nteams
            return false
        end
    end
    for r=1:nteams
        if sort(S[r,:]) != vcat(-nteams:-r-1, -r+1: -1, 1:r-1, r+1:nteams)
            return false
        end
    end
    return true
end

#=
evaluate constraint on given schedule S

returns tuple of hard and soft penalty induced by constraint
=#
function evalConstraint(S, c::CAcon)
    count = 0
    for i in c.teams1
        for s in c.slots
            match = S[i,s]
            if c.mode === MODE_HA || sign(match) == c.mode
                if abs(match) in c.teams2
                    count+=1
                end
            end
        end
    end
    if c.pen == 0
        return max(count-c.ub, 0), 0
    else
        return 0, max(count-c.ub, 0)*c.pen
    end
end

function evalConstraint(S, c::BRcon)
    count = 0
    for i in c.teams
        for s in c.slots
            s==1 && continue
            if c.mode === MODE_HA || sign(S[i,s]) == c.mode
                if sign(S[i,s]) == sign(S[i,s-1])
                    count += 1
                end
            end
        end
    end
    if c.pen == 0
        return max(count-c.ub, 0), 0
    else
        return 0, max(count-c.ub, 0)*c.pen
    end
end

function evalConstraint(S, c::GAcon)
    count = 0
    for s in c.slots
        for i = 1:size(S,1)
            if sign(S[i,s]) == 1 && (i, S[i,s]) in c.meetings
                count += 1
            end
        end
    end
    dif = max(max(c.lb-count, 0), max(count-c.ub, 0))
    if c.pen == 0
        return dif, 0
    else
        return 0, c.pen*dif
    end
end

function evalConstraint(S, c::FAcon)
    return 0,0
end

function evalConstraint(S, c::SEcon)
    count = 0
    games = Array{Int64}(undef, size(S,1), size(S,1))
    for s = 1:size(S,2)
        for i = 1:size(S,1)
            if sign(S[i,s]) == 1
                games[i, S[i,s]] = s
            end
        end
    end
    for i in c.teams
        for j in c.teams
            j <= i && continue
            if abs(games[i,j] - games[j,i]) < c.lb
                count += c.lb - abs(games[i,j] - games[j,i])
            end
        end
    end
    if c.pen == 0
        return count, 0
    else
        return 0, count*c.pen
    end
end
