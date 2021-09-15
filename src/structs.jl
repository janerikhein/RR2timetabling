global const MODE_H = Int8(1)
global const MODE_A = Int8(-1)
global const MODE_HA = Int8(0)

global const GRB_ENV = Gurobi.Env()
global const GRB_ENV_2 = Gurobi.Env()
global const GRB_ENV_3 = Gurobi.Env()

struct IdxSet
    binrepr::Int64
end

struct RR2
    nteams::Int
    nslots::Int
    isphased::Bool
end

abstract type Constraint end

struct CAcon <: Constraint
    teams1::IdxSet
    teams2::IdxSet
    slots::IdxSet
    ub::Int
    mode::Int
    pen::Int
end

struct BRcon <: Constraint
    teams::IdxSet
    slots::IdxSet
    ub::Int
    mode::Int
    pen::Int
end

struct GAcon <: Constraint
    slots::IdxSet
    lb::Int
    ub::Int
    meetings::Vector{Tuple{Int,Int}}
    pen::Int
end

struct FAcon <: Constraint
    teams::IdxSet
    slots::IdxSet
    ub::Int
    pen::Int
end

struct SEcon <: Constraint
    teams::IdxSet
    lb::Int
    pen::Int
end

struct Instance
    rr2::RR2
    constrCA::Vector{CAcon}
    constrBR::Vector{BRcon}
    constrGA::Vector{GAcon}
    constrFA::Vector{FAcon}
    constrSE::Vector{SEcon}
end

struct TabuSearchOptimizer
    maxiter::Int
    tabulength::Int
    rr2::RR2
    constraints::Vector{Constraint}
    conentries::BitArray{3}
end

const Move = Tuple{Symbol, Int64, Int64, Vararg{Int64, N} where N}

mutable struct SearchState
    it::Int # iteration count
    val::Float64 # current value
    sched::Matrix{Int}
    optval::Float64 # best known value
    optsched::Matrix{Int} # best known feasible schedule
    lastimprov::Int # last iteration in which global improvement of obj function was fou
    offset::Vector{Int}
    moves::Vector{Move}
    tabu::Dict{Move,Int}
    blocked::Dict{Move,Int}
    pen::Vector{Float64} # penalty for each constraint
    fr::Matrix{Int}
    diversify::Bool
end

function IdxSet(A::Vector{T}) where {T<:Integer}
    binrepr = 0
    for i in A
        if binrepr & (1 << (i - 1)) !== 0
            throw(ArgumentError("Invalid Array: has duplicates"))
        end
        if i < 1 || i > 64
            throw(ArgumentError("Invalid Array: has values not in range 1:64"))
        end
        binrepr += 2^(i - 1)
    end
    return IdxSet(binrepr)
end

Base.eltype(::IdxSet) = Int
@inline Base.in(n::Int64, s::IdxSet) = s.binrepr & (1 << (n - 1)) !== 0
@inline Base.length(s::IdxSet) = count_ones(s.binrepr)

function Base.iterate(s::IdxSet, state = 1)
    trunc = s.binrepr >>> (state - 1)
    trunc === 0 && return nothing
    c = trailing_zeros(trunc)
    return state + c, state + c + 1
end

show(io::IO, x::IdxSet) = join(string.(i for i in x), ' ')

function Instance(rr2::RR2, cons::Vector{Constraint})
    D = Dict(CAcon => 1, BRcon => 2, GAcon => 3, FAcon => 4, SEcon => 5)
    groupedcons = (
        Vector{CAcon}(),
        Vector{BRcon}(),
        Vector{GAcon}(),
        Vector{FAcon}(),
        Vector{SEcon}(),
    )
    for c in cons
        push!(groupedcons[D[typeof(c)]], c)
    end
    return Instance(rr2, groupedcons...)
end

struct PGBACut
    ub::Int
    assignments::Vector{Tuple{Int,Int}}
end
