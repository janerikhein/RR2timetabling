
global const MODE_H = Int8(1)
global const MODE_A = Int8(-1) #dont change values, sign in importnat
global const MODE_HA = Int8(0)

struct IdxSet
    binrepr::Int64
end

function IdxSet(A::Vector{Int})
    binrepr = 0
    for i in A
        if binrepr & (1<<(Int64(i)-1)) !== Int64(0)
            throw(ArgumentError("Invalid Array: has duplicates"))
        end
        if i < 1 || i>64
            throw(ArgumentError("Invalid Array: has values not in range 1:64"))
        end
        binrepr += 2^(i-1)
    end
    return IdxSet(binrepr)
end


Base.eltype(::IdxSet) = Int64

@inline Base.in(n::Int64, s::IdxSet) = s.binrepr & (1<<(n-1)) !== Int64(0)
@inline Base.length(s::IdxSet) = count_ones(s.binrepr)


function Base.iterate(s::IdxSet, state = 1)
    trunc = s.binrepr >>> (state-1)
    trunc === 0 && return nothing
    c = trailing_zeros(trunc)
    state + c, state + c + 1
end

show(io::IO, x::IdxSet) = join(string.(i for i in x), ' ')


struct CAcon
    teams1::IdxSet
    teams2::IdxSet
    slots::IdxSet
    ub::Int8
    mode::Int8
    pen::Int8
end

struct BRcon
    teams::IdxSet
    slots::IdxSet
    ub::Int8
    mode::Int8
    pen::Int8
end

struct GAcon
    slots::IdxSet
    lb::Int8
    ub::Int8
    meetings::Vector{Tuple{Int64,Int64}}
    pen::Int8
end

struct FAcon
    teams::IdxSet
    slots::IdxSet
    ub::Int8
    pen::Int8
end

struct SEcon
    teams::IdxSet
    lb::Int8
    pen::Int8
end

struct RR2
    nteams::UInt8
    nslots::UInt8
    isphased::Bool
end
