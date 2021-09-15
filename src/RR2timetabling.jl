module RR2timetabling

using JuMP, Gurobi
using LinearAlgebra
using Combinatorics
using Random
using StructArrays
using ImportMacros
using Setfield
@import LightXML as XML


include("structs.jl")
include("util.jl")
include("xmlio.jl")
include("ipmodels.jl")
include("constraints.jl")
include("decomposition.jl")
include("localsearch.jl")
include("neighborhood.jl")
include("heuristic.jl")


export parse_robinx, toxml
export ip_model, base_ip_model
export heur_solution, pgba_heuristic

end # module
