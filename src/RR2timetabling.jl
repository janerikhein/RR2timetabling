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

export optimize!
export Instance, RR2
export parse_robinx, toxml, tofile, ishard
export ip_model, base_ip_model
export parse_solution, heur_solution
export heur_solution, pgba_heuristic
export eval_constraint, constraints
export is_pattern_constraint, is_pattern_set_constraint


end # module
