"""
    PDPCR

Heuristic cho bài toán Pickup and Delivery with Cooperative Robots (PDP-CR):
tabu search trên đồ thị luồng với các lân cận reroute, reverse và redistribute
(Nguyen, Quilliot, Gayon, Nguyen, Ta — "Redistribute neighborhood in the PDP-CR").
"""
module PDPCR

using Random
using Printf

include("instance.jl")
include("solution.jl")
include("initial.jl")
include("neighborhoods.jl")
include("tabu.jl")

export Instance, example1, read_lilim, random_instance, write_instance, read_instance,
       source, sink, nnodes,
       Solution, add_flow!, flow, arcs, Workspace, makespan!, tails!, critical_arcs,
       topo_by_start, check_feasible, robot_routes, print_schedule,
       greedy_solution,
       Move, reroute_moves, reverse_moves, redistribute, apply!, undo!,
       TabuParams, tabu_search, TabuResult

end
