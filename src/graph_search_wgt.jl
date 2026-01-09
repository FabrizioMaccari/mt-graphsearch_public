"""
Module containing the problem definition structure and the function to find tours within the hopping graph
"""

"""
    WeightedHGProb

Structure to define a weighted hopping graph search problem with departure and target conditions.
"""
struct WeightedHGProb
    dep_cond::SVector{2,Float64}  # departure condition [level_id, alpha]
    tar_cond::SVector{2,Float64}  # target condition [level_id, alpha]
    des_space :: DesignSpace      # design space

    flag_in2out::Bool            # allow in-to-out moon changes
    pump_down_only::Bool         # only pump down flybys (increasing alpha)
    
    p_ΔV_array:: Vector{Float64}  # cost weights: 0=time optimal, 1=fuel optimal
end



"""
    WeightedHGProb(dep_cond::SVector{2, Float64}, tar_cond::SVector{2, Float64}, des_space :: DesignSpace, p_ΔV_array::Vector{Float64}, flag_in2out::Bool, pump_down_only::Bool)
    
Create a WeightedHGProb structure defining the hopping graph search problem.

# Returns
- `WeightedHGProb`: search problem definition
"""
function WeightedHGProb(dep_cond::SVector{2, Float64}, tar_cond::SVector{2, Float64}, des_space :: DesignSpace, p_ΔV_array::Vector{Float64}, flag_in2out::Bool, pump_down_only::Bool)
    
    return WeightedHGProb(dep_cond, tar_cond, des_space, flag_in2out, pump_down_only ,p_ΔV_array)
end

"""
   WeightedHGPath
   
Struct to contain the information of a shortest path solution
"""
struct WeightedHGPath
    ΔV_tot :: Float64
    ToF_tot :: Float64
    path_edges :: Vector{SimpleWeightedEdge}
    seq :: Vector{KepTransfer}
    seq_name :: Vector{String}
    p_ΔV :: Float64 
end

"""
    WeightedHGPath(ΔV_tot::Float64, ToF_tot::Float64, path_edges::Vector{Edge} , seq::Vector{KepTransfer}, seq_name::Vector{String}, p_ΔV::Float64)

Create a WeightedHGPath structure containing the information of a shortest path solution.

# Arguments
- `ΔV_tot::Float64`: total ΔV of the path [km/s]
- `ToF_tot::Float64`: total time of flight of the path [days or appropriate unit]
- `seq::Vector{KepTransfer}`: sequence of KepTransfer objects along the path
- `seq_name::Vector{String}`: sequence of transfer names along the path
- `p_ΔV::Float64`: cost weight used for this path (0=time optimal, 1=fuel optimal)

# Returns
- `WeightedHGPath`: path solution
"""
function WeightedHGPath(ΔV_tot::Float64, ToF_tot::Float64, path_edges::Vector{SimpleWeightedEdge} ,seq::Vector{KepTransfer}, seq_name::Vector{String}, p_ΔV::Float64)
    return WeightedHGPath(ΔV_tot, ToF_tot, path_edges, seq, seq_name, p_ΔV);
end

"""
    WeightedHGPath(sp_weighted::Vector, fb_to_node::Dict, id_to_transfer::Dict, p_ΔV::Float64)

Create a WeightedHGPath from a shortest path solution (e.g., from A* algorithm).

# Arguments
- `sp_weighted::Vector`: shortest path as vector of edges
- `fb_to_node::Dict`: mapping from node ID to tuple of (entry_transfer, exit_transfer)
- `id_to_transfer::Dict{Int, KepTransfer}`: mapping from transfer ID to KepTransfer object
- `p_ΔV::Float64`: cost weight used (0=time optimal, 1=fuel optimal)

# Returns
- `WeightedHGPath`: constructed path solution
"""
function WeightedHGPath(sp_weighted::Vector, fb_to_node::Dict, id_to_transfer::Dict{Int, KepTransfer}, p_ΔV::Float64)
    ΔV_tot = 0.0
    ToF_tot = 0.0
    seq = KepTransfer[]
    seq_name = String[]
    
    # Skip first edge (from virtual departure node)
    for (i,ed) in enumerate(sp_weighted[2:end])
        sr = src(ed)
        ds = dst(ed)
        
        # # Skip last edge (to virtual target node)
        # if ds == 2
        #     continue
        # end

        # store first transfer
        if i == 1
            first_transfer = fb_to_node[sr][1]
            
            # add to sequence
            push!(seq, first_transfer)
            push!(seq_name, first_transfer.name)

            # consider costs
            ΔV_tot  += first_transfer.ΔV
            ToF_tot += first_transfer.ToF
        end
        

        # Get the exit transfer from the source node
        node_transfers = fb_to_node[sr]
        exit_transfer = node_transfers[2]
        
        # Add to sequence
        push!(seq, exit_transfer)
        push!(seq_name, exit_transfer.name)
        
        # Accumulate costs
        ΔV_tot += exit_transfer.ΔV
        ToF_tot += exit_transfer.ToF
    end
    
    return WeightedHGPath(ΔV_tot, ToF_tot, sp_weighted, seq, seq_name, p_ΔV)
end



"""
    WeightedHGSol

struct to solve and store the solutions of a WeightedHGProb problem
"""
struct WeightedHGSol
    prb :: WeightedHGProb
    Nsol :: Float64

    paths :: Vector{WeightedHGPath}
    ΔV_array :: Vector{Float64}
    ToF_array :: Vector{Float64}

end

"""
    WeightedHGSol(prob::WeightedHGProb)

Solves the given problem, storing the array of shortest path solutions for each p_ΔV value.

# Arguments
- `prob::WeightedHGProb`: the hopping graph search problem to solve

# Returns
- `WeightedHGSol`: solution structure containing all paths and their costs
"""
function WeightedHGSol(prob::WeightedHGProb)
    # Unpack problem parameters
    dep_cond = prob.dep_cond
    tar_cond = prob.tar_cond
    des_space = prob.des_space
    p_ΔV_array = prob.p_ΔV_array
    flag_in2out = prob.flag_in2out
    pump_down_only = prob.pump_down_only
    
    # Initialize solution storage
    Nsol = length(p_ΔV_array)
    paths = Vector{WeightedHGPath}(undef, Nsol)
    ΔV_array = zeros(Float64, Nsol)
    ToF_array = zeros(Float64, Nsol)

    # for the first iteration the nodes and edges must be defined
    Hop_Graph, src_array, dst_array, wgt_array, dv_tof_array = build_hopping_graph(
        dep_cond, 
        tar_cond, 
        p_ΔV_array[1], 
        flag_in2out, 
        pump_down_only, 
        des_space
    );
    graph = Hop_Graph.graph;
    weights_gs = weights(graph)
    
    # Solve for each p_ΔV value
    for (i, p_ΔV) in enumerate(p_ΔV_array)
        println("Solving for p_ΔV = $p_ΔV ($i of $Nsol)")
        
        # Build hopping graph for this cost weight
        if i > 1
            # from second run onwards only the edge weights need to be updated
            wgt_array = update_wgts(p_ΔV, dv_tof_array)
            graph = SimpleWeightedDiGraph(src_array, dst_array, wgt_array);
            weights_gs = weights(graph)
        end
        
        # Find shortest path using A* with no heuristic - Dijkstra
        sp_result = a_star(graph, 1, 2, weights_gs)
        
        # Extract path
        path = WeightedHGPath(sp_result, Hop_Graph.fb_to_node, des_space.id_to_transfer, p_ΔV)
        
        # Store results
        paths[i] = path
        ΔV_array[i] = path.ΔV_tot
        ToF_array[i] = path.ToF_tot
        
        if length(sp_result) == 0
            println("  search failed for p_ΔV =$(p_ΔV)")
        else
            println("  ΔV = $(path.ΔV_tot) m/s, ToF = $(path.ToF_tot) days")
        end
    end
    
    return WeightedHGSol(prob, Nsol, paths, ΔV_array, ToF_array);
end