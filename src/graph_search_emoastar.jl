"""
Module containing the problem definition structure and the function to find tours within the hopping graph
"""

"""
    MOHGProb

Structure to define a Multi-Objective hopping graph search problem with departure and target conditions.
"""
struct MOHGProb
    dep_cond::SVector{2,Float64}  # departure condition [level_id, alpha]
    tar_cond::SVector{2,Float64}  # target condition [level_id, alpha]
    des_space :: DesignSpace      # design space

    flag_in2out::Bool            # allow in-to-out moon changes
    pump_down_only::Bool         # only pump down flybys (increasing alpha)
end



# """
#     MOHGSearch(dep_cond::SVector{2, Float64}, tar_cond::SVector{2, Float64}, des_space :: DesignSpace, flag_in2out::Bool, pump_down_only::Bool)
    
# Create a MOHGSearch structure defining the hopping graph search problem.

# # Returns
# - `MOHGSearch`: search problem definition
# """
# function MOHGProb(dep_cond::SVector{2, Float64}, tar_cond::SVector{2, Float64}, des_space :: DesignSpace, flag_in2out::Bool, pump_down_only::Bool)
    
#     return MOHGProb(dep_cond, tar_cond, des_space, flag_in2out, pump_down_only)
# end

"""
   MOHGPath
   
Struct to contain the information of a shortest path solution
"""
struct MOHGPath
    ΔV_tot :: Float64
    ToF_tot :: Float64
    path_edges :: Vector{SimpleWeightedEdge}
    seq :: Vector{KepTransfer}
    seq_name :: Vector{String}
end

# """
#     MOHGPath(ΔV_tot::Float64, ToF_tot::Float64, path_edges::Vector{Edge} , seq::Vector{KepTransfer}, seq_name::Vector{String})

# Create a MOHGPath structure containing the information of a shortest path solution.

# # Arguments
# - `ΔV_tot::Float64`: total ΔV of the path [km/s]
# - `ToF_tot::Float64`: total time of flight of the path [days or appropriate unit]
# - `seq::Vector{KepTransfer}`: sequence of KepTransfer objects along the path
# - `seq_name::Vector{String}`: sequence of transfer names along the path

# # Returns
# - `MOHGPath`: path solution
# """
# function MOHGPath(ΔV_tot::Float64, ToF_tot::Float64, path_edges::Vector{SimpleWeightedEdge} ,seq::Vector{KepTransfer}, seq_name::Vector{String})
#     return MOHGPath(ΔV_tot, ToF_tot, path_edges, seq, seq_name);
# end

"""
    prepare_MOHG_search(prob::WeightedHGProb, path::String)

Creates the .gr files for the given problem
# Arguments
- `prob::WeightedHGProb`: the hopping graph search problem to solve
- `prob::String`: the path where to save the .gr files

# Returns
- `Hop_Graph`: the hopping graph of the problem
"""
function prepare_MOHG_search(prob::MOHGProb, path::String)
    # Unpack problem parameters
    dep_cond = prob.dep_cond
    tar_cond = prob.tar_cond
    des_space = prob.des_space
    flag_in2out = prob.flag_in2out
    pump_down_only = prob.pump_down_only
    
    # create hopping graph
    Hop_Graph, src_array, dst_array, wgt_array, dv_tof_array = build_hopping_graph(
        dep_cond, 
        tar_cond, 
        0.0, 
        flag_in2out, 
        pump_down_only, 
        des_space
    );
    
    # Extract ΔV and ToF arrays from dv_tof_array
    ΔV_array = [dvtof[1] for dvtof in dv_tof_array]
    ToF_array = [dvtof[2] for dvtof in dv_tof_array]

    # build .gr files
    dv_filepath, tof_filepath = wrap_edgedef_ΔV_ToF(path, src_array, dst_array, ΔV_array, ToF_array)
    return Hop_Graph

end


"""
    wrap_edgedef_ΔV_ToF(des_space_path::String, src_array::Vector{Int64}, dst_array::Vector{Int64}, ΔV_array::Vector{Float64}, ToF_array::Vector{Float64})

Given the edge definition arrays (src ids, dst ids and ΔV and ToF costs) it creates 2 .gr files in DIMACS format, ready to be fed to EMOA* c++ implementation.

# Arguments
- `des_space_path::String`: Path to the design space folder where .gr files will be saved
- `src_array::Vector{Int64}`: Source node IDs for each edge
- `dst_array::Vector{Int64}`: Destination node IDs for each edge
- `ΔV_array::Vector{Float64}`: ΔV costs for each edge [m/s]
- `ToF_array::Vector{Float64}`: Time of flight for each edge [days]

# Output Files
- `graph_dv.gr`: Graph file with ΔV as edge weights
- `graph_tof.gr`: Graph file with ToF as edge weights

# Format
The files follow the DIMACS shortest path format:
- Problem line: `p sp <nodes> <edges>`
- Arc lines: `a <src> <dst> <weight>`
"""
function wrap_edgedef_ΔV_ToF(des_space_path::String, src_array::Vector{Int64}, dst_array::Vector{Int64}, ΔV_array::Vector{Float64}, ToF_array::Vector{Float64})
    
    # Validate input arrays have same length
    n_edges = length(src_array)
    if length(dst_array) != n_edges || length(ΔV_array) != n_edges || length(ToF_array) != n_edges
        error("All input arrays must have the same length")
    end
    
    # Determine number of nodes (maximum node ID in src and dst arrays)
    n_nodes = max(maximum(src_array), maximum(dst_array))
    
    println("Creating DIMACS graph files with $n_nodes nodes and $n_edges edges")
    
    # Create output file paths
    dv_filepath = joinpath(des_space_path, "graph_dv.gr")
    tof_filepath = joinpath(des_space_path, "graph_tof.gr")
    
    # Write ΔV graph file
    open(dv_filepath, "w") do io
        # Write problem line
        println(io, "c DIMACS shortest path graph - ΔV objective [m/s]")
        println(io, "p sp $n_nodes $n_edges")
        
        # Write arc lines
        for i in 1:n_edges
            # Format: a <src> <dst> <weight>
            println(io, "a $(src_array[i]) $(dst_array[i]) $(ΔV_array[i])")
        end
    end
    
    println("Created ΔV graph file: $dv_filepath")
    
    # Write ToF graph file
    open(tof_filepath, "w") do io
        # Write problem line
        println(io, "c DIMACS shortest path graph - ToF objective [days]")
        println(io, "p sp $n_nodes $n_edges")
        
        # Write arc lines
        for i in 1:n_edges
            # Format: a <src> <dst> <weight>
            println(io, "a $(src_array[i]) $(dst_array[i]) $(ToF_array[i])")
        end
    end
    
    println("Created ToF graph file: $tof_filepath")
    
    return dv_filepath, tof_filepath
end



"""
    MOHGSol

struct to collect and store the solutions of a MOHGProb problem
"""
struct MOHGSol
    prb :: MOHGProb
    Nsol :: Int64

    paths :: Vector{MOHGPath}
    ΔV_array :: Vector{Float64}
    ToF_array :: Vector{Float64}

end

"""
    parse_EMOA_results(results_path::String, Hop_Graph::HoppingGraphSingle, des_space::DesignSpace, prob::MOHGProb)

Parses the EMOA* results file and creates a MOHGSol structure containing all Pareto-optimal solutions.

# Arguments
- `results_path::String`: Path to the EMOA_results.txt file
- `Hop_Graph::HoppingGraphSingle`: The hopping graph used for the search
- `prob::MOHGProb`: The multi-objective hopping graph problem definition

# Returns
- `MOHGSol`: Solution structure containing all Pareto-optimal paths with their ΔV and ToF costs

# File Format
The EMOA_results.txt file has the following structure:
- Header lines with search statistics
- For each solution (starting from line with "Label:"):
  - Label line: "Label: <id>"
  - Cost line: "[<ΔV>,<ToF>,]"
  - Path line: sequence of node IDs separated by spaces
"""
function parse_EMOA_results(results_path::String, Hop_Graph::HoppingGraphSingle, prob::MOHGProb)
    
    # Read the entire file
    lines = readlines(results_path)
    
    # Find where the solutions start (after the header)
    solution_start_idx = findfirst(contains("Label:"), lines)
    
    if isnothing(solution_start_idx)
        error("No solutions found in EMOA results file")
    end
    
    # Parse solutions
    paths = MOHGPath[]
    ΔV_array = Float64[]
    ToF_array = Float64[]
    
    i = solution_start_idx
    while i <= length(lines)
        line = lines[i]
        
        # Check if this is a label line (start of a new solution)
        if startswith(line, "Label:")
            # Next line contains costs [ΔV, ToF,]
            i += 1
            cost_line = strip(lines[i])
            
            # Parse costs - format: [476.119021,694.833418,]
            cost_line = replace(cost_line, "[" => "")
            cost_line = replace(cost_line, "]" => "")
            # Split on comma FIRST, before removing commas
            cost_strings = strip.(split(cost_line, ","))
            # Filter out empty strings
            cost_strings = filter(!isempty, cost_strings)
            costs = parse.(Float64, cost_strings)
            
            ΔV_tot = costs[1]
            ToF_tot = costs[2]
            
            # Next line contains the path (node IDs)
            i += 1
            path_line = strip(lines[i])
            node_ids = parse.(Int, split(path_line))
            
            # Convert node sequence to MOHGPath
            # Skip first node (virtual departure) and last node (virtual target)
            path_nodes = node_ids[2:end-1]
            
            # Build the path
            path = build_MOHGPath_from_nodes(path_nodes, ΔV_tot, ToF_tot, Hop_Graph)
            
            # Store solution
            push!(paths, path)
            push!(ΔV_array, ΔV_tot)
            push!(ToF_array, ToF_tot)
        end
        
        i += 1
    end
    
    Nsol = length(paths)
    
    println("Parsed $Nsol Pareto-optimal solutions from EMOA* results")
    
    # Sort solutions by increasing ΔV
    sort_idx = sortperm(ΔV_array)
    paths = paths[sort_idx]
    ΔV_array = ΔV_array[sort_idx]
    ToF_array = ToF_array[sort_idx]

    return MOHGSol(prob, Nsol, paths, ΔV_array, ToF_array)
end

"""
    build_MOHGPath_from_nodes(node_ids::Vector{Int}, ΔV_tot::Float64, ToF_tot::Float64, Hop_Graph::HoppingGraphSingle)

Constructs a MOHGPath from a sequence of node IDs.

# Arguments
- `node_ids::Vector{Int}`: Sequence of node IDs in the path
- `ΔV_tot::Float64`: Total ΔV cost of the path
- `ToF_tot::Float64`: Total time of flight of the path
- `Hop_Graph::HoppingGraphSingle`: Hopping graph containing node-to-transfer mappings

# Returns
- `MOHGPath`: Path object with transfer sequence and costs
"""
function build_MOHGPath_from_nodes(node_ids::Vector{Int}, ΔV_tot::Float64, ToF_tot::Float64, Hop_Graph::HoppingGraphSingle)
    
    fb_to_node = Hop_Graph.fb_to_node
    
    seq = KepTransfer[]
    seq_name = String[]
    path_edges = SimpleWeightedEdge[]
    
    # Process each node in the sequence
    for (i, node_id) in enumerate(node_ids)
        
        # Get transfers associated with this node
        if haskey(fb_to_node, node_id)
            node_transfers = fb_to_node[node_id]
            entry_transfer = node_transfers[1]
            exit_transfer = node_transfers[2]
            
            # For first node, add the entry transfer
            if i == 1
                push!(seq, entry_transfer)
                push!(seq_name, entry_transfer.name)
            end
            
            # Add the exit transfer (this is the transfer leaving this node)
            push!(seq, exit_transfer)
            push!(seq_name, exit_transfer.name)
            
            # Create edge if not the last node
            if i < length(node_ids)
                next_node_id = node_ids[i + 1]
                edge = SimpleWeightedEdge(node_id, next_node_id, 0.0)
                push!(path_edges, edge)
            end
        else
            @warn "Node $node_id not found in fb_to_node mapping"
        end
    end
    
    return MOHGPath(ΔV_tot, ToF_tot, path_edges, seq, seq_name)
end

