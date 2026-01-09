"""
Module containing the functions to build the hopping graph, taking as input the design space
"""

"""
    HoppingGraphSingle
Single objective hopping graph
"""
struct HoppingGraphSingle
    id :: String
    graph :: SimpleWeightedDiGraph
    fb_to_node :: Dict{Int, Vector{Any}}

    dep_cond::SVector{2, Float64}  # level_id, α
    tar_cond::SVector{2, Float64}  # level_id, α - if α is -1 all the tours getting to levels[level_id] are considered
    p_ΔV :: Float64

end

"""
    build_hopping_graph(dep_cond::SVector{2, Float64}, tar_cond::SVector{2, Float64}, p_ΔV :: Float64, pump_down_only::Bool, des_space::DesignSpace)
Given the input design space, it assembles the hopping graph, ready to go through the optimization
If flag_in2out is true, it includes also moon changes going from an inner to an outer moon
If pump_down_only is true, it includes only pump down flybys
"""
function build_hopping_graph(dep_cond::SVector{2,Float64}, tar_cond::SVector{2, Float64}, p_ΔV :: Float64, flag_in2out::Bool, pump_down_only::Bool, des_space::DesignSpace)
    id = des_space.id  # recover id to link the design space with the corresponding graph

    levels = des_space.vinf_levels  # recover vinf levels
    println(length(levels))
    # intialize edge definition list
    src_array = Vector{Int}()
    dst_array = Vector{Int}()
    wgt_array = Vector{Float64}()
    dvtof_list = Vector{Vector{Float64}}()  # Vector of [ΔV, ToF] pairs

    # create source and end node
    fb_to_node = Dict{Int, Vector{Any}}()
    fb_to_node[1] = dep_cond
    fb_to_node[2] = tar_cond
    curr_id = 3


    # cycle through the levels and compile the edge definition lists
    src_levels, dst_levels, wgt_levels, dvtof_list_levels, fb_to_node_levels, curr_id = build_level_edges(curr_id, p_ΔV, levels, pump_down_only)
    append!(src_array, src_levels)
    append!(dst_array, dst_levels)
    append!(wgt_array, wgt_levels)
    append!(dvtof_list, dvtof_list_levels)
    merge!(fb_to_node, fb_to_node_levels)

    println("Level nodes and edges computed: "*string(curr_id-2)*" nodes and "*string(length(src_array))*" edges")
    
    # compute edges for moon changes
    # create moon change nodes
    fb_to_node_mc, curr_id = build_moon_change_nodes(levels, curr_id, flag_in2out)
    merge!(fb_to_node, fb_to_node_mc)
    
    # create edges with the nodes within each level
    
    # Precompute lookup tables for transfer objects to node keys
    tr2_to_keys = Dict{Any, Vector{Int}}()
    tr1_to_keys = Dict{Any, Vector{Int}}()
    for (k, v) in fb_to_node
        push!(get!(tr2_to_keys, v[2], Int[]), k)
        push!(get!(tr1_to_keys, v[1], Int[]), k)
    end
    
    
    mc_node_keys = collect(keys(fb_to_node_mc))
    println("creating edges for "*string(length(mc_node_keys))*" moon change nodes")
    for (i,mc_key) in enumerate(mc_node_keys)
        # recover node
        mc_node = fb_to_node_mc[mc_key]
        tr1 = mc_node[1]
        tr2 = mc_node[2]

        # recover nodes leading to the moon change
        nodes_to_mc = get(tr2_to_keys, tr1, Int[])

        # add the corresponding edges
        append!(src_array, nodes_to_mc)
        append!(dst_array, mc_key*ones(length(nodes_to_mc)))
        append!(wgt_array, (p_ΔV*tr1.ΔV + tr1.ToF*(1 - p_ΔV))*ones(length(nodes_to_mc)))
        append!(dvtof_list, [[tr1.ΔV, tr1.ToF] for _ in 1:length(nodes_to_mc)])

        # recover nodes that can be visited after the moon change
        mc_to_nodes = get(tr1_to_keys, tr2, Int[])
        # add the corrsponding edges
        append!(src_array, mc_key*ones(length(mc_to_nodes)))
        append!(dst_array, mc_to_nodes)
        append!(wgt_array, (p_ΔV*tr2.ΔV + tr2.ToF*(1 - p_ΔV))*ones(length(mc_to_nodes)))
        append!(dvtof_list, [[tr2.ΔV, tr2.ToF] for _ in 1:length(mc_to_nodes)])
    end
    println("Added moon change nodes and edges, total graph size: "*string(curr_id-2)*" nodes and "*string(length(src_array))*" edges")
    
    # intialize specific problem
    # identify nodes connected to departure condition
    println("initializing problem")
    level_dep_id, α_dep = dep_cond
    level_dep = levels[Int(level_dep_id)]
    δmax_dep = level_dep.δmax
    dep_nodes = [k for (k, v) in fb_to_node if k > 2 && (abs(v[1].arc1.α  - α_dep) <= δmax_dep ) && (v[1].arc1.v∞ == level_dep.v∞)]
    dep_dvs = [v[1].ΔV for (k, v) in fb_to_node if k > 2 && (abs(v[1].arc1.α  - α_dep) <= δmax_dep ) && (v[1].arc1.v∞ == level_dep.v∞)]
    dep_tofs = [v[1].ToF for (k, v) in fb_to_node if k > 2 && (abs(v[1].arc1.α  - α_dep) <= δmax_dep ) && (v[1].arc1.v∞ == level_dep.v∞)]
    
    # create edges to start the tour
    println("number of departure nodes: "*string(length(dep_nodes)))

    append!(src_array, ones(length(dep_nodes)))
    append!(dst_array, dep_nodes)
    # append!(wgt_array, 0.0*zeros(length(dep_nodes)))
    # append!(dvtof_list, [[0.0, 0.0] for _ in 1:length(dep_nodes)])
    append!(wgt_array, p_ΔV.*dep_dvs .+ (1 - p_ΔV).*dep_tofs)
    append!(dvtof_list, [[dep_dvs[i], dep_tofs[i]] for i in 1:length(dep_nodes)])

    # identify nodes connected to target condition
    level_tar_id, α_tar = tar_cond
    level_tar = levels[Int(level_tar_id)]
    δmax_tar = level_dep.δmax
    if α_tar != -1.0
        tar_nodes = [k for (k, v) in fb_to_node if k > 2 &&(abs(v[2].arc2.α  - α_tar) <= δmax_tar ) && (v[2].arc2.v∞ == level_tar.v∞)]
        tar_dvs = [v[2].ΔV for (k, v) in fb_to_node if k > 2 && (abs(v[2].arc2.α  - α_tar) <= δmax_tar ) && (v[2].arc2.v∞ == level_tar.v∞)]
        tar_tofs = [v[2].ToF for (k, v) in fb_to_node if k > 2 && (abs(v[2].arc2.α  - α_tar) <= δmax_tar ) && (v[2].arc2.v∞ == level_tar.v∞)]
    else
        tar_nodes = [k for (k, v) in fb_to_node if k > 2 && (v[2].arc2.v∞ == level_tar.v∞)]
        tar_dvs = [v[2].ΔV for (k, v) in fb_to_node if k > 2  && (v[2].arc2.v∞ == level_tar.v∞)]
        tar_tofs = [v[2].ToF for (k, v) in fb_to_node if k > 2 && (v[2].arc2.v∞ == level_tar.v∞)]
    end
    # create null weight edges to start the tour
    println("number of target nodes: "*string(length(tar_nodes)))
    append!(src_array, tar_nodes)
    append!(dst_array, 2*ones(length(tar_nodes)))
    # append!(wgt_array, 0.0*zeros(length(tar_nodes)))
    # append!(dvtof_list, [[0.0, 0.0] for _ in 1:length(tar_nodes)])
    append!(wgt_array, p_ΔV.*tar_dvs .+ (1 - p_ΔV).*tar_tofs)
    append!(dvtof_list, [[tar_dvs[i], tar_tofs[i]] for i in 1:length(tar_nodes)])


    # build graph and provide metrics
    println("creating graph")
    
    graph = SimpleWeightedDiGraph(src_array, dst_array, wgt_array);  
    println("Hopping Graph built, $(curr_id - 2) vertices and $(length(src_array)) edges")
    
    return HoppingGraphSingle(id, graph, fb_to_node, dep_cond, tar_cond, p_ΔV), src_array, dst_array, wgt_array, dvtof_list

end

"""
    build_level_edges(curr_id::Int,  p_ΔV::Float64, flag_in2out::Bool, levels::Vector{VinfLevel},pump_down_only::Bool=true)
Given a vinf level list, it produces the hopping graph definition variables (src,dst, wgt) and the translation
dictionary. The nodes are named with id starting from curr_id +1.
"""
function build_level_edges(curr_id::Int, p_ΔV::Float64, levels::Vector{VinfLevel}, pump_down_only::Bool=true)
    # intialize edge definition list
    src_array = Vector{Int}()
    dst_array = Vector{Int}()
    wgt_array = Vector{Float64}()
    dvtof_list = Vector{Vector{Float64}}()  # Vector of [ΔV, ToF] pairs
    fb_to_node = Dict{Int, Vector{Any}}()
    n_levels = length(levels)
    for level_id in eachindex(levels)
        level = levels[level_id]
        # assemble resonance - resonance nodes
        δmax = level.δmax
        println(rad2deg(δmax))
        fb_to_node_comb, curr_id = build_nodes_from_transfers(curr_id, δmax, level.res_transfers, p_ΔV,pump_down_only)
        merge!(fb_to_node, fb_to_node_comb)
     
        # assemble resonance - non-resonance nodes
        fb_to_node_comb, curr_id = build_nodes_from_transfers(curr_id, δmax, level.res_transfers, level.nonres_transfers_I, p_ΔV,pump_down_only)
        merge!(fb_to_node, fb_to_node_comb)
        
        fb_to_node_comb, curr_id = build_nodes_from_transfers(curr_id, δmax, level.res_transfers, level.nonres_transfers_O, p_ΔV,pump_down_only)
        merge!(fb_to_node, fb_to_node_comb)

        # assemble resonance - vilt nodes
        fb_to_node_comb, curr_id = build_nodes_from_transfers(curr_id, δmax, level.res_transfers, level.vilts_outgoing, p_ΔV , pump_down_only, false)
        merge!(fb_to_node, fb_to_node_comb)

        fb_to_node_comb, curr_id = build_nodes_from_transfers(curr_id, δmax, level.vilts_incoming, level.res_transfers, p_ΔV , pump_down_only, false)
        merge!(fb_to_node, fb_to_node_comb)
        
        # assemble non-resonant - non-resonant nodes
        fb_to_node_comb, curr_id = build_nodes_from_transfers(curr_id, δmax, level.nonres_transfers_I, level.nonres_transfers_O, p_ΔV,pump_down_only)
        merge!(fb_to_node, fb_to_node_comb)
        
        # assemble non-resonant - vilt nodes
        fb_to_node_comb, curr_id = build_nodes_from_transfers(curr_id, δmax, level.nonres_transfers_O, level.vilts_outgoing, p_ΔV , pump_down_only, false)
        merge!(fb_to_node, fb_to_node_comb)
        
        fb_to_node_comb, curr_id = build_nodes_from_transfers(curr_id, δmax, level.nonres_transfers_I, level.vilts_outgoing, p_ΔV , pump_down_only, false)
        merge!(fb_to_node, fb_to_node_comb)

        fb_to_node_comb, curr_id = build_nodes_from_transfers(curr_id, δmax, level.vilts_incoming, level.nonres_transfers_I, p_ΔV , pump_down_only, false)
        merge!(fb_to_node, fb_to_node_comb)

        fb_to_node_comb, curr_id = build_nodes_from_transfers(curr_id, δmax, level.vilts_incoming, level.nonres_transfers_O, p_ΔV , pump_down_only, false)
        merge!(fb_to_node, fb_to_node_comb)
        
        # assemble vilt - vilt nodes 
        fb_to_node_comb, curr_id = build_nodes_from_transfers(curr_id, δmax, level.vilts_incoming, level.vilts_outgoing, p_ΔV , pump_down_only, false)
        merge!(fb_to_node, fb_to_node_comb)
        
        println("Processed level "*string(level_id)*" of "*string(n_levels))

    end

    # now create edges between the nodes (an edge exists when the exit transfer of a flyby is the entry transfer of another)
    ids = collect(keys(fb_to_node))
    n_keys = length(ids)

    # Precompute lookups - in terms of transfer id
    exit_id_to_nodes = Dict{Int, Vector{Int}}()
    entry_id_to_nodes = Dict{Int, Vector{Int}}()
    exit_id_to_weight = Dict{Int, Float64}()
    exit_id_to_dvtof = Dict{Int, Vector{Float64}}()
    for k in ids
        exit_tr = fb_to_node[k][2]
        exit_id = fb_to_node[k][2].id
        entry_id = fb_to_node[k][1].id

        exit_id_to_weight[exit_id] = exit_tr.ΔV * p_ΔV + exit_tr.ToF * (1 - p_ΔV)
        exit_id_to_dvtof[exit_id] = [exit_tr.ΔV, exit_tr.ToF]
        
        push!(get!(exit_id_to_nodes, exit_id, Int[]), k)
        push!(get!(entry_id_to_nodes, entry_id, Int[]), k)
    end

    # Efficient edge generation
    for id in intersect(keys(exit_id_to_nodes), keys(entry_id_to_nodes))
        wgt = exit_id_to_weight[id]
        dvtof = exit_id_to_dvtof[id]  # [ΔV, ToF]

        for src in exit_id_to_nodes[id], dst in entry_id_to_nodes[id]
            # add edge src -> dst
            push!(src_array, src)
            push!(dst_array, dst)
            push!(wgt_array, wgt)
            push!(dvtof_list, dvtof)  # Much faster!
        end
    end
    
    println(string(length(src_array))*" level edges generated")    
    return src_array, dst_array, wgt_array, dvtof_list, fb_to_node, curr_id
end

"""
    build_nodes_from_transfers(curr_id::Int,δ_max::Float64, KepTransfers1::Vector{KepTransfer}, KepTransfers2::Vector{KepTransfer}, pump_down_only::Bool=True, flag_j2i::Bool=true)

Given two lists of Keplerian transfers, it provides the hopping graph nodes and the translation 
dictionary, for the nodes linking with a flyby of maximum turning angle delta_max one element of KepTransfers1 to one of KepTransfers2 and if flag_j2i=true viceversa.The nodes are named with id starting from curr_id +1 
"""
function build_nodes_from_transfers(curr_id::Int, δ_max::Float64, KepTransfers1::Vector{KepTransfer}, KepTransfers2::Vector{KepTransfer},p_ΔV::Float64, pump_down_only::Bool=true, flag_j2i::Bool=true)
    # intialize edge definition list
    fb_to_node = Dict{Int, Vector{Any}}()

    # choose whether to include pump down flybys only
    if pump_down_only
        check_fb = (α1, α2, δ_max) -> check_feas_fb_pumpdown(α1, α2, δ_max)
    else
        check_fb = (α1, α2, δ_max) -> check_feas_fb(α1, α2, δ_max)
    end
    # cycle through the transfers and build nodes
    
    for (i, transf_i) in enumerate(KepTransfers1)
        
        # check i to j flyby
        α_i = transf_i.arc2.α
        k_io_i = transf_i.k_io2
        
        # check j to i fb
        α_i_ji = transf_i.arc1.α
        k_io_i_ji = transf_i.k_io1

        # check from current pump angle onwards
        KepTransfers2_tocheck = KepTransfers2
        for (j, transf_j) in enumerate(KepTransfers2_tocheck)
            # check i to j fb
            α_j = transf_j.arc1.α
            k_io_j = transf_j.k_io1
            if check_fb(α_i, α_j, δ_max) && (k_io_i == k_io_j)
                fb_to_node[curr_id] = [transf_i, transf_j]
                curr_id += 1
            end

            # check j to i fb
            if flag_j2i
                α_j_ji = transf_j.arc2.α
                k_io_j_ji = transf_j.k_io2
                if check_fb(α_j_ji, α_i_ji, δ_max) && (k_io_i_ji == k_io_j_ji)
                    fb_to_node[curr_id] = [transf_j, transf_i]
                    curr_id += 1
                end
            end
        end
    end
    
    return fb_to_node, curr_id
end

"""
    build_nodes_from_transfers(curr_id::Int,δ_max::Float64, KepTransfers::Vector{KepTransfers}, p_ΔV::Float64,pump_down_only::Bool=True, flag_j2i=true)
    
Given a list of Keplerian transfers, it provides the hopping graph nodes and the translation 
dictionary, for the nodes linking with a flyby of maximum turning angle delta_max one element of KepTransfers to another.The nodes are named with id starting from curr_id +1  
"""
function build_nodes_from_transfers(curr_id::Int,  δ_max::Float64, KepTransfers::Vector{KepTransfer}, p_ΔV::Float64,pump_down_only::Bool=true, flag_j2i::Bool=true)
    
    fb_to_node, curr_id = build_nodes_from_transfers(curr_id, δ_max, KepTransfers, KepTransfers, p_ΔV,pump_down_only, flag_j2i)
    return fb_to_node, curr_id
end

"""
    build_moon_change_edges(levels::Vector{VinfLevel}, curr_id::Int, flag_in2out::false)

Given a list of levels, it provides the nodes definitions linking transfers from levels belonging to two different moons with moon changes, if possible.
if flag_in2out is true, it includes also moon changes going from an inner to an outer moon
"""
function build_moon_change_nodes(levels::Vector{VinfLevel}, curr_id::Int, flag_in2out::Bool=false)

    fb_to_node = Dict() 

    for level_id in eachindex(levels)
        level = levels[level_id]
        # collect the maximum turning angle
        δmax_dep = level.δmax

        # collect the moon changes
        moon_change_set = level.moon_change_set
        
        for moon_change in moon_change_set
            α_dep = moon_change[1]
            
            level_arr = levels[Int(moon_change[2])]

            α_arr = moon_change[3]
            δmax_arr = level_arr.δmax

            # collect all transfers of level within δ_max from α_dep
            all_transfers_dep = vcat(level.res_transfers,
                                     level.nonres_transfers_I,
                                     level.nonres_transfers_O,
                                     level.vilts_incoming)

            dep_idxs = findall(tr -> abs(tr.arc2.α - α_dep) <= δmax_dep, all_transfers_dep)
            
            # collect all transfers of level_arr within δ_max from α_arr
            all_transfers_arr = vcat(level_arr.res_transfers,
                                     level_arr.nonres_transfers_I,
                                     level_arr.nonres_transfers_O,
                                     level_arr.vilts_outgoing)

            arr_idxs = findall(tr -> abs(tr.arc1.α - α_arr) <= δmax_arr, all_transfers_arr)

            # create nodes
            for idx_dep in dep_idxs
                tr_dep = all_transfers_dep[idx_dep]
                moon_dep_id = tr_dep.arc2.moon.id
                
                for idx_arr in arr_idxs
                    tr_arr = all_transfers_arr[idx_arr]
                    moon_arr_id = tr_arr.arc1.moon.id
                    if moon_arr_id > moon_dep_id 
                        fb_to_node[curr_id] = [tr_dep, tr_arr]
                        curr_id += 1
                    elseif (moon_arr_id < moon_dep_id) && flag_in2out
                        fb_to_node[curr_id] = [tr_dep, tr_arr]
                        curr_id += 1
                    end
               end
            end
        end
    end

    return fb_to_node, curr_id
end

"""
    update_wgts(p_ΔV::Float64, dvtof_list::Vector{Vector{Float64}}) -> Vector{Float64}

Recompute edge weights for a hopping graph using a new cost weight parameter.

This function allows re-weighting the graph edges without rebuilding the entire graph structure.
The weights are computed as a linear combination of ΔV (fuel) and time of flight (ToF).

# Arguments
- `p_ΔV::Float64`: Cost weight parameter between 0 and 1
  - `p_ΔV = 1.0`: Fuel-optimal (minimizes ΔV)
  - `p_ΔV = 0.0`: Time-optimal (minimizes ToF)
  - `0 < p_ΔV < 1`: Trade-off between fuel and time
- `dvtof_list::Vector{Vector{Float64}}`: Vector of [ΔV, ToF] pairs for each edge

# Returns
- `wgt_array::Vector{Float64}`: Vector of N edge weights computed as:
  ```
  wgt[i] = ΔV[i] * p_ΔV + ToF[i] * (1 - p_ΔV)
  ```

# Example
```julia
dvtof_list = [[2.5, 150.0], [1.8, 120.0], [3.2, 200.0]]
weights_fuel = update_wgts(0.95, dvtof_list)
weights_time = update_wgts(0.1, dvtof_list)
```
"""
function update_wgts(p_ΔV::Float64, dvtof_list::Vector{Vector{Float64}})
    wgt_array = [dv_tof[1] * p_ΔV + dv_tof[2] * (1 - p_ΔV) for dv_tof in dvtof_list]
    return wgt_array
end

# utility functions

function check_feas_fb_pumpdown(α₁, α₂, δ_max)
    return (α₂ != α₁) && (α₂ > α₁) && (α₂ - α₁ - δ_max <= 0)
end

function check_feas_fb(α₁, α₂, δ_max)
    return (α₂ != α₁) && (abs(α₂ - α₁) - δ_max <= 0)
end


