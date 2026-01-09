"""
Module containing the data structures and functions to discretize the design space
"""

mutable struct VinfLevel 
    id :: Integer
    moon :: Moon
    v∞ :: Float64
    v∞_dim :: Float64  # km/s
    δmax :: Float64  # rad

    res_list :: Vector{Any}
    α_list :: Vector{Any}

    res_transfers :: Vector{KepTransfer}
    nonres_transfers_I :: Vector{KepTransfer}
    nonres_transfers_O :: Vector{KepTransfer}
    vilts_outgoing :: Vector{KepTransfer}
    vilts_incoming :: Vector{KepTransfer}

    moon_change_set :: Vector{Any}
end

function VinfLevel_init(id::Integer, moon::Moon, vinf_dim::Float64)
    
    vinf = vinf_dim / moon.vga_dim
    δmax = 2*asin( moon.μ_dim / (moon.μ_dim + moon.min_rp_fb_dim*vinf_dim^2) )  # max turning angle
    res_list = Vector{Any}()
    α_list = Vector{Float64}()
    res_transfers = Vector{KepTransfer}()
    nonres_transfers_I = Vector{KepTransfer}()
    nonres_transfers_O = Vector{KepTransfer}()
    vilts_outgoing = Vector{KepTransfer}()
    vilts_incoming = Vector{KepTransfer}()
    moon_change_set = Vector{Any}()
    return VinfLevel(id, moon, vinf, vinf_dim, δmax,res_list, α_list, res_transfers, nonres_transfers_I, nonres_transfers_O, vilts_outgoing, vilts_incoming, moon_change_set)
end


struct DesignSpace
    id :: String
    moons :: Vector{Moon}             # list of the Moon objects in the system
    vinf_grid :: SMatrix{N,3,Float64} where N   # vinf grid, N moons × 3 columns (vinf1 [km/s], vinfamax [km/s], step [m/s])
    vinf_levels :: Vector{VinfLevel}  # list of the VinfLevel objects in the system

    Res_selection :: String  # "Auto" or "Manual"
    Res_list :: Vector{Vector{Vector{Int64}}}  # list of resonances for the manual selection [[N1,M1; N2,M2; ...], ... ], leave empty for automatic selection

    ΔV_max :: Float64  # maximum ΔV allowed for the vilts
    only_external_vilts :: Bool  # if true, only external vilts are computed
    choose_vilts :: Bool  # if true for any N-M (K) combination between two vinf levels only the lowest ΔV vilt is stored
    only_vinf_down :: Bool  # if true  only vilts going from a higher to a lower vinf are considered
    toll :: Float64  # tolerance for the non-resonant transfer and vilt computation

    total_bal :: Int
    total_vilts :: Int
    total_moonchanges :: Int

    id_to_transfer :: Dict{Int, KepTransfer}  # mapping from transfer id to KepTransfer object
end

function create_designspace(id, moons::Vector{Moon}, vinf_grid::SMatrix{N,3,Float64} where N, Res_selection::String, Res_list::Vector{Vector{Vector{Int64}}}, ΔV_max::Float64, only_external_vilts::Bool, choose_vilts::Bool, only_vinf_down::Bool, toll::Float64)

    # initialize id assingment
    level_id = 1
    transf_id = 1

    vinf_levels = Vector{VinfLevel}()

    total_bal = 0  # total number of ballistic transfers
    for ind in eachindex(moons)
        # create the vinf grid for the moon
        moon = moons[ind]
        vinf_span = vinf_grid[ind,:]
        vinf_min = vinf_span[1]
        vinf_max = vinf_span[2]
        vinf_step = vinf_span[3]
        vinf_values = vinf_min:vinf_step/1000:vinf_max
        
        # create the levels
        for v∞ in vinf_values
            level = VinfLevel_init(level_id, moon, v∞)

            if Res_selection == "Manual"
                cand_res_list =  Res_list[ind]  # add manually the resonances
                # compute the transfers and store them in the level
                level, transf_id = compute_resonant_transfers!(level, cand_res_list,transf_id)
                n_res = length(level.res_transfers)
            elseif Res_selection == "Auto"
                level = select_resonances!(level)  # automatic resonance selection and computation - to implemet
            end
                 
            level, transf_id = compute_nonresonant_transfers!(level, toll, transf_id )
            n_nonres = length(level.nonres_transfers_I) + length(level.nonres_transfers_O)
            # store the level and increase the id
            push!(vinf_levels, level)
            level_id += 1

            total_bal = total_bal + n_res + n_nonres 
        end

    end
    println("computed "*string(total_bal)*" ballistic transfers")
    # divide the levels per moon
    moon_to_levels = Dict{Moon, Vector{Int}}()
    for (i, level) in enumerate(vinf_levels)
        moon = level.moon
        if haskey(moon_to_levels, moon)
            push!(moon_to_levels[moon], i)
        else
            moon_to_levels[moon] = [i]
        end
    end

    moons = collect(keys(moon_to_levels))

    # compute vilts
    vinf_levels, total_vilts, transf_id = compute_vilts!(moons, moon_to_levels, vinf_levels, ΔV_max, only_external_vilts, choose_vilts, only_vinf_down, toll,transf_id + 1)
    # compute moon changes
    println("computed "*string(total_vilts)*" vilts")
    vinf_levels, total_moonchanges = compute_moon_changes!(moons, moon_to_levels, vinf_levels) 
    println("computed "*string(total_moonchanges)*" moon changes")
    
    # Build dictionary
    id_to_transfer = Dict{Int, KepTransfer}()
    for level in vinf_levels
        for tr in level.res_transfers
            id_to_transfer[tr.id] = tr
        end
        for tr in level.nonres_transfers_I
            id_to_transfer[tr.id] = tr
        end
        for tr in level.nonres_transfers_O
            id_to_transfer[tr.id] = tr
        end
        for tr in level.vilts_outgoing
            id_to_transfer[tr.id] = tr
        end
        for tr in level.vilts_incoming
            id_to_transfer[tr.id] = tr
        end
    end
    println("built id_to_transfer dictionary with "*string(length(id_to_transfer))*" transfers")
    
    
    
    return DesignSpace(id, moons, vinf_grid, vinf_levels, Res_selection, Res_list, ΔV_max, only_external_vilts, choose_vilts, only_vinf_down,toll, total_bal, total_vilts, total_moonchanges, id_to_transfer)
    # return DesignSpace(moons, vinf_grid, vinf_levels, Res_selection, Res_list, ΔV_max, only_external_vilts, choose_vilts, toll, total_bal,0,0)
end

function select_resonances!(level::VinfLevel)
    # automatic resonance selection based on the moon orbital period and the vinf
    # to implement
    return level
end

"""
    compute_resonant_transfers!(level::VinfLevel, cand_res_list::Vector{Vector{Int64}}, id::Integer = 0)
computes and stores the resonant transfers for the given level using its resonance list. Resonant transfers are stored as KepTransfer objects in the level.res_transfers field.
the corresponding pump angles are stored in the level.α_list field.

"""
function compute_resonant_transfers!(level::VinfLevel, cand_res_list::Vector{Vector{Int64}}, id::Integer = 0)
    # compute the resonant transfers for the level
    v∞ = level.v∞
    res_list = Vector{Any}()
    α_list = Vector{Float64}()
    res_transfers = Vector{Any}()

    for cand_res in cand_res_list
        N = cand_res[1]
        M = cand_res[2]
        # compute pump angle
        α = compute_pump_angle(N, M, v∞, level.moon.vga)
        
        if !isnan(α)
            # build the inbound arc
            arc_I = keparc_vinf(level.moon, [v∞, α, 0.0, 0.0])
            # build the inbound resonant transfer
            res_I = res(id, N, M, -1, arc_I)
            id = id + 1
            # build the outbound arc
            arc_O = keparc_vinf(level.moon, [v∞, α, π, 0.0])
            # build the outbound resonant transfer
            res_O = res(id, N, M, 1, arc_O)
            id = id + 1

            # store the transfers and the corresponding pump angle
            push!(res_transfers, res_I)
            push!(res_transfers, res_O)
            push!(α_list, α)
            push!(res_list, cand_res)
        end
    end
    # sort the pump angle list in ascending order and reorder res_list accordingly
    α_pairs = collect(zip(α_list, res_list))
    sort!(α_pairs, by = x -> x[1])                  
    α_list_sorted = [p[1] for p in α_pairs]
    res_list_sorted = [p[2] for p in α_pairs]

    sort!(res_transfers, by = x -> x.arc1.α)

    level.res_transfers = res_transfers   
    level.res_list = res_list_sorted
    level.α_list = α_list_sorted

    return level, id
end

"""
    compute_nonresonant_transfers!(level::VinfLevel, toll::Float64, id::Integer = 0)
computes and stores the non-resonant transfers for the given level using its resonance list. Non-resonant transfers are stored as KepTransfer objects in the level.nonres_transfers field.
"""
function compute_nonresonant_transfers!(level::VinfLevel, toll::Float64, id::Integer = 0)
    # compute the resonant transfers for the level
    v∞ = level.v∞
    for cand_res in level.res_list
        N = cand_res[1]
        M = cand_res[2]
        # compute IO non resonant transfer
        nonres_I = compute_nonres_transfer(level.moon, 0.0, v∞, M, N, -1, toll, id)
        if nonres_I.name != "Fail"
            # store the transfers and the corresponding pump angle
            push!(level.nonres_transfers_I, nonres_I)
            id = id + 1
        end

        # compute OI non resonant transfer
        nonres_O = compute_nonres_transfer(level.moon, 0.0, v∞, M, N, 1, toll, id)
        if nonres_O.name != "Fail"
            # store the transfers and the corresponding pump angle
            push!(level.nonres_transfers_O, nonres_O)
            id = id + 1
        end
    end
    return level, id
end

"""
    compute_vilts!(moons::Vector{Moon}, moon_to_levels::Vector{Any},levels::Vector{VinfLevel}, ΔV_max::Float64, only_external_vilts::Bool, choose_vilts::Bool, only_vinf_down::Bool,toll::Float64, id::Integer = 0)
computes and stores the VILT transfers for the levels. VILT transfers are stored as KepTransfer objects in the level.vilts_outgoing and level.vilts_incoming and level.vilts_outgoing fields.
For now, only N-M (K=0) vilts are computed.
"""

function compute_vilts!(moons::Vector{Moon}, moon_to_levels::Dict{Moon, Vector{Int}}, levels::Vector{VinfLevel}, ΔV_max::Float64, only_external_vilts::Bool, choose_vilts::Bool, only_vinf_down::Bool,toll::Float64, id::Integer = 0)
    K = 0  # only N-M (0) vilts for now

    total_vilts = 0
    curr_id = id
    # levels is your Vector{VinfLevel}
    # moon_to_levels = Dict{Moon, Vector{Int}}()
    # for (i, level) in enumerate(levels)
    #     moon = level.moon
    #     if haskey(moon_to_levels, moon)
    #         push!(moon_to_levels[moon], i)
    #     else
    #         moon_to_levels[moon] = [i]
    #     end
    # end
    # moons = collect(keys(moon_to_levels))
    for moon in moons
        println("looking for vilts around " * moon.name)
        level_indices = moon_to_levels[moon]
        for i in level_indices
            level_1 = levels[i]
            
            for j in i+1:maximum(level_indices)
                level_2 = levels[j]
                common_res = intersect(level_1.res_list, level_2.res_list)
                for res in common_res
                    N = res[1]
                    M = res[2]
                    # compute vilts between levels i and j
                    # external vilts
                    level_1, level_2, n_vilts, curr_id = combination_vilts!(level_1, level_2, N, M, K, 1, ΔV_max, choose_vilts, only_vinf_down,toll, curr_id + 1)
                    total_vilts += n_vilts
                    # internal vilts
                    if !only_external_vilts
                        level_1, level_2, n_vilts, curr_id = combination_vilts!(level_1, level_2, N, M, K, -1, ΔV_max, choose_vilts, only_vinf_down, toll, curr_id + 1)
                        total_vilts += n_vilts
                    end
                end

            end
        end
    end

    # sort the vilts_list in each level
    for level in levels
        sort!(level.vilts_outgoing, by = vilt -> vilt.arc1.α)
        sort!(level.vilts_incoming, by = vilt -> vilt.arc1.α)
    end
    return levels, total_vilts, curr_id
end

"""
    combination_vilts!(level1::VinfLevel, level2::VinfLevel,N::Integer, M::Integer, K::Integer, k_ei::Integer, ΔV_max::Float64, choose_vilts::Bool, only_vinf_down::Bool, toll::Float64, id::Integer=0)

Computes the 4 possible vilts between two vinf levels and then stores them. If choose_vilts is true, only the lowest ΔV vilt is stored and its direction is inverted to have both incoming and outgoing vilts.
"""

function combination_vilts!(level_1::VinfLevel, level_2::VinfLevel, N::Integer, M::Integer, K::Integer, k_ei::Integer, ΔV_max::Float64, choose_vilts::Bool, only_vinf_down::Bool,toll::Float64, id::Integer = 0)
    v∞_1 = level_1.v∞
    v∞_2 = level_2.v∞

    n_vilts = 0
    # external vilts
    # compute II vilt
    ΔVs = Vector{Float64}()
    combination_vilts_12 = Vector{KepTransfer}()

    vilt_II_12 = compute_vilt(level_1.moon, v∞_1, v∞_2, N, M, K, -1, -1, k_ei, 0.0, toll, id)
    if vilt_II_12.name != "Fail" && vilt_II_12.ΔV <= ΔV_max
        push!(ΔVs, vilt_II_12.ΔV)
        push!(combination_vilts_12, vilt_II_12)
        id = id + 1
    end

    # compute IO vilt
    vilt_IO_12 = compute_vilt(level_1.moon, v∞_1, v∞_2, N, M, K, -1, 1, k_ei, 0.0, toll, id)
    if vilt_IO_12.name != "Fail" && vilt_IO_12.ΔV <= ΔV_max
        push!(ΔVs, vilt_IO_12.ΔV)
        push!(combination_vilts_12, vilt_IO_12)
        id = id + 1
    end

    # compute OI vilt
    vilt_OI_12 = compute_vilt(level_1.moon, v∞_1, v∞_2, N, M, K, 1, -1, k_ei, 0.0, toll)
    if vilt_OI_12.name != "Fail" && vilt_OI_12.ΔV <= ΔV_max
        push!(ΔVs, vilt_OI_12.ΔV)
        push!(combination_vilts_12, vilt_OI_12)
        id = id + 1
    end

    # compute OO vilt
    vilt_OO_12 = compute_vilt(level_1.moon, v∞_1, v∞_2, N, M, K, 1, 1, k_ei, 0.0, toll, id)
    if vilt_OO_12.name != "Fail" && vilt_OO_12.ΔV <= ΔV_max
        push!(ΔVs, vilt_OO_12.ΔV)
        push!(combination_vilts_12, vilt_OO_12)
        id = id + 1
    end

    # if choose_vilts is true, only store the lowest ΔV vilt and invert it
    if choose_vilts && length(ΔVs) > 0
        min_index = argmin(ΔVs)
        min_vilt_12 = combination_vilts_12[min_index]
        if only_vinf_down
            if v∞_1 > v∞_2
                push!(level_1.vilts_outgoing, min_vilt_12)
                push!(level_2.vilts_incoming, min_vilt_12)
                n_vilts = 1
            end
        else
            push!(level_1.vilts_outgoing, min_vilt_12)
            push!(level_2.vilts_incoming, min_vilt_12)
            n_vilts = 1
        end
        

        min_vilt_21 = invert_vilt(min_vilt_12, toll, id)
        if min_vilt_21.name != "Fail"
            if only_vinf_down
                if v∞_1 < v∞_2
                    push!(level_1.vilts_incoming, min_vilt_21)
                    push!(level_2.vilts_outgoing, min_vilt_21)
                    n_vilts += 1
                end
            else
                push!(level_1.vilts_incoming, min_vilt_21)
                push!(level_2.vilts_outgoing, min_vilt_21)
                n_vilts += 1
            end
            
        end

    elseif length(ΔVs) > 0
        # store them all
        n_vilts = length(combination_vilts_12)
        for vilt in combination_vilts_12
            if only_vinf_down
                if v∞_1 > v∞_2
                    push!(level_1.vilts_outgoing, vilt)
                    push!(level_2.vilts_incoming, vilt)
                end
            else
                push!(level_1.vilts_outgoing, vilt)
                push!(level_2.vilts_incoming, vilt)
            end
            
            vilt_inv = invert_vilt(vilt, toll, id)
            id = id + 1
            if vilt_inv.name != "Fail"
                if only_vinf_down
                    if v∞_1 < v∞_2
                        push!(level_1.vilts_incoming, vilt_inv)
                        push!(level_2.vilts_outgoing, vilt_inv)
                        n_vilts += 1
                    end
                else
                    push!(level_1.vilts_incoming, vilt_inv)
                    push!(level_2.vilts_outgoing, vilt_inv)
                    n_vilts += 1
                end
                
            end
        end
    end

    return level_1, level_2, n_vilts, id
end

"""
    compute_moon_changes!(moons::Vector{Moon}, moon_to_levels::Vector{Any}, vinf_levels::Vector{VinfLevel})
Creates the ballistic moon change transfers between the levels.
"""
function compute_moon_changes!(moons::Vector{Moon}, moon_to_levels::Dict{Moon, Vector{Int}}, vinf_levels::Vector{VinfLevel})
    total_moonchanges = 0

    for t_moon in eachindex(moons)
        moon1 = moons[t_moon]
        levels1_idx = moon_to_levels[moon1]
        # cycle through the vinf levels of moon1
        for t_level1 in levels1_idx
            level1 = vinf_levels[t_level1]
            v∞1 = level1.v∞
            rga1_dim = moon1.rga_dim
            
            δmax₁  = level1.δmax
            α_list₁ = level1.α_list

            # cycle through the successive moons
            for t_Change in t_moon+1:length(moons)
                moon2 = moons[t_Change]
                rga2_dim = moon2.rga_dim
                levels2_idx = moon_to_levels[moon2]

                # cycle through the vinf levels of moon2
                for t_level2 in levels2_idx
                    level2 = vinf_levels[t_level2]
                    v∞2 = level2.v∞
                    δmax₂ = level2.δmax  # max turning angle
                    α_list₂ = level2.α_list  # resonances

                    # nonlinear system solution - non dimensional a1, a2 as variables
                    flag = false  # false if the transfer does not exist, true if it exists
                    # Tisserand constants
                    C1 = 3 - v∞1^2
                    C2 = 3 - v∞2^2
                    # 2nd order equation for a1
                    A =  rga1_dim*C1^2 - rga2_dim*C2^2 
                    B = 2 * ( rga2_dim^2*C2 - rga1_dim^2*C1 ) / rga1_dim
                    C = ( rga1_dim^3 - rga2_dim^3 ) / rga1_dim^2

                    a11 = ( -B  + sqrt(B^2  - 4*A*C) )/ (2*A)
                    a12 = ( -B  - sqrt(B^2  - 4*A*C) )/ (2*A)

                    # selection of the reasonable root and check if the transfer exists
                    if a11 > 0 && a12 <= 0
                        a1 = a11
                        a2 = rga1_dim*a1/rga2_dim
                        # pump angles computation
                        α₁ = compute_pump_angle(v∞1,a1)
                        α₂ = compute_pump_angle(v∞2,a2)

                        if !isnan(α₁) && !isnan(α₂)
                            flag = true
                        end

                    elseif a11 <= 0 && a12 > 0
                        a1 = a12
                        a2 = rga1_dim*a1/rga2_dim

                        α₁ = compute_pump_angle(v∞1,a1)
                        α₂ = compute_pump_angle(v∞2,a2)
                        if !isnan(α₁) && !isnan(α₂)
                            flag = true
                        end

                    else
                        a21 = rga1_dim*a11/rga2_dim
                        a22 = rga1_dim*a12/rga2_dim

                        α11 =compute_pump_angle(v∞1, a11) 
                        α12 =compute_pump_angle(v∞1, a12)
                        α21 =compute_pump_angle(v∞2, a21) 
                        α22 =compute_pump_angle(v∞2, a22) 

                        α11_check = isnan(α11)
                        α12_check = isnan(α12)
                        α21_check = isnan(α21)
                        α22_check = isnan(α22)

                        # check if the first solution makes sense
                        if !α11_check && !α21_check
                            flag = true
                            α₁ = α11
                            α₂ = α21
                            a1 = a11
                            a2 = a21

                        # check if the second solution makes sense
                        elseif !α12_check && !α22_check
                            flag = true
                            α₁ = α12
                            α₂ = α22
                            a1 = a12
                            a2 = a22
                        end
                    end
                    # store the transfer if it exists and if it is reachable by the resonances
                    if flag
                        δ_to_res₁ = [abs(α₁ - α_res) for α_res in α_list₁]
                        δ_to_res₂ = [abs(α₂ - α_res) for α_res in α_list₂]

                        if minimum(δ_to_res₁) <= δmax₁ && minimum(δ_to_res₂) <= δmax₂
                            total_moonchanges += 1
                            # moon changes are stored as [departure pump angle, arrival level id, arrival pump angle]
                            # store the transfer for the first moon
                            push!(level1.moon_change_set, [α₁,level2.id, α₂])

                            # store the transfer for the second moon
                            push!(level2.moon_change_set, [α₂,level1.id, α₁])
                            
                        end
                    end
                end
            end
        end
    end
    return vinf_levels, total_moonchanges
end