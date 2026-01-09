"""
Post-processing functions for shortest path solutions of the weigthed single objective graph search
"""



"""
    compute_fb_altitude(v∞_dim::Float64, α_in::Float64, α_out::Float64, moon::Moon)

Computes the fb altitude around the given moon for a planar fly by linking the input pump angles at the speficied v∞
"""
function compute_fb_altitude(v∞_dim::Float64, α_in::Float64, α_out::Float64, moon::Moon)
    
    # Turning angle
    δ = abs(α_out - α_in)
    

    # Pericenter radius and altitude
    if δ > 0
        rp_fb_dim = moon.μ_dim * (1/sin(δ/2) - 1) / v∞_dim^2
        FB_alt = rp_fb_dim - moon.radius_dim
        if FB_alt < moon.min_altfb_dim
            println("crash on $(moon.name), flyby at $(FB_alt) km of altitude")
        end
    else
        FB_alt = Inf
    end
    return FB_alt
end

"""
    compute_fb_altitude(arc_in::KepArc, arc_out::KepArc) -> Float64, Bool

Computes flyby altitude between two keplerian arcs.
"""
function compute_fb_altitude(arc_in::KepArc, arc_out::KepArc)
    
    # Use departure moon for moon changes
    moon = arc_in.moon
    moon_out = arc_out.moon
    
    if moon_out.id != moon.id
        mc_flag = true
        FB_alt = 0.0
    else
        mc_flag = false
        # Recover moon characteristics
        radius_dim = moon.radius_dim
        vga_dim = moon.vga_dim
        μ_dim = moon.μ_dim
        
        # Get dimensional v∞
        v∞_dim = arc_in.v∞ * vga_dim
        
        # Get pump angles
        α_in = arc_in.α
        α_out = arc_out.α
        
        FB_alt = compute_fb_altitude(v∞_dim, α_in, α_out, moon)

    end

    return FB_alt, mc_flag
end


"""
    compute_fb_altitude(transf_in::KepTransfer, transf_out::KepTransfer) -> Float64, Bool

Computes flyby altitude between two transfers.
"""
function compute_fb_altitude(transf_in::KepTransfer, transf_out::KepTransfer)
    # Get the arcs
    arc_in = transf_in.arc2
    arc_out = transf_out.arc1
   
    FB_alt, mc_flag = compute_fb_altitude(arc_in, arc_out)
    
    return FB_alt, mc_flag
end


"""
    compute_tour_fbs(dep_cond::SVector{2, Float64}, tar_cond::SVector{2, Float64}, path_seq::Vector{KepTransfer}, levels::Vector{VinfLevel}) -> Vector{KepTransfer}, Vector{Float64}

Computes flyby altitudes of a complete tour. Moon change transfers are also recovered and added to the tour
"""
function compute_tour_fbs(dep_cond::SVector{2, Float64}, tar_cond::SVector{2, Float64}, path_seq::Vector{KepTransfer}, levels::Vector{VinfLevel})
    
    # unpack departure and arrival conditions
    level_dep = levels[Int(dep_cond[1])]
    v∞_dim_dep = level_dep.v∞_dim
    α_dep = dep_cond[2] 
    moon_dep = level_dep.moon

    level_tar = levels[Int(tar_cond[1])]
    v∞_dim_tar = level_tar.v∞_dim
    α_tar = tar_cond[2] 
    moon_tar = level_tar.moon
    # initialize collectors
    tour_seq = Vector{KepTransfer}()
    FB_alt_array = Vector{Float64}()
    

    # compute first flyby
    
    transf_out = path_seq[1]
    moon_out = transf_out.arc1.moon

    if moon_dep.id == moon_out.id  # the tour starts with a same-moon transfer
        FB_alt = compute_fb_altitude(v∞_dim_dep, α_dep, transf_out.arc1.α, moon_dep)
        push!(tour_seq, transf_out)
        push!(FB_alt_array, FB_alt)

    else  # the tour starts with a moon change
        v∞_out_dim = transf_out.arc1.v∞*moon_out.vga_dim
        level_out_id = findfirst(l -> abs(l.v∞_dim - v∞_out_dim) < 1e-10 && l.moon.id == moon_out.id, levels)
        
        level_out = levels[level_out_id]
        k_io_dep = 1  # outbound from departure moon
        k_io_out = transf_out.k_io1  # same as first transfer

        mc_transf, FB_alt_in, FB_alt_out =  add_moon_change(level_dep, level_out, α_dep, k_io_dep, transf_out.arc1.α, k_io_out)
        
        # add fly bys and transfer to collectors
        push!(tour_seq, mc_transf)
        push!(FB_alt_array, FB_alt_in)

        push!(tour_seq, transf_out)
        push!(FB_alt_array, FB_alt_out)

    end

    # now iterate in the sequence
    transf_in = transf_out
    moon_in = transf_in.arc2.moon
    α_in = transf_in.arc2.α
    k_io_in = transf_in.k_io2
    v∞_dim_in = transf_in.arc2.v∞*moon_in.vga_dim

    for transf_out in path_seq[2:end]

        # recover arrival transfer characteristics
        moon_out = transf_out.arc1.moon
        α_out = transf_out.arc1.α
        k_io_out = transf_out.k_io1
        v∞_dim_out = transf_out.arc1.v∞*moon_out.vga_dim

        if moon_in.id == moon_out.id  # same-moon transfer
            FB_alt = compute_fb_altitude(v∞_dim_in, α_in, α_out, moon_out)
            push!(tour_seq, transf_out)
            push!(FB_alt_array, FB_alt)
        else
            # recover vinf levels
            level_in_id = findfirst(l -> abs(l.v∞_dim - v∞_dim_in) < 1e-10 && l.moon.id == moon_in.id, levels)
            level_in = levels[level_in_id]

            level_out_id = findfirst(l -> abs(l.v∞_dim - v∞_dim_out) <1e-10 && l.moon.id == moon_out.id, levels)
            level_out = levels[level_out_id]

            mc_transf, FB_alt_in, FB_alt_out =  add_moon_change(level_in, level_out, α_in, k_io_in, α_out, k_io_out)

            # add fly bys and transfer to collectors
            push!(tour_seq, mc_transf)
            push!(FB_alt_array, FB_alt_in)

            push!(tour_seq, transf_out)
            push!(FB_alt_array, FB_alt_out)
        end

        transf_in = transf_out
        # recover departure transfer characteristics
        moon_in = transf_in.arc2.moon
        α_in = transf_in.arc2.α
        k_io_in = transf_in.k_io2
        v∞_dim_in = transf_in.arc2.v∞*moon_in.vga_dim

    end

    # compute last flybys if needed
    if moon_in.id == moon_tar.id  # the tour ends with a same-moon transfer
        if α_tar != -1.0  # add last flyby if the tour ends with a specific pump angle
            FB_alt = compute_fb_altitude(v∞_dim_in, α_in, α_tar, moon_tar)
            push!(FB_alt_array, FB_alt)
        end

    else  # the tour ends with a moon change

        k_io_tar = k_io_in  # same as last transfer
        mc_transf, FB_alt_in, FB_alt_out =  add_moon_change(level_in, level_tar, α_in, k_io_in, α_tar, k_io_tar)
        
        # add fly bys and mc transfer to collectors
        push!(tour_seq, mc_transf)
        push!(FB_alt_array, FB_alt_in)
        
        # check if the tour ends with a specific pump angle or just at the level
        if isnan(FB_alt_out) == false
            push!(FB_alt_array, FB_alt_out)
        end

    end

    return tour_seq, FB_alt_array
end



"""
    compute_leg_data(seq::Vector{KepTransfer}, levels::Vector{VinfLevel}) -> Dict{String, Dict}

Computes leg data for each moon encountered in the transfer sequence.

# Arguments
- `seq::Vector{KepTransfer}`: Sequence of Keplerian transfers (including moon changes)
- `levels::Vector{VinfLevel}`: All v∞ levels in the design space

# Returns
- `Dict{String, Dict}`: Dictionary with moon names as keys, each containing:
  - `ΔV`: Total ΔV of transfers within that moon [m/s]
  - `ToF`: Total time of flight within that moon [days]
  - `vinf_values`: Array of v∞ values encountered [km/s]
  - `fb_altitudes`: Array of flyby altitudes [km]
  - `sequence`: Vector of KepTransfer objects for that moon

# Notes
- Moon change transfers appear in both legs (at the end of departure moon leg and at the start of arrival moon leg)
"""
function compute_leg_data(seq::Vector{KepTransfer}, levels::Vector{VinfLevel})
    
    # Initialize dictionary to store leg data per moon
    leg_data = Dict{String, Dict{Symbol, Any}}()
    
    if isempty(seq)
        return leg_data
    end
    
    # Helper function to initialize a new leg
    function init_leg()
        return Dict{Symbol, Any}(
            :ΔV => 0.0,
            :ToF => 0.0,
            :vinf_values => Float64[],
            :fb_altitudes => Float64[],
            :sequence => KepTransfer[]
        )
    end
    
    # Track current moon leg
    current_moon = seq[1].arc1.moon
    current_moon_name = current_moon.name
    
    # Initialize first leg
    if !haskey(leg_data, current_moon_name)
        leg_data[current_moon_name] = init_leg()
    end
    
    # Add first transfer's entry v∞
    push!(leg_data[current_moon_name][:vinf_values], seq[1].arc1.v∞ * current_moon.vga_dim)
    
    # Process the sequence
    for i in 1:length(seq)
        transf = seq[i]
        
        # Check if this is a moon change transfer
        is_moon_change = transf.arc1.moon.id != transf.arc2.moon.id
        
        if is_moon_change
            # Moon change - add to both departure and arrival moon legs
            moon_dep = transf.arc1.moon
            moon_arr = transf.arc2.moon
            moon_dep_name = moon_dep.name
            moon_arr_name = moon_arr.name
            
            # Add to departure moon leg (as last transfer)
            push!(leg_data[current_moon_name][:sequence], transf)
            leg_data[current_moon_name][:ΔV] += transf.ΔV
            leg_data[current_moon_name][:ToF] += transf.ToF
            push!(leg_data[current_moon_name][:vinf_values], transf.arc1.v∞ * moon_dep.vga_dim)
            
            # Initialize new leg for arrival moon if needed
            if !haskey(leg_data, moon_arr_name)
                leg_data[moon_arr_name] = init_leg()
            end
            
            # Add to arrival moon leg (as first transfer)
            push!(leg_data[moon_arr_name][:sequence], transf)
            leg_data[moon_arr_name][:ΔV] += transf.ΔV
            leg_data[moon_arr_name][:ToF] += transf.ToF
            push!(leg_data[moon_arr_name][:vinf_values], transf.arc2.v∞ * moon_arr.vga_dim)
            
            # Update current moon
            current_moon = moon_arr
            current_moon_name = moon_arr_name
            
        else
            # Regular transfer - belongs to a single moon
            moon = transf.arc1.moon
            moon_name = moon.name
            
            # Add transfer data
            push!(leg_data[current_moon_name][:sequence], transf)
            leg_data[current_moon_name][:ΔV] += transf.ΔV
            leg_data[current_moon_name][:ToF] += transf.ToF
            push!(leg_data[current_moon_name][:vinf_values], transf.arc2.v∞ * moon.vga_dim)
            
            # Compute flyby altitude if there's a next transfer on the same moon
            if i < length(seq)
                next_transf = seq[i+1]
                next_is_moon_change = next_transf.arc1.moon.id != next_transf.arc2.moon.id
                
                # Only compute flyby if next transfer is on the same moon (and not a moon change)
                # if !next_is_moon_change && next_transf.arc1.moon.id == moon.id
                FB_alt, _ = compute_fb_altitude(transf, next_transf)
                push!(leg_data[current_moon_name][:fb_altitudes], FB_alt)
                # end
            end
        end
    end
    
    return leg_data
end

function compute_orbit_insertion_ΔV()

end

function patch_tour()  # to make L continuous and select the moon changes

end



"""
    find_moon_change(level_dep::VinfLevel, level_arr_id::Int) -> Vector{Tuple{Float64, Float64}}

Returns moon change pump angle pair between two levels.

# Arguments
- `level_dep_id::Int`: Departure level ID
- `level_arr_id::Int`: Arrival level ID
- `levels::Vector{VinfLevel}`: All v∞ levels in the design space

# Returns
- Tuple{Float64, Float64}: (α_dep, α_arr) pair for the moon change between the levels
"""
function find_moon_change(level_dep::VinfLevel, level_arr_id::Int)
    
    # recover departure levels' moon changes
    moon_changes = level_dep.moon_change_set
    
    pairs = Tuple{Float64, Float64}[]
    npairs = 0
    # identify used moon change
    for moon_change in moon_changes

        α_dep = moon_change[1]
        mc_level_arr_id = Int(moon_change[2])
        α_arr = moon_change[3]
        
        if  mc_level_arr_id == level_arr_id
            push!(pairs, (α_dep, α_arr))
            npairs += 1
        end
    end

    # check if the definition was ambigous (analytically for any vinf, vinf2 only a moon change exists)
    if npairs  > 1
        error("moon changes incorrectly defined")
    elseif npairs == 0
        println(level_dep.id)
        println(level_arr_id)
        error("moon change does not exist")
    else
        return pairs[1]
    end
    
end


"""
    add_moon_change(level_1::VinfLevel, level_2::VinfLevel, α_1::Float64, k_io1::Int64, α_2::Float64, k_io2::Int64)

adds a moon change transfer and computes the needed flybys altitudes given the flybys before and after the moon change
"""
function add_moon_change(level_1::VinfLevel, level_2::VinfLevel, α_1::Float64, k_io1::Int64, α_2::Float64, k_io2::Int64)
    
    α_tup_mc = find_moon_change(level_1, level_2.id)

    # recover pump angles
    α_mc1 = α_tup_mc[1]
    α_mc2 = α_tup_mc[2]

    # build moon change model 
    # the tour starts with a moon change
    # Build moon change arcs

    # crank angles
    κ_1 = k_io1 == 1 ? π : 0.0
    κ_2 = k_io2 == 1 ? π : 0.0

    # Build first arc
    v∞_dim1 = level_1.v∞_dim
    moon1 = level_1.moon
    arc1 = keparc_vinf(moon1, [level_1.v∞, α_mc1, κ_1, 0.0])

    # Build second arc
    v∞_dim2 = level_2.v∞_dim
    moon2 = level_2.moon
    arc2 = keparc_vinf(moon2, [level_2.v∞, α_mc2, κ_2, 0.0])

    # Build moon change transfer
    mc_transf = moonchange(0, k_io1, k_io2, arc1, arc2, 0.0)

    # compute fb altitudes going in and out of the moon change
    FB_alt_in = compute_fb_altitude(v∞_dim1, α_1, α_mc1, moon1)
    if α_2 != -1.0 
        FB_alt_out = compute_fb_altitude(v∞_dim2, α_mc2, α_2, moon2)
    else
        FB_alt_out = NaN
    end

    return mc_transf, FB_alt_in, FB_alt_out
end


"""
    patch_tour(input_leg::Vector{KepTransfer}) -> Vector{KepTransfer}, Float64, Float64

Given a tour leg from the optimizer, checks dimensional ΔV and ToF and adjusts orbital parameters to make the leg continuous.

# Arguments
- `input_leg::Vector{KepTransfer}`: List of KepTransfer objects building the trajectory
- `levels::Vector{VinfLevel}`: All v∞ levels in the design space

# Returns
- `output_leg::Vector{KepTransfer}`: Continuous trajectory in Keplerian phase-free motion
- `ΔV_leg::Float64`: Total ΔV of the leg [m/s]
- `ToF_leg::Float64`: Total ToF of the leg [days]
"""
function patch_tour(input_leg::Vector{KepTransfer})
    
    if isempty(input_leg)
        return KepTransfer[], 0.0, 0.0
    end
    
    # Initialize output
    ΔV_leg = sum(transf.ΔV for transf in input_leg)
    ToF_leg = sum(transf.ToF for transf in input_leg)
    
    output_leg = Vector{KepTransfer}(undef, length(input_leg))
    
    # Process each transfer
    for t_model in 1:length(input_leg)
        
        curr_model_in = input_leg[t_model]
        curr_moon = curr_model_in.arc1.moon
        
        
        if t_model > 1
            L1 = output_leg[t_model-1].arc2.L
        else
            L1 = 0.0 
        end
        # Recover the arcs
        arc1_old = curr_model_in.arc1
        arc2_old = curr_model_in.arc2
        
        # Check for moon change
        if curr_moon.id == curr_model_in.arc2.moon.id

            # Same moon - rebuild the model for continuity
            α1_old = arc1_old.α
            κ1_old = arc1_old.κ
            v∞1 = arc1_old.v∞
            
            arc1_curr = keparc_vinf(curr_moon, [v∞1, α1_old, κ1_old, L1])
            
            α2_old = arc2_old.α
            κ2_old = arc2_old.κ
            v∞2 = arc2_old.v∞

            L2 = (L1 + curr_model_in.ToF*86400/(curr_moon.Tga_dim) * 2π) % (2π)
            
            arc2_curr = keparc_vinf(curr_moon, [v∞2, α2_old, κ2_old, L2])
            
            # Rebuild transfer with updated arcs
            if curr_model_in.N != 0 && curr_model_in.M != 0
                if curr_model_in.ΔV == 0
                    # Resonant or non-resonant transfer
                    if curr_model_in.k_io1 == curr_model_in.k_io2
                        curr_model = res(curr_model_in.id, curr_model_in.N, 
                                            curr_model_in.M, curr_model_in.k_io1,
                                            arc1_curr)
                    else
                        curr_model = nonres(curr_model_in.id, curr_model_in.N,
                                          curr_model_in.M, curr_model_in.k_io1,
                                          arc1_curr, arc2_curr, curr_model_in.ToF)
                    end
                else
                    # VILT transfer
                    curr_model = vilt(curr_model_in.id, curr_model_in.N, 
                                    curr_model_in.M, curr_model_in.K,
                                    curr_model_in.k_ei, curr_model_in.k_io1,
                                    curr_model_in.k_io2, arc1_curr, arc2_curr,
                                    curr_model_in.ΔV, curr_model_in.ToF)
                end
            else
                curr_model = curr_model_in  # Keep as is if special transfer
            end
            
            output_leg[t_model] = curr_model
            
        else
            # the transfer is a moon change - rebuild it for continuity and compute tof
            prev_moon = curr_moon
            curr_moon = curr_model_in.arc2.moon

            # recover i/o info
            k_io_dep = curr_model_in.k_io1
            k_io_arr = curr_model_in.k_io2

            # Build departure arc
            α1 = arc1_old.α
            κ1 = arc1_old.κ
            v∞1 = arc1_old.v∞
            
            arc1 = keparc_vinf(prev_moon, [v∞1, α1, κ1, L1])
                       
            # recover mockup arrival arc to get sma and e
            arc2_mockup = arc2_old
            α2 = arc2_mockup.α
            κ2 = arc2_mockup.κ
            v∞2 = arc2_mockup.v∞            
            # Moon change ToF computation
            a1 = arc1.sma
            e1 = arc1.ecc
            E1 = -cos(κ1) * acos((a1 - 1) / (e1 * a1))
            τ1_dim = a1^(3/2) / (2π) * (E1 - e1 * sin(E1)) * prev_moon.Tga_dim
            
            a2 = arc2_mockup.sma
            e2 = arc2_mockup.ecc
            E2 = -cos(κ2) * acos((a2 - 1) / (e2 * a2))
            τ2_dim = a2^(3/2) / (2π) * (E2 - e2 * sin(E2)) * curr_moon.Tga_dim
            
            # True anomaly of 2nd encounter

            p = a1 * (1 - e1^2) * prev_moon.rga_dim
            cos_f = 1/e1 * (p / curr_moon.rga_dim - 1)
            f2 = k_io_arr * acos(cos_f)
            
            # Time of flight
            ToF_change_dim = τ2_dim - τ1_dim
            if ToF_change_dim <= 0
                ToF_change_dim = ToF_change_dim + a2^(3/2) * curr_moon.Tga_dim
            end
            # Build arc2 with actual f
            Kep_arr = [a1 * prev_moon.rga_dim / curr_moon.rga_dim, e1, 
                      arc1.inc, arc1.Ω, arc1.ω, f2]
            
            # Position recovery
            r2, v2 = Kep_To_Car(Kep_arr...)
            
            # Longitude of encounter
            L2_arr = atan(r2[2], r2[1])
            
            # Build arrival arc with actual longitude of encounter
            arc2 = keparc_vinf(curr_moon, [v∞2, α2, κ2, L2_arr])
            
            # Build moon change transfer
            moon_change_model = moonchange(0, k_io_dep, k_io_arr, arc1, arc2, ToF_change_dim)
            
            output_leg[t_model] =  moon_change_model
            
            
        end
    end
    
    return output_leg, ΔV_leg, ToF_leg
end

"""
    patch_solution(sol::WeightedHGSol)

Returns the continous tour sequences of a solution to a graph search problem
"""
function patch_solution(sol::WeightedHGSol)

    dep_cond = sol.prb.dep_cond
    tar_cond = sol.prb.tar_cond
    vinf_levels = sol.prb.des_space.vinf_levels
    seqs_array = [path.seq for path in sol.paths]

    n_seqs = length(seqs_array)
    cont_tours = Vector{Vector{KepTransfer}}(undef, n_seqs) 
    for t = 1:n_seqs
        tour_seq, FB_alt_array = compute_tour_fbs(dep_cond, tar_cond, seqs_array[t], vinf_levels)
        cont_tour, dv_check, tof_check = patch_tour(tour_seq)
        cont_tours[t] = cont_tour
    end
    return cont_tours
end


"""
    patch_solution(sol::WeightedHGSol)

Returns the continous tour sequences of a solution to a graph search problem
"""
function patch_solution(sol::MOHGSol)

    dep_cond = sol.prb.dep_cond
    tar_cond = sol.prb.tar_cond
    vinf_levels = sol.prb.des_space.vinf_levels
    seqs_array = [path.seq for path in sol.paths]

    n_seqs = length(seqs_array)
    cont_tours = Vector{Vector{KepTransfer}}(undef, n_seqs) 
    for t = 1:n_seqs
        tour_seq, FB_alt_array = compute_tour_fbs(dep_cond, tar_cond, seqs_array[t], vinf_levels)
        cont_tour, dv_check, tof_check = patch_tour(tour_seq)
        cont_tours[t] = cont_tour
    end
    return cont_tours
end