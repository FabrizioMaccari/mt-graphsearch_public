"""
Plotting functions to visualize tours on the v-inf circles of the encountered moons
"""

"""
    plot_transfer_circle(moon::Moon, transfers::Vector{KepTransfer}, resonances::Vector{Vector{Int}}, filepath::String; flag::Int=0)

Plots the sequence of same body transfers on a planar section of the v∞ sphere.

# Arguments
- `transfers::Vector{KepTransfer}`: List containing the transfer objects to represent
- `resonances::Vector{Vector{Int}}`: List of resonances locii to plot [[N1,M1], ..., [Nn,Mn]]
- `filepath::String`: Path where to save the figure
- `flag::Int=0`: Moon change handling flag
  - 0: plot the moon change transfers on the last vinf level
  - 1: plot the moon change transfers on the first vinf level

# Returns
- `fig`: Figure object representing the evolution of the v∞ vector

# Notes
- Represents the v∞ vector evolution with respect to the selected Moon
- Shows locii of resonances encountered in the transfers
- Links to departure/target angles if provided
- Legend is placed outside the plot area
- v∞ circles are color-coded according to their dimensional value
"""
function plot_transfer_circle(transfers::Vector{KepTransfer}, resonances::Vector{Vector{Int}}, filepath::String, α_dep::Float64=NaN, α_tar::Float64=NaN)
    
    # Create figure
    star_marker_size = 12
    cmap = :managua
    cmap2 = :glasbey_bw_minc_20_maxl_70_n256 
    fig = Figure(size=(1000, 1200))
    ax = Axis(fig[1, 1],
        xlabel="p₁ [-], rga",
        ylabel="p₂ [-], vga",
        aspect=DataAspect(),
        xlabelsize=25,
        ylabelsize=25,
        xticklabelsize=25,
        yticklabelsize=25,
        
    )
    
    first_model = transfers[1]
    last_model = transfers[end]

    # Moon velocity recovery
    moon = first_model.arc1.moon
    vga = moon.vga
    vga_dim = moon.vga_dim
    
    # Initialize variables
    t = 0
    V∞_old = 0.0
    flag_arr = 0
    flag_dep = 0
    
    # Moon change handling    
    # Determine transfer sequence based on moon changes
    if first_model.arc1.moon.id == first_model.arc2.moon.id && 
       last_model.arc1.moon.id == last_model.arc2.moon.id
        # No moon change
        transfers_hop = transfers
        
    elseif first_model.arc1.moon.id != first_model.arc2.moon.id && 
           last_model.arc1.moon.id == last_model.arc2.moon.id
        # Moon change only at beginning
        transfers_hop = transfers[2:end]
       
        V∞_vec_in = v∞_vec_from_α_κ(first_model.arc2)

        V∞_in = first_model.arc2.v∞
        
        V∞_vec_1 = v∞_vec_from_α_κ(transfers[2].arc1)

        V∞_1 = transfers[2].arc1.v∞
        
        scatter!(ax, [V∞_vec_in[1]], [V∞_vec_in[2]], 
                marker=:rtriangle, markersize=12, color=:green,
                label="Arrival at Moon")
        
        flag_arr = 1
        
    elseif first_model.arc1.moon.id == first_model.arc2.moon.id && 
           last_model.arc1.moon.id != last_model.arc2.moon.id
        # Moon change only at end
        transfers_hop = transfers[1:end-1]
        
        V∞_vec_final = v∞_vec_from_α_κ(last_model.arc1)
        
        scatter!(ax, [V∞_vec_final[1]], [V∞_vec_final[2]],
                marker=:ltriangle, markersize=12, color=:red,
                label="Departure from Moon")
        
        flag_dep = 1
        
    else
        # Moon change at both ends
        transfers_hop = transfers[2:end-1]
        
        V∞_vec_in = v∞_vec_from_α_κ(first_model.arc2)
        V∞_in = first_model.arc2.v∞
        
        V∞_vec_1 = v∞_vec_from_α_κ(transfers[2].arc1)
        V∞_1 = transfers[2].arc1.v∞
        
        scatter!(ax, [V∞_vec_in[1]], [V∞_vec_in[2]],
                marker=:rtriangle, markersize=12, color=:green,
                label="Arrival at Moon")
        
        # Plot flyby arc using Rodrigues formula
        h_hyp_dir = cross([V∞_vec_in..., 0], [V∞_vec_1..., 0])
        h_hyp_dir = h_hyp_dir ./ norm(h_hyp_dir)
        
        δ_a = acos(dot(V∞_vec_in, V∞_vec_1) / (V∞_in * V∞_1))
        deltas = 0:0.5π/180:δ_a
        
        V∞_plot_arr = zeros(length(deltas), 2)
        for (i, δ) in enumerate(deltas)
            v_rot = [V∞_vec_in..., 0] .* cos(δ) .+ 
                    cross(h_hyp_dir, [V∞_vec_in..., 0]) .* sin(δ) .+
                    h_hyp_dir .* dot(h_hyp_dir, [V∞_vec_in..., 0]) .* (1 - cos(δ))
            V∞_plot_arr[i, :] = v_rot[1:2]
        end
        
        lines!(ax, V∞_plot_arr[:, 1], V∞_plot_arr[:, 2], color=:black, linewidth=1)
        
        flag_arr = 1
        
        V∞_vec_final = v∞_vec_from_α_κ(last_model.arc1)
        scatter!(ax, [V∞_vec_final[1]], [V∞_vec_final[2]],
                marker=:ltriangle, markersize=12, color=:red,
                label="Departure from Moon")
        
        flag_dep = 1
    end
    
    # Process same-moon sequence
    V∞_plot = Vector{Vector{Float64}}()
    V∞_levels = Float64[]
    V∞_vec_old = zeros(2)
    V∞_vec_2 = zeros(2)  # Initialize
    V∞_2 = 0.0  # Initialize
    
    for (idx, transfer) in enumerate(transfers_hop)
        # Check same moon
        if transfer.arc1.moon.id != transfer.arc2.moon.id
            error("Trajectories referred to different moons")
        end
        
        # v∞ vectors
        V∞_vec_1 = v∞_vec_from_α_κ(transfer.arc1)
        V∞_vec_2 = v∞_vec_from_α_κ(transfer.arc2)
        
        V∞_1 = transfer.arc1.v∞
        V∞_2 = transfer.arc2.v∞
        
        L1 = transfer.arc1.L
        L2 = transfer.arc2.L
        
        # Link to previous arc
        if idx > 1 && V∞_old ≈ V∞_1
            h_hyp_dir = cross([V∞_vec_old..., 0], [V∞_vec_1..., 0])
            norm_h = norm(h_hyp_dir)
            
            if norm_h > 1e-10
                h_hyp_dir = h_hyp_dir ./ norm_h
                
                δ_a = acos(clamp(dot(V∞_vec_old, V∞_vec_1) / (V∞_old * V∞_1), -1.0, 1.0))
                deltas = 0:0.5π/180:δ_a
                
                for δ in deltas
                    v_rot = [V∞_vec_old..., 0] .* cos(δ) .+
                            cross(h_hyp_dir, [V∞_vec_old..., 0]) .* sin(δ) .+
                            h_hyp_dir .* dot(h_hyp_dir, [V∞_vec_old..., 0]) .* (1 - cos(δ))
                    push!(V∞_plot, v_rot[1:2])
                end
            end
        end
        
        # Plot transfer
        if V∞_1 ≈ V∞_2 && L1 ≈ L2  # Ballistic resonant flyby
            if idx == 1
                isempty(V∞_levels) && push!(V∞_levels, V∞_1)
                
                if flag_arr == 0 && isnan(α_dep)
                    scatter!(ax, [V∞_vec_1[1]], [V∞_vec_1[2]],
                            marker=:x, markersize=10, color=:black,
                            label="First v∞")
                else
                    scatter!(ax, [V∞_vec_1[1]], [V∞_vec_1[2]],
                            marker=:star5, markersize=star_marker_size, color=:black)
                end
                
                scatter!(ax, [V∞_vec_2[1]], [V∞_vec_2[2]],
                        marker=:star5, markersize=star_marker_size, color=:black,
                        label="v∞ sequence")
                push!(V∞_plot, V∞_vec_2)
                
            else
                scatter!(ax, [V∞_vec_1[1]], [V∞_vec_1[2]],
                        marker=:star5, markersize=star_marker_size, color=:black)
                scatter!(ax, [V∞_vec_2[1]], [V∞_vec_2[2]],
                        marker=:star5, markersize=star_marker_size, color=:black)
            end
            
        else  # Leveraging or non-resonant transfer
            if idx == 1
                isempty(V∞_levels) && push!(V∞_levels, V∞_1)
                
                if flag_arr == 0 && isnan(α_dep)
                    scatter!(ax, [V∞_vec_1[1]], [V∞_vec_1[2]],
                            marker=:x, markersize=10, color=:black,
                            label="First v∞")
                else
                    scatter!(ax, [V∞_vec_1[1]], [V∞_vec_1[2]],
                            marker=:star5, markersize=star_marker_size, color=:black)
                end
                
                scatter!(ax, [V∞_vec_2[1]], [V∞_vec_2[2]],
                        marker=:star5, markersize=star_marker_size, color=:black,
                        label="v∞ sequence")
                push!(V∞_plot, V∞_vec_1)
                push!(V∞_plot, V∞_vec_2)
                
            else
                scatter!(ax, [V∞_vec_1[1]], [V∞_vec_1[2]],
                        marker=:star5, markersize=star_marker_size, color=:black)
                scatter!(ax, [V∞_vec_2[1]], [V∞_vec_2[2]],
                        marker=:star5, markersize=star_marker_size, color=:black)
                push!(V∞_plot, V∞_vec_1)
                push!(V∞_plot, V∞_vec_2)
                
                if !isempty(V∞_levels) && !(V∞_1 ≈ V∞_levels[end])
                    push!(V∞_levels, V∞_1)
                end
            end
            
            if !(V∞_1 ≈ V∞_2)
                push!(V∞_levels, V∞_2)
            end
        end
        
        V∞_vec_old = V∞_vec_2
        V∞_old = V∞_2
    end
    
    # Final point (only if no target condition and transfers_hop not empty)
    if !isempty(transfers_hop)
        if flag_dep == 0 && isnan(α_tar)
            scatter!(ax, [V∞_vec_2[1]], [V∞_vec_2[2]],
                    marker=:circle, markersize=10, color=:black,
                    label="Last v∞")
        elseif flag_dep == 1 && isnan(α_tar)
            h_hyp_dir = cross([V∞_vec_2..., 0], [V∞_vec_final..., 0])
            norm_h = norm(h_hyp_dir)
            
            if norm_h > 1e-10
                h_hyp_dir = h_hyp_dir ./ norm_h
                
                δ_a = acos(clamp(dot(V∞_vec_2, V∞_vec_final) / (V∞_2 * V∞_2), -1.0, 1.0))
                deltas = 0:0.5π/180:δ_a
                
                for δ in deltas
                    v_rot = [V∞_vec_2..., 0] .* cos(δ) .+
                            cross(h_hyp_dir, [V∞_vec_2..., 0]) .* sin(δ) .+
                            h_hyp_dir .* dot(h_hyp_dir, [V∞_vec_2..., 0]) .* (1 - cos(δ))
                    push!(V∞_plot, v_rot[1:2])
                end
            end
        end
    end
    
    # v∞ circle preparation
    α_vec = range(-π, π, length=360)
    X_vec = cos.(α_vec)
    Y_vec = sin.(α_vec)
    
    # κ angles
    κ_vec = [0.0, π]
    
    # Determine if this is a single-level sequence
    single_level = (length(V∞_levels) == 1)
    
    # Determine colormap range based on dimensional v∞ values
    if !isempty(V∞_levels)
        V∞_dim_min = minimum(V∞_levels) * vga_dim
        V∞_dim_max = maximum(V∞_levels) * vga_dim
        
        # Plot v∞ circles with color coding or label
        for V∞ in V∞_levels
            V∞_dim = V∞ * vga_dim
            
            if single_level
                color_val = 0.01
                # Single level: use label with v∞ value
                lines!(ax, V∞ .* X_vec, V∞ .* Y_vec,
                      color=get(colorschemes[cmap], color_val), linewidth=2,
                      label="v∞ = $(round(V∞_dim, digits=3)) km/s")
            else
                # Multiple levels: use color coding
                # Normalize to [0, 1] for colormap
                if V∞_dim_max - V∞_dim_min != 0
                    color_val = (V∞_dim - V∞_dim_min) / (V∞_dim_max - V∞_dim_min)
                else
                    color_val = 0.5
                end
                circle_color = get(colorschemes[cmap], color_val)
                
                lines!(ax, V∞ .* X_vec, V∞ .* Y_vec,
                      color=circle_color, linewidth=2)
            end
        end
        
        # Add colorbar only for multi-level sequences
        if !single_level
            Colorbar(fig[1, 2],
                    colormap=cmap,
                    limits=(V∞_dim_min, V∞_dim_max),
                    label="v∞ [km/s]",
                    labelsize=25,
                    ticklabelsize=25,
                    height=Relative(0.7))
        end
    end

    # Link to departure condition if provided
    if !isnan(α_dep) && !isempty(transfers_hop)
        # Get first transfer's v∞
        first_transfer = transfers_hop[1]
        V∞_first = first_transfer.arc1.v∞
        V∞_vec_first = v∞_vec_from_α_κ(first_transfer.arc1)
        # Departure v∞ vector
        V∞_vec_dep = v∞_vec_from_α_κ(V∞_first,α_dep, first_transfer.arc1.κ)
        
        # Plot departure point
        scatter!(ax, [V∞_vec_dep[1]], [V∞_vec_dep[2]],
                marker=:diamond, markersize=12, color=:blue,
                label="Departure condition")
        
        # Plot arc from departure to first transfer using Rodrigues formula
        h_hyp_dir = cross([V∞_vec_dep..., 0], [V∞_vec_first..., 0])
        norm_h = norm(h_hyp_dir)
        
        if norm_h > 1e-10  # Only plot if there's a rotation
            h_hyp_dir = h_hyp_dir ./ norm_h
            
            δ_a = acos(clamp(dot(V∞_vec_dep, V∞_vec_first) / (V∞_first * V∞_first), -1.0, 1.0))
            deltas = range(0, δ_a, length=100)
            
            dep_arc = zeros(length(deltas), 2)
            for (i, δ) in enumerate(deltas)
                v_rot = [V∞_vec_dep..., 0] .* cos(δ) .+
                        cross(h_hyp_dir, [V∞_vec_dep..., 0]) .* sin(δ) .+
                        h_hyp_dir .* dot(h_hyp_dir, [V∞_vec_dep..., 0]) .* (1 - cos(δ))
                dep_arc[i, :] = v_rot[1:2]
            end
            
            lines!(ax, dep_arc[:, 1], dep_arc[:, 2], 
                  color=:black, linewidth=2, linestyle=:dash,
                  label="Departure flyby")
        end
    end
    
    # Link to target condition if provided
    if !isnan(α_tar) && !isempty(transfers_hop) && α_tar != -1.0
        # Get last transfer's v∞
        V∞_last = V∞_2
        
        # Target v∞ vector
        if V∞_vec_2[2] > 0
            V∞_vec_tar = v∞_vec_from_α_κ(V∞_last, α_tar, π)
        else
            V∞_vec_tar = v∞_vec_from_α_κ(V∞_last, α_tar, 0.0)
        end
        # Plot target point
        scatter!(ax, [V∞_vec_tar[1]], [V∞_vec_tar[2]],
                marker=:diamond, markersize=12, color=:orange,
                label="Target condition")
        
        # Plot arc from last transfer to target using Rodrigues formula
        h_hyp_dir = cross([V∞_vec_2..., 0], [V∞_vec_tar..., 0])
        norm_h = norm(h_hyp_dir)
        
        if norm_h > 1e-10  # Only plot if there's a rotation
            h_hyp_dir = h_hyp_dir ./ norm_h
            
            δ_a = acos(clamp(dot(V∞_vec_2, V∞_vec_tar) / (V∞_last * V∞_last), -1.0, 1.0))
            deltas = range(0, δ_a, length=100)
            
            tar_arc = zeros(length(deltas), 2)
            for (i, δ) in enumerate(deltas)
                v_rot = [V∞_vec_2..., 0] .* cos(δ) .+
                        cross(h_hyp_dir, [V∞_vec_2..., 0]) .* sin(δ) .+
                        h_hyp_dir .* dot(h_hyp_dir, [V∞_vec_2..., 0]) .* (1 - cos(δ))
                tar_arc[i, :] = v_rot[1:2]
            end
            
            lines!(ax, tar_arc[:, 1], tar_arc[:, 2],
                  color=:black, linewidth=2, linestyle=:dash,
                  label="Arrival flyby")
        end
    end
    

    # Plot resonances from actual transfers using N:M values
    # Use hawaiiS colormap for resonances
    n_resonances = length(resonances)
    if n_resonances > 0
        # Get colors from the hawaiiS colormap
        hawaii_colors = [get(colorschemes[cmap2], i / max(n_resonances, 2)) for i in 1:n_resonances]
    end
    
    for (res_idx, res) in enumerate(resonances)
        N = res[1]
        M = res[2]
        
        V∞_vec_res = Vector{Vector{Float64}}()
        
        for V∞ in V∞_levels
            # Compute pump angle for this N:M resonance at this v∞ level
            α_res = compute_pump_angle(N, M, V∞, vga)
            
            if !isnan(α_res)
                # Compute v∞ vectors for both κ values
                for κ in κ_vec
                    V∞_vec = v∞_vec_from_α_κ(V∞, α_res, κ)
                    push!(V∞_vec_res, V∞_vec)
                end
            end
        end
        
        if !isempty(V∞_vec_res)
            res_x = [v[1] for v in V∞_vec_res]
            res_y = [v[2] for v in V∞_vec_res]
            scatter!(ax, res_x, res_y, 
                    marker=:circle, 
                    markersize=8,
                    color=hawaii_colors[res_idx],
                    label="$N:$M Res.")
        end
    end
    
    # Plot v∞ history
    if !isempty(V∞_plot)
        plot_x = [v[1] for v in V∞_plot]
        plot_y = [v[2] for v in V∞_plot]
        lines!(ax, plot_x, plot_y, color=:black, linewidth=1)
    end
    
    # Axis setup
    if !isempty(V∞_levels)
        radius = maximum(V∞_levels) * 1.05
    else
        radius = 1.0
    end
    xlims!(ax, -radius, radius)
    ylims!(ax, -radius, radius)
    
    # Legend outside the plot area
    if single_level
        # Single level: legend below plot
        Legend(fig[2, 1], ax, 
               orientation=:horizontal,
               tellwidth=false,
               tellheight=true,
               labelsize=25,
               nbanks=3)
    else
        # Multiple levels: legend below plot and colorbar
        Legend(fig[2, 1:2], ax, 
               orientation=:horizontal,
               tellwidth=false,
               tellheight=true,
               labelsize=25,
               nbanks=3)
    end
    
    # Save figure
    save(filepath, fig, px_per_unit=2)
    
    return fig
end


"""
    plot_tour_circles(dep_cond::SVector{2,Float64}, tar_cond::SVector{2,Float64}, tour::Vector{KepTransfer}, levels::Vector{VinfLevel}, folder::String)

Plots the Tour sequences of same body transfers on planar section of the v∞ sphere.
Works only with tours connected with moon change transfers.
Saves one figure per moon leg in the specified folder.

# Arguments
- `dep_cond::SVector{2,Float64}`: Departure condition [level_id, α]
- `tar_cond::SVector{2,Float64}`: Target condition [level_id, α]
- `tour::Vector{KepTransfer}`: List containing the transfer objects to represent
- `levels::Vector{VinfLevel}`: All v∞ levels in the design space
- `folder::String`: Folder path where to save the figures

# Notes
- Creates separate sequences for the legs of the tour at different moons
- Saves figures with moon names as filenames
- Each figure shows the evolution of the v∞ vector and locii of resonances
- Links first leg to departure condition and last leg to target condition
"""
function plot_tour_circles(dep_cond::SVector{2,Float64}, tar_cond::SVector{2,Float64}, tour::Vector{KepTransfer}, levels::Vector{VinfLevel}, folder::String)
    
    # Create folder if it doesn't exist
    if !isdir(folder)
        mkpath(folder)
    end
    
    # Get departure and target information
    level_dep = levels[Int(dep_cond[1])]
    α_dep = dep_cond[2]
    moon_dep = level_dep.moon
    
    level_tar = levels[Int(tar_cond[1])]
    α_tar = tar_cond[2]
    moon_tar = level_tar.moon
    
    # Create separate sequences for the legs of the tour at different moons
    tour_split_legs = Vector{Vector{KepTransfer}}()
    curr_leg = KepTransfer[]
    
    for transfer in tour
        curr_moon_id1 = transfer.arc1.moon.id
        curr_moon_id2 = transfer.arc2.moon.id
        
        if curr_moon_id1 == curr_moon_id2  # no moon change, continue current leg
            push!(curr_leg, transfer)
        else  # moon change detected
            push!(curr_leg, transfer)  # add moon change to current leg to end it
            push!(tour_split_legs, curr_leg)  # save leg
            curr_leg = KepTransfer[]  # empty curr_leg collector to save next moon leg
        end
    end
    
    push!(tour_split_legs, curr_leg)  # append last leg
    
    # Create v∞ circle plot for each leg
    for (leg_idx, leg) in enumerate(tour_split_legs)
        if isempty(leg)
            continue
        end
        
        curr_moon = leg[1].arc2.moon
        
        # Extract resonances from the leg (excluding last transfer which might be moon change)
        resonances = Vector{Vector{Int}}()
        for model in leg[1:end-1]
            if model.N != 0 && model.M != 0  # Only add if valid resonance
                push!(resonances, [model.N, model.M])
            end
        end
        
        # Remove duplicates
        unique!(resonances)
        
        # Create filepath with moon name
        filepath = joinpath(folder, "$(curr_moon.name)_vinf_circle.png")
        
        # Determine if this is first or last leg and pass boundary conditions
        is_first_leg = (leg_idx == 1) && (curr_moon.id == moon_dep.id)
        is_last_leg = (leg_idx == length(tour_split_legs)) && (curr_moon.id == moon_tar.id)
        
        dep_α = is_first_leg ? α_dep : NaN
        tar_α = is_last_leg ? α_tar : NaN
        
        # Plot and save
        fig = plot_transfer_circle(leg, resonances, filepath, dep_α, tar_α)
        
        println("Saved v∞ circle plot for $(curr_moon.name) to $filepath")
    end
    
    return nothing
end