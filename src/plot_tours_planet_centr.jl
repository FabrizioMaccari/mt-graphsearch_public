"""
Plotting functions for tour sequences in the planet centred inertial ref. frame
"""

"""
    plot_moon_tour_2d(moons::Vector{Moon}, transfers::Vector{KepTransfer}, filepath::String; points::Int=300)

Given a list of transfers and a list of moons, it assembles and saves the plot of the Tour in the planet centred inertial frame.

# Arguments
- `moons::Vector{Moon}`: List of Moons to insert in the plot
- `transfers::Vector{KepTransfer}`: List of transfer objects to insert in the plot
- `filepath::String`: Path where to save the figure
- `points::Int=300`: Number of points to plot per transfer

# Notes
- Plots the full trajectory in gray/blue color
- Moon orbits are plotted with their respective colors
- Flyby locations are marked with dots
- Figure is saved to the specified filepath
"""
function plot_moon_tour_2d(moons::Vector{Moon}, transfers::Vector{KepTransfer}, filepath::String; points::Int=300)
    
    # Create figure
    fig = Figure(size = (700, 700))
    ax = Axis(fig[1, 1],
        xlabel="x [km]",
        ylabel="y [km]",
        aspect=DataAspect(),
        xlabelsize=16,
        ylabelsize=16,
        xticklabelsize=16,
        yticklabelsize=16
    )
    
    r_a_max = 0.0  # max distance from primary for scaling
    
    # Arrays to store trajectory points for continuous plotting
    all_x = Float64[]
    all_y = Float64[]
    fb_x = Float64[]  # flyby locations
    fb_y = Float64[]
    
    # Process each transfer
    for (t, transfer) in enumerate(transfers)
        
        # Check if moon change
        is_moon_change = transfer.arc1.moon.id != transfer.arc2.moon.id
        
        # Recover Keplerian parameters
        Kep1 = copy(transfer.arc1.coe)
        Kep2 = copy(transfer.arc2.coe)
        f1 = Kep1[6]
        f2 = Kep2[6]
        
        # Max apocenter for scaling
        r_a_model = max(Kep1[1] * (1 + Kep1[2]), Kep2[1] * (1 + Kep2[2]))
        
        # Arc limits setup based on transfer type
        if transfer.k_ei == 1 && !is_moon_change  # exterior leveraging
            f1_limits = [f1, π + 2 * transfer.K * π]
            f2_limits = [π, 2 * (1 + transfer.M - transfer.K) * π + f2]
            
        elseif transfer.k_ei == -1 && !is_moon_change  # interior leveraging
            if transfer.k_io1 == -1  # inbound departure
                f1_limits = [f1, 2 * transfer.K * π]
            else  # outbound departure
                f1_limits = [f1, 2 * (transfer.K + 1) * π]
            end
            
            if transfer.k_io2 == 1
                f2_limits = [0.0, f2 + 2π * (transfer.M - transfer.K)]
            else
                f2_limits = [0.0, f2 + 2π * (transfer.M - transfer.K + 1)]
            end
            
        else  # moon change
            f1_limits = [f1, π]
            f2_limits = [π, f2 + 4π]
        end
        
        # Generate arc points
        f1_vec = range(f1_limits[1], f1_limits[2], length=points)
        f2_vec = range(f2_limits[1], f2_limits[2], length=points)
        
        r1_array = zeros(3, points)
        r2_array = zeros(3, points)
        
        # Compute positions for first arc
        for (i, f) in enumerate(f1_vec)
            Kep1_temp = [Kep1[1], Kep1[2], Kep1[3], Kep1[4], Kep1[5], f]
            r, v = Kep_To_Car(Kep1_temp...)
            r1_array[:, i] = r
        end
        
        # Compute positions for second arc
        for (i, f) in enumerate(f2_vec)
            Kep2_temp = [Kep2[1], Kep2[2], Kep2[3], Kep2[4], Kep2[5], f]
            r, v = Kep_To_Car(Kep2_temp...)
            r2_array[:, i] = r
        end
        
        # Dimension recovery
        LU1 = transfer.arc1.moon.rga_dim
        LU2 = transfer.arc2.moon.rga_dim
        
        # Re-scaling
        r1_array = r1_array .* LU1
        r2_array = r2_array .* LU2
        
        # Update max distance
        if r_a_model == Kep1[1] * (1 + Kep1[2])
            r_a_model = r_a_model * LU1
        else
            r_a_model = r_a_model * LU2
        end
        
        if r_a_model > r_a_max
            r_a_max = r_a_model
        end
        
        # Concatenate arcs
        r_array_model = hcat(r1_array, r2_array)
        
        # Store trajectory points
        append!(all_x, r_array_model[1, :])
        append!(all_y, r_array_model[2, :])
        
        # Store flyby location (start of this transfer)
        push!(fb_x, r_array_model[1, 1])
        push!(fb_y, r_array_model[2, 1])
    end
    
    # Plot trajectory
    lines!(ax, all_x, all_y, 
        linewidth=0.5, 
        color=(:gray, 0.7),
        label="Trajectory"
    )
    
    # Plot flyby locations
    scatter!(ax, fb_x, fb_y,
        marker=:circle,
        markersize=8,
        color=(:gray, 0.7),
        label="Fly-bys"
    )
    
    # Plot moon orbits
    for moon in moons
        # Check for max distance
        # if moon.rga_dim > r_a_max
        #     r_a_max = moon.rga_dim
        # end
        
        α_vec = range(-π, π, length=points)
        X_vec = moon.rga_dim .* cos.(α_vec)
        Y_vec = moon.rga_dim .* sin.(α_vec)
        
        lines!(ax, X_vec, Y_vec,
            linewidth=2,
            color=moon.color,
            label=moon.name
        )
    end
    
    # Axis formatting
    axislegend(ax, position=:rt)
    xlims!(ax, -1.2 * r_a_max, 1.2 * r_a_max)
    ylims!(ax, -1.2 * r_a_max, 1.2 * r_a_max)
    
    # Save figure
    save(filepath, fig, px_per_unit=2)
    
    return fig
end


function plot_solution_2d(sol::WeightedHGSol, test_name::String, img_folder::String)
    # recover moons
    moons = sol.prb.des_space.moons

    cont_tours_array = patch_solution(sol)
    # pdv_array = [path.p_ΔV for path in sol.paths]
    n_sol = length(cont_tours_array)

    for t = 1:n_sol

        fig_tour_name = "tour_pdv_$(t).png"
        fig_tour = plot_moon_tour_2d(
            moons,
            cont_tours_array[t],
            joinpath(img_folder, fig_tour_name),
            points=300
        )

    end


end


function plot_solution_2d(sol::MOHGSol, test_name::String, img_folder::String)
    # recover moons
    moons = sol.prb.des_space.moons

    cont_tours_array = patch_solution(sol)
    # pdv_array = [path.p_ΔV for path in sol.paths]
    n_sol = length(cont_tours_array)

    for t = 1:n_sol

        fig_tour_name = "tour_$(t).png"
        fig_tour = plot_moon_tour_2d(
            moons,
            cont_tours_array[t],
            joinpath(img_folder, fig_tour_name),
            points=300
        )

    end


end