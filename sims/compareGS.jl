include(joinpath("..", "universe.jl"))
# Define des space
begin
    working_folder = "sims/"
    test_name = "TetEnc"
    test_path = working_folder*test_name

    moons = [Tethys, Enceladus];

    v∞_step = 50.0

    v∞_grid = @SMatrix[ 
            0.5 0.9 v∞_step;
            0.3 0.8 v∞_step];  # v∞ min km/s,v∞ max km/s, step m/s

    ΔV_max =  50.0;  # m/s
    
    graph_node_no = 0
    graph_edge_no = 0
    num_solutions_emoa = 0
    
    only_external_vilts = false;
    choose_vilts = false;
    only_vinf_down = true;
    Res_selection = "Manual";
    Res_list = [[[1,1], [11,9], [6,5], [7,6], [8,7], [9,8], [10,9], [11,10], [12,11], [13,12], [14,13], [15,14], [19,18], [25,24], [35,34], [34,35], [24,25], [18,19], [14,15], [13,14], [10,11], [9,10], [7,8], [13,15], [34,33], [33,32], [33,32], [32,33], [32,31], [31,32], [31,30], [30,31], [30,29], [29,30], [29,28], [28,29], [28,27], [27,28]],  # Tethys
                [[1,1], [7,6], [20,17], [15,13], [8,7], [17,15], [9,8], [19,17], [10,9], [21,19], [11,10], [12,11], [13,12], [14,13], [15,14], [16,15], [19,18], [24,23], [17,16], [21,20], [13,11], [22,19], [15,13], [25,22], [18,17]]  # Enceladus
        ];  # desired list of resonances per moon, to be provided if res_selection is "custom"
    
    toll = 1e-10;
    
    # Placeholder for the EMOA* local directory - repo at https://github.com/rap-lab-org/public_emoa
    emoa_local_dir = "/home/fmaccari/Codes/Resources/public_emoa/build";
end

# define design space and problem
begin 

    moons = [Titan, Rhea, Dione, Tethys, Enceladus];

    des_space = save_designspace_full(test_path, moons, v∞_grid, Res_selection, Res_list, ΔV_max, only_external_vilts, choose_vilts, only_vinf_down,toll)

end

## define and solve problem - EMOA*
begin

    prob_emoa = MOHGProb(dep_cond, tar_cond, des_space, flag_in2out, pump_down_only)

    # create graph definition .gr files
    graph_emoa_filename = string("def")
    hop_graph_emoa = prepare_MOHG_search(prob_emoa, test_path, graph_emoa_filename)
    
    ### Setup paths and parameters for the EMOA* call
    dv_graph_path = joinpath(test_path, "$(graph_emoa_filename)_dv.gr")
    tof_graph_path = joinpath(test_path, "$(graph_emoa_filename)_tof.gr")
    res_path = joinpath(test_path, "EMOA_results_$(k).txt")
    
    # Assuming start and target nodes are extracted from prob_emoa (adjust as needed)
    start_node = 1 
    target_node = 2 
    max_tcomp = 600 
    num_objs = 2
    ## solve problem - EMOA*
    # Create the command.
    emoa_cmd = Cmd(`./run_emoa $start_node $target_node $max_tcomp $num_objs $(abspath(dv_graph_path)) $(abspath(tof_graph_path)) $(abspath(res_path))`, dir=emoa_local_dir)
    
    ### Execute and capture output
    terminal_output = read(ignorestatus(emoa_cmd), String)
    ### Parse the execution time
    # This targets: "with [NUM] solutions in [NUM](for heu) + [NUM](for search)"
    time_match = match(r"with\s+(\d+)\s+solutions\s+in\s+([\d.]+)\s*\(for heu\)\s*\+\s*([\d.]+)\s*\(for search\)", terminal_output)
    if time_match !== nothing
        # Parse captures in order of appearance
        num_sols = parse(Int, time_match.captures[1])
        time_heu = parse(Float64, time_match.captures[2])
        time_search = parse(Float64, time_match.captures[3])
        
        num_solutions_emoa = num_sols
        comp_times_emoa = time_heu + time_search
    else
        println("Warning: Could not parse execution time for iteration.")
        println("Raw terminal output was: \n", terminal_output)
    end

    # after having solved the problem, post process graph solutions
    sol_emoa = parse_EMOA_results(res_path, hop_graph_emoa, prob_emoa)
end
## define problem - wgt Dijkstra
begin 
    prob_wgt = WeightedHGProb(dep_cond, tar_cond, des_space, p_ΔV_array, flag_in2out, pump_down_only)

    # solve it
    sol_wgt= WeightedHGSol(prob_wgt)
end

# plot solutions
begin
    img_folder = working_folder*"/images/"
    fig_dvtof_GS_name = test_name*"-dvtof_fronts-$(length(p_ΔV_array)).png"
    fig_dvtof_GS = Figure(size=(600, 500))
    ax = Axis(fig_dvtof_GS[1, 1],
    xlabel="ΔV [m/s]",
    ylabel="Time of Flight [days]",
    # title="Pareto Front: ΔV vs ToF"
    )

    cmap = :managua
    scatter!(ax, sol_wgt.ΔV_array, sol_wgt.ToF_array,
    color=p_ΔV_array,
    colormap=cmap,
    colorrange=(0, 1),
    markersize=12,
    marker = :diamond,
    label = "Weighted Dijkstra"
    )

    lines!(ax, sol_wgt.ΔV_array , sol_wgt.ToF_array,
    color=:gray,
    linestyle=:dash,
    linewidth=1.5
    )

    Colorbar(fig_dvtof_GS[1, 2], 
    limits=(0, 1), 
    colormap=cmap,
    label="p_ΔV"
    )

    scatter!(ax, sol_emoa.ΔV_array, sol_emoa.ToF_array,
    color=:black,
    markersize=6,
    marker=:dtriangle,
    label = "EMOA*"
    )

    lines!(ax, sol_emoa.ΔV_array , sol_emoa.ToF_array,
    color=:gray,
    linestyle=:dash,
    linewidth=1.5
    )

    save(joinpath(img_folder, fig_dvtof_GS_name), fig_dvtof_GS)

end

# post processing
# select a solution and study it
sol = sol_wgt;
# Plot and save all tours in planet centred inertial frame
begin

    if typeof(sol) == MOHGSol
        tour_img_folder = joinpath(img_folder, test_name*"-emoa-tours")
    else
        tour_img_folder = joinpath(img_folder, test_name*"-wgt-tours")
    end
        if !isdir(tour_img_folder)
        mkpath(tour_img_folder)
    end
    plot_solution_2d(sol, test_name, tour_img_folder)
end


# compute flyby altitudes
path_id = 90  # select solution index
path_seq = sol.paths[path_id].seq;
tour_seq, fb_alts = compute_tour_fbs(dep_cond, tar_cond, path_seq, des_space.vinf_levels);

leg_data = compute_leg_data(tour_seq, des_space.vinf_levels)
for (moon_name, data) in leg_data
    
    println("Moon: $(moon_name)")
    println("  ΔV: $(data[:ΔV]) m/s")
    println("  ToF: $(data[:ToF]) days")
    println("  Flybys: $(data[:fb_altitudes]) km")
    println("  v∞: $(data[:vinf_values]) km/s")
end

cont_tour, dv_check, tof_check = patch_tour(tour_seq);
plot_tour_circles(dep_cond, tar_cond, cont_tour, des_space.vinf_levels, tour_img_folder)


# Print detailed information for each leg in the continuous tour
println("\n=== Tour Leg Details ===")
println("Total ΔV: $dv_check m/s")
println("Total ToF: $tof_check days")
println("\nLeg-by-leg breakdown:\n")

# First flyby: from departure condition to first transfer
first_leg = cont_tour[1]
level_dep = des_space.vinf_levels[Int(dep_cond[1])]
moon_dep = level_dep.moon
v∞_dim_dep = level_dep.v∞_dim
α_dep_cond = dep_cond[2]
α_first = first_leg.arc1.α

if moon_dep.id == first_leg.arc1.moon.id
    FB_alt_dep = compute_fb_altitude(v∞_dim_dep, α_dep_cond, α_first, moon_dep)
    println("Initial flyby (dep → Leg 1): $(round(FB_alt_dep, digits=2)) km")
else
    println("Initial: Moon change from $(moon_dep.name) to $(first_leg.arc1.moon.name)")
end
println()

for (n, leg) in enumerate(cont_tour)
    println("Leg $n: $(leg.name)")
    println("  ΔV: $(round(leg.ΔV, digits=3)) m/s, ToF: $(round(leg.ToF, digits=2)) days")
    println("  $(leg.arc1.moon.name) (v∞ = $(round(leg.arc1.v∞ * leg.arc1.moon.vga_dim, digits=3)) km/s) → $(leg.arc2.moon.name) (v∞ = $(round(leg.arc2.v∞ * leg.arc2.moon.vga_dim, digits=3)) km/s)")
    
    # Flyby altitude to next leg
    if n < length(cont_tour)
        next_leg = cont_tour[n + 1]
        FB_alt, mc_flag = compute_fb_altitude(leg, next_leg)
        
        if mc_flag
            println("  → Moon change to $(next_leg.arc1.moon.name)")
        else
            println("  → Flyby to Leg $(n+1): $(round(FB_alt, digits=2)) km")
        end
    else
        # Last leg - check if there's a flyby to target condition
        level_tar = des_space.vinf_levels[Int(tar_cond[1])]
        moon_tar = level_tar.moon
        α_tar_cond = tar_cond[2]
        
        if α_tar_cond != -1.0 && moon_tar.id == leg.arc2.moon.id
            v∞_dim_tar = leg.arc2.v∞ * leg.arc2.moon.vga_dim
            α_last = leg.arc2.α
            FB_alt_tar = compute_fb_altitude(v∞_dim_tar, α_last, α_tar_cond, moon_tar)
            println("  → Final flyby (Leg $n → target): $(round(FB_alt_tar, digits=2)) km")
        end
    end
    println()
end