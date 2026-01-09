include(joinpath("..", "universe.jl"))

# Define des space
begin
    working_folder = "sims/"
    test_name = "enceladus_vinfdown"
    create_des = true  # true to create the design space, false to load it
    test_path = working_folder*test_name
end
if create_des  # create a design space with given characteristics

    # specify moons to include
    moons = [Titan, Rhea, Dione, Tethys, Enceladus];
    # specify v infinite grid - 1 row per moon 
    v∞_grid = @SMatrix[1.3 1.5 100; 
                0.75 1.81 50; 
                0.65 1.2 50; 
                0.55 1.00 50;
                0.20 0.8 50];  # v∞ min km/s,v∞ max km/s, step m/s
    ΔV_max =  50.0;  # maximum m/s

    # specify design space settings
    only_external_vilts = false;
    choose_vilts = false;
    only_vinf_down = true;
    Res_selection = "Manual";
    Res_list = [[[3,1], [2,1], [1,1]],  # Titan
                [[1,1], [13,7], [9,7], [4,5], [5,3], [3,1], [7,6], [3,2], [7,5], [4,3], [5,4], [6,5], [15,14], [14,15], [6,7], [2,1], [5,3], [8,5], [7,4], [15,8], [17,8], [9,8], [8,9], [8,7], [7,8], [10,9], [9,10], [11,10], [10,11], [12,11], [11,12], [13,12], [12,13]],  # Rhea
                [[1,1], [4,3], [9,7], [5,4], [6,5], [7,6], [8,7], [9,8], [10,9], [11,10], [12,11], [13,12], [19,18], [18,19], [14,15], [12,13], [11,12], [10,11], [9,10], [8,9], [7,8], [6,7], [27,26], [26,27], [26,25], [25,26], [25,24], [24,25]],  # Dione
                [[1,1], [11,9], [6,5], [7,6], [8,7], [9,8], [10,9], [11,10], [12,11], [13,12], [14,13], [15,14], [19,18], [25,24], [35,34], [34,35], [24,25], [18,19], [14,15], [13,14], [10,11], [9,10], [7,8], [13,15], [34,33], [33,32], [33,32], [32,33], [32,31], [31,32], [31,30], [30,31], [30,29], [29,30], [29,28], [28,29], [28,27], [27,28]],  # Tethys
                [[1,1], [7,6], [20,17], [15,13], [8,7], [17,15], [9,8], [19,17], [10,9], [21,19], [11,10], [12,11], [13,12], [14,13], [15,14], [16,15], [19,18], [24,23], [17,16], [21,20], [13,11], [22,19], [15,13], [25,22], [18,17]]  # Enceladus
        ];  # desired list of resonances per moon, to be provided if res_selection is "Manual"
    toll = 1e-10;

    # create and saves design space
    des_space = save_designspace_full(test_path, moons, v∞_grid, Res_selection, Res_list, ΔV_max, only_external_vilts, choose_vilts, only_vinf_down,toll)

else  # load a pre-computed design space
    moons = [Titan, Rhea, Dione, Tethys, Enceladus];
    des_space = load_designspace_full(test_path, moons)
end

# Define problem
begin
    titan_vinf = 1.5
    dep_level_id = [level.id for level in des_space.vinf_levels if (level.v∞_dim == titan_vinf) && level.moon.id == 1] # Titan max vinf level
    dep_alpha = deg2rad(50) # initial pump angle 
    dep_cond = @SVector[dep_level_id[1], dep_alpha]
    
    tar_level_id = [level.id for level in des_space.vinf_levels if level.v∞_dim == 0.21]  # enceladus min vinf level
    
    tar_alpha = -1.0 # nan to just get to the level
    tar_cond = @SVector[tar_level_id[1], tar_alpha]

    p_ΔV_array = collect(0.0:0.0025:0.99) # vector of weights to test. Each one must be between 0 and 1: 0- time optimal, 1-fuel optimal 
    
    # select graph edges generation preferences
    flag_in2out = false  # no in to out moon change
    pump_down_only = true  # only pump down flybys
end


## define problem - EMOA*
begin
    prob_emoa = MOHGProb(dep_cond, tar_cond, des_space, flag_in2out, pump_down_only)

    # create graph definition .gr files
    hop_graph_emoa = prepare_MOHG_search(prob_emoa, test_path)

end

## solve problem - EMOA*
# create call to use in emoa_public - https://github.com/rap-lab-org/public_emoa
# command line to call: ./run_emoa  1 2 600 2 ../graph_dv.gr ../graph_tof.gr ../mttest/EMOA_results.txt
# after having solved the problem, post process graph solutions
sol_emoa = parse_EMOA_results(test_path*"/EMOA_results.txt", hop_graph_emoa, prob_emoa)

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

    # scatter!(ax, sol_emoa.ΔV_array, sol_emoa.ToF_array,
    # color=:black,
    # markersize=6,
    # marker=:dtriangle,
    # label = "EMOA*"
    # )

    # lines!(ax, sol_emoa.ΔV_array , sol_emoa.ToF_array,
    # color=:gray,
    # linestyle=:dash,
    # linewidth=1.5
    # )


    # comparison with literature - code on https://github.com/andreabellome/saturn_moon_tours
    ΔV_isae, ToF_isae = read_dvtof_csv(working_folder*"bellome_pareto_dvtof.csv")
    scatter!(ax, ΔV_isae, ToF_isae,
    color=:gray,
    markersize=6,
    marker=:utriangle,
    label = "Bellome et al."
    )
    axislegend(ax, position=:rt)

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