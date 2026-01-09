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

# define problem
prob = WeightedHGProb(dep_cond, tar_cond, des_space, p_ΔV_array, flag_in2out, pump_down_only)

# solve it
sol= WeightedHGSol(prob)
# println(sol.ΔV_array)
# println(sol.ToF_array)

# plot solutions
begin
    img_folder = working_folder*"/images/"
    fig_dvtof_name = test_name*"-dvtof_front.png"
    fig_dvtof = Figure(size=(800, 600))
    ax = Axis(fig_dvtof[1, 1],
    xlabel="ΔV [m/s]",
    ylabel="Time of Flight [days]",
    title="Pareto Front: ΔV vs ToF"
    )

    scatter!(ax, sol.ΔV_array, sol.ToF_array,
    color=p_ΔV_array,
    colormap=:viridis,
    colorrange=(0, 1),
    markersize=12
    )

    lines!(ax, sol.ΔV_array , sol.ToF_array,
    color=:gray,
    linestyle=:dash,
    linewidth=1.5
    )

    Colorbar(fig_dvtof[1, 2], 
    limits=(0, 1), 
    colormap=:viridis,
    label="p_ΔV"
    )

    save(joinpath(img_folder, fig_dvtof_name), fig_dvtof)

end

# post processing
# Plot and save all tours in planet centred inertial frame
begin
    tour_img_folder = joinpath(img_folder, test_name*"-tours")
    if !isdir(tour_img_folder)
        mkpath(tour_img_folder)
    end
    plot_solution_2d(sol, test_name, tour_img_folder)
end

# select a solution and study it
# compute flyby altitudes
path_id = 97  # select a solution 
path_seq = sol.paths[path_id].seq
tour_seq, fb_alts = compute_tour_fbs(dep_cond, tar_cond, path_seq, des_space.vinf_levels)

leg_data = compute_leg_data(tour_seq, des_space.vinf_levels)
for (moon_name, data) in leg_data
    println("Moon: $(moon_name)")
    println("  ΔV: $(data[:ΔV]) m/s")
    println("  ToF: $(data[:ToF]) days")
    println("  Flybys: $(data[:fb_altitudes]) km")
    println("  v∞: $(data[:vinf_values]) km/s")
end

cont_tour, dv_check, tof_check = patch_tour(tour_seq)
plot_tour_circles(dep_cond, tar_cond, cont_tour, des_space.vinf_levels, tour_img_folder)


