"""
Module containing the saving and loading function for DesignSpace objects into .json and .csv files
"""

"""
    save_designspace_full(test_name::String, moons::Vector{Moon}, vinf_grid, Res_selection, Res_list, ΔV_max, only_external_vilts, choose_vilts, only_vinf_down,toll)

Creates a DesignSpace struct, a folder at `test_name`, a .json file with input info, and per-level .json/.csv files.
"""
function save_designspace_full(test_name::String, moons::Vector{Moon}, vinf_grid, Res_selection, Res_list, ΔV_max, only_external_vilts,choose_vilts, only_vinf_down,toll)
    # 1. Create DesignSpace
    ds = create_designspace(test_name, moons, vinf_grid, Res_selection, Res_list, ΔV_max, only_external_vilts, choose_vilts, only_vinf_down,toll)
    
    # 2. Create folder
    if !isdir(test_name)
        mkpath(test_name)
    end

    # 3. Save input info as .json
    input_dict = Dict(
        "moons" => [Dict("name" => moon.name, "id" => moon.id) for moon in moons],
        "vinf_grid" => collect(vinf_grid),
        "Res_selection" => Res_selection,
        "Res_list" => Res_list,
        "ΔV_max" => ΔV_max,
        "only_external_vilts" => only_external_vilts,
        "choose_vilts" => choose_vilts,
        "only_vinf_down" => only_vinf_down,
        "toll" => toll
    )
    JSON3.write(joinpath(test_name, "00_inputs.json"), input_dict; indent=4);

    # 3b. Save data info (totals and id_to_transfer mapping)
    data_dict = Dict(
        "total_bal" => ds.total_bal,
        "total_vilts" => ds.total_vilts,
        "total_moonchanges" => ds.total_moonchanges,
        "id_to_transfer" => Dict(string(k) => v for (k, v) in ds.id_to_transfer);  # keys as strings for JSON
    )
    JSON3.write(joinpath(test_name, "01_data.json"), data_dict; indent=4);

    # 4. Per-level files
    println(length(ds.vinf_levels))
    for level in ds.vinf_levels
        # println(length(ds.vinf_levels))
        # println(level.id)
        # --- Save level definition as JSON ---
        level_def = Dict(
            "id" => level.id,
            "moon" => level.moon.name,
            "moon_id" => level.moon.id,
            "v∞" => level.v∞,
            "v∞_dim" => level.v∞_dim,
            "δmax" => level.δmax,
            "res_list" => level.res_list,
            "α_list" => level.α_list
        )
        JSON3.write(joinpath(test_name, @sprintf("%d-def.json", level.id)), level_def; indent=4);

        # --- Helper to flatten a KepTransfer for CSV (without enc_par and coe) ---
        function flatten_keptransfer_csv(tr)
            Dict(
                "id" => tr.id,
                "N" => tr.N,
                "M" => tr.M,
                "K" => tr.K,
                "k_ei" => tr.k_ei,
                "k_io1" => tr.k_io1,
                "k_io2" => tr.k_io2,
                "ΔV" => tr.ΔV,
                "ToF" => tr.ToF,
                "res_ratio" => tr.res_ratio,
                "name" => tr.name,
                "arc1_moon_id" => tr.arc1.moon.id,
                "arc1_v∞" => tr.arc1.v∞,
                "arc1_α" => tr.arc1.α,
                "arc1_κ" => tr.arc1.κ,
                "arc1_L" => tr.arc1.L,
                "arc1_sma" => tr.arc1.sma,
                "arc1_ecc" => tr.arc1.ecc,
                "arc1_inc" => tr.arc1.inc,
                "arc1_Ω" => tr.arc1.Ω,
                "arc1_ω" => tr.arc1.ω,
                "arc1_f" => tr.arc1.f,
                "arc2_moon_id" => tr.arc2.moon.id,
                "arc2_v∞" => tr.arc2.v∞,
                "arc2_α" => tr.arc2.α,
                "arc2_κ" => tr.arc2.κ,
                "arc2_L" => tr.arc2.L,
                "arc2_sma" => tr.arc2.sma,
                "arc2_ecc" => tr.arc2.ecc,
                "arc2_inc" => tr.arc2.inc,
                "arc2_Ω" => tr.arc2.Ω,
                "arc2_ω" => tr.arc2.ω,
                "arc2_f" => tr.arc2.f
            )
        end

        # --- Save transfer lists as CSV ---
        transfer_types = [
            ("res_transfers", level.res_transfers),
            ("nonres_transfers_I", level.nonres_transfers_I),
            ("nonres_transfers_O", level.nonres_transfers_O),
            ("vilts_outgoing", level.vilts_outgoing),
            ("vilts_incoming", level.vilts_incoming)
        ]
        for (name, transfers) in transfer_types
            if !isempty(transfers)
                rows = [flatten_keptransfer_csv(tr) for tr in transfers];
                file = joinpath(test_name, @sprintf("%d-%s.csv", level.id, name));
                CSV.write(file, Tables.columntable(rows))
            end
        end

        # --- Save moon_change_set as CSV ---
        moon_change_file = joinpath(test_name, @sprintf("%d-moon_change_set.csv", level.id))
        if !isempty(level.moon_change_set)
            rows = [Dict("entry1" => mc[1], "entry2" => mc[2], "entry3" => mc[3]) for mc in level.moon_change_set];
            CSV.write(moon_change_file, Tables.columntable(rows))
        else
            # Write empty CSV with headers
            CSV.write(moon_change_file, Tables.columntable([Dict("entry1"=>missing, "entry2"=>missing, "entry3"=>missing)]));
        end
    end
    
    return ds
end



"""
    load_designspace_full(folder_path::String, moons::Vector{Moon}) -> DesignSpace

Loads a DesignSpace structure from a folder containing inputs.json and per-level .json/.csv files.
You must provide the moons vector, as only moon names/IDs are stored in the files.
"""
function load_designspace_full(folder_path::String, moons::Vector{Moon})
    # Helper to get moon object by name or id
    moon_by_name = Dict(moon.name => moon for moon in moons)
    moon_by_id = Dict(moon.id => moon for moon in moons)

    # 1. Load input info from inputs.json
    input_dict = JSON3.read(joinpath(folder_path, "00_inputs.json"))
    
    # Reconstruct vinf_grid
    vinf_grid_arr = Array(input_dict["vinf_grid"])
    
    if length(vinf_grid_arr) > 3
        N = Int(length(vinf_grid_arr)/3)
        
        vinf_grid = SMatrix{N, 3, Float64}(vinf_grid_arr)
    elseif length(vinf_grid_arr) == 3
        N = 1
        vinf_grid = SMatrix{1, 3, Float64}(reshape(vinf_grid_arr, 1, 3))
    else
        error("vinf_grid_arr has unexpected dimensions")
    end

    # 2. Find all level definition files
    level_files = filter(f -> endswith(f, "-def.json"), readdir(folder_path))
    level_ids = [parse(Int, split(f, "-")[1]) for f in level_files]
    sort!(level_ids)

    # 3. Load each VinfLevel
    vinf_levels = VinfLevel[]
    for level_id in level_ids
        # Load level definition
        level_def = JSON3.read(joinpath(folder_path, @sprintf("%d-def.json", level_id)))
        
        moon = moon_by_id[level_def["moon_id"]]
        
        # Helper to load KepTransfer from CSV
        function load_transfers_csv(filename)
            if !isfile(filename)
                return KepTransfer[]
            end
            df = CSV.read(filename, DataFrame)
            if nrow(df) == 0
                return KepTransfer[]
            end
            
            transfers = KepTransfer[]
            for row in eachrow(df)
                # Reconstruct enc_par and coe from saved parameters
                arc1_enc_par = [row.arc1_v∞, row.arc1_α, row.arc1_κ, row.arc1_L]
                arc1_coe = [row.arc1_sma, row.arc1_ecc, row.arc1_inc, row.arc1_Ω, row.arc1_ω, row.arc1_f]
                arc2_enc_par = [row.arc2_v∞, row.arc2_α, row.arc2_κ, row.arc2_L]
                arc2_coe = [row.arc2_sma, row.arc2_ecc, row.arc2_inc, row.arc2_Ω, row.arc2_ω, row.arc2_f]
                
                arc1 = KepArc(
                    moon_by_id[row.arc1_moon_id],
                    row.arc1_v∞,
                    row.arc1_α,
                    row.arc1_κ,
                    row.arc1_L,
                    arc1_enc_par,
                    row.arc1_sma,
                    row.arc1_ecc,
                    row.arc1_inc,
                    row.arc1_Ω,
                    row.arc1_ω,
                    row.arc1_f,
                    arc1_coe
                )
                arc2 = KepArc(
                    moon_by_id[row.arc2_moon_id],
                    row.arc2_v∞,
                    row.arc2_α,
                    row.arc2_κ,
                    row.arc2_L,
                    arc2_enc_par,
                    row.arc2_sma,
                    row.arc2_ecc,
                    row.arc2_inc,
                    row.arc2_Ω,
                    row.arc2_ω,
                    row.arc2_f,
                    arc2_coe
                )
                tr = KepTransfer(
                    row.id,
                    row.N,
                    row.M,
                    row.K,
                    row.k_ei,
                    row.k_io1,
                    row.k_io2,
                    arc1,
                    arc2,
                    row.ΔV,
                    row.ToF,
                    row.res_ratio,
                    row.name
                )
                push!(transfers, tr)
            end
            return transfers
        end

        # Load transfer lists
        res_transfers = load_transfers_csv(joinpath(folder_path, @sprintf("%d-res_transfers.csv", level_id)))
        nonres_transfers_I = load_transfers_csv(joinpath(folder_path, @sprintf("%d-nonres_transfers_I.csv", level_id)))
        nonres_transfers_O = load_transfers_csv(joinpath(folder_path, @sprintf("%d-nonres_transfers_O.csv", level_id)))
        vilts_outgoing = load_transfers_csv(joinpath(folder_path, @sprintf("%d-vilts_outgoing.csv", level_id)))
        vilts_incoming = load_transfers_csv(joinpath(folder_path, @sprintf("%d-vilts_incoming.csv", level_id)))

        # Load moon_change_set
        mc_file = joinpath(folder_path, @sprintf("%d-moon_change_set.csv", level_id))
        moon_change_set = Vector{Any}()
        if isfile(mc_file)
            df_mc = CSV.read(mc_file, DataFrame)
            for row in eachrow(df_mc)
                if !ismissing(row.entry1)
                    push!(moon_change_set, [row.entry1, row.entry2, row.entry3])
                end
            end
        end

        # Create VinfLevel
        level = VinfLevel(
            level_def["id"],
            moon,
            level_def["v∞"],
            level_def["v∞_dim"],
            level_def["δmax"],
            level_def["res_list"],
            level_def["α_list"],
            res_transfers,
            nonres_transfers_I,
            nonres_transfers_O,
            vilts_outgoing,
            vilts_incoming,
            moon_change_set
        )
        push!(vinf_levels, level)
    end

    # 4. Load data info (totals and id_to_transfer)
    data_dict = JSON3.read(joinpath(folder_path, "01_data.json"))
    total_bal = data_dict["total_bal"]
    total_vilts = data_dict["total_vilts"]
    total_moonchanges = data_dict["total_moonchanges"]
    
    # Rebuild id_to_transfer dictionary
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

    # 5. Create DesignSpace
    ds = DesignSpace(
        basename(folder_path),
        moons,
        vinf_grid,
        vinf_levels,
        input_dict["Res_selection"],
        input_dict["Res_list"],
        input_dict["ΔV_max"],
        input_dict["only_external_vilts"],
        input_dict["choose_vilts"],
        input_dict["only_vinf_down"],
        input_dict["toll"],
        total_bal,
        total_vilts,
        total_moonchanges,
        id_to_transfer
    )

    return ds
end