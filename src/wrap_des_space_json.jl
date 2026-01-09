"""
Module containing the functions to write a des space into a json and then to create a des_space from it
"""

"""
    save_designspace_json(ds::DesignSpace, filename::String)

Stores a DesignSpace structure into a JSON file, including all its parameters and all VinfLevel objects and their transfers.
"""
function save_designspace_json(ds::DesignSpace, filename::String)

    # Serialize the DesignSpace
    ds_dict = Dict(
        "id" => ds.id,
        "moons" => [moon.name for moon in ds.moons],
        "vinf_grid" => ds.vinf_grid |> collect,
        "vinf_levels" => [serialize_vinflevel(level) for level in ds.vinf_levels],
        "Res_selection" => ds.Res_selection,
        "Res_list" => ds.Res_list,
        "ΔV_max" => ds.ΔV_max,
        "only_external_vilts" => ds.only_external_vilts,
        "choose_vilts" => ds.choose_vilts,
        "toll" => ds.toll,
        "total_bal" => ds.total_bal,
        "total_vilts" => ds.total_vilts,
        "total_moonchanges" => ds.total_moonchanges
    )

    open(filename, "w") do io
        JSON3.write(io, ds_dict; indent = 4)
    end
end

# Helper to serialize a KepArc
function serialize_keparc(arc)
    return Dict(
        "moon" => arc.moon.name,
        "v∞" => arc.v∞,
        "α" => arc.α,
        "κ" => arc.κ,
        "L" => arc.L,
        "enc_par" => arc.enc_par,
        "sma" => arc.sma,
        "ecc" => arc.ecc,
        "inc" => arc.inc,
        "Ω" => arc.Ω,
        "ω" => arc.ω,
        "f" => arc.f,
        "coe" => arc.coe
    )
end

# Helper to serialize a KepTransfer
function serialize_keptransfer(tr)
    return Dict(
        "id" => tr.id,
        "N" => tr.N,
        "M" => tr.M,
        "K" => tr.K,
        "k_ei" => tr.k_ei,
        "k_io1" => tr.k_io1,
        "k_io2" => tr.k_io2,
        "arc1" => serialize_keparc(tr.arc1),
        "arc2" => serialize_keparc(tr.arc2),
        "ΔV" => tr.ΔV,
        "ToF" => tr.ToF,
        "res_ratio" => tr.res_ratio,
        "name" => tr.name
    )
end


# Helper to serialize a VinfLevel
function serialize_vinflevel(level)
    return Dict(
        "id" => level.id,
        "moon" => level.moon.name,
        "v∞" => level.v∞,
        "v∞_dim" => level.v∞_dim,
        "δmax" => level.δmax,
        "res_list" => level.res_list,
        "α_list" => level.α_list,
        "res_transfers" => [serialize_keptransfer(tr) for tr in level.res_transfers],
        "nonres_transfers_I" => [serialize_keptransfer(tr) for tr in level.nonres_transfers_I],
        "nonres_transfers_O" => [serialize_keptransfer(tr) for tr in level.nonres_transfers_O],
        "vilts_outgoing" => [serialize_keptransfer(tr) for tr in level.vilts_outgoing],
        "vilts_incoming" => [serialize_keptransfer(tr) for tr in level.vilts_incoming],
        "moon_change_set" => level.moon_change_set
    )
end



"""
    load_designspace_json(filename::String, moons::Vector{Moon}) -> DesignSpace

Loads a DesignSpace structure from a JSON file, reconstructing all VinfLevel and KepTransfer objects.
You must provide the moons vector, as only moon names are stored in the JSON.
"""
function load_designspace_json(filename::String, moons::Vector{Moon})
    # Helper to get moon object by name
    moon_by_name = Dict(moon.name => moon for moon in moons)

    # Load JSON
    ds_dict = JSON3.read(filename)

    # Reconstruct DesignSpace
    vinf_levels = [deserialize_vinflevel(level, moon_by_name) for level in ds_dict["vinf_levels"]]
    vinf_grid_arr = Array(ds_dict["vinf_grid"])  # convert JSON3.Array to Array

    if ndims(vinf_grid_arr) == 2
        N = size(vinf_grid_arr, 1)
        vinf_grid = SMatrix{N, 3, Float64}(vinf_grid_arr)
    elseif ndims(vinf_grid_arr) == 1
        N = 1
        vinf_grid = SMatrix{1, 3, Float64}(reshape(vinf_grid_arr, 1, 3))
    end

    vinf_grid = SMatrix{N, 3, Float64}(vinf_grid_arr)
    Res_list = ds_dict["Res_list"]

    ds = DesignSpace(
        ds_dict["id"],
        moons,
        vinf_grid,
        vinf_levels,
        ds_dict["Res_selection"],
        Res_list,
        ds_dict["ΔV_max"],
        ds_dict["only_external_vilts"],
        ds_dict["choose_vilts"],
        ds_dict["toll"],
        ds_dict["total_bal"],
        ds_dict["total_vilts"],
        ds_dict["total_moonchanges"]
    )
    return ds
end


# Helper to deserialize KepArc
function deserialize_keparc(dict, moon_by_name)
    moon = moon_by_name[dict["moon"]]
    return KepArc(
        moon,
        dict["v∞"],
        dict["α"],
        dict["κ"],
        dict["L"],
        dict["enc_par"],
        dict["sma"],
        dict["ecc"],
        dict["inc"],
        dict["Ω"],
        dict["ω"],
        dict["f"],
        dict["coe"]
    )
end


# Helper to deserialize KepTransfer
function deserialize_keptransfer(dict, moon_by_name)
    return KepTransfer(
        dict["id"],
        dict["N"],
        dict["M"],
        dict["K"],
        dict["k_ei"],
        dict["k_io1"],
        dict["k_io2"],
        deserialize_keparc(dict["arc1"], moon_by_name),
        deserialize_keparc(dict["arc2"], moon_by_name),
        dict["ΔV"],
        dict["ToF"],
        dict["res_ratio"],
        dict["name"]
    )
end

# Helper to deserialize VinfLevel
function deserialize_vinflevel(dict, moon_by_name)
    moon = moon_by_name[dict["moon"]]
    level = VinfLevel(
        dict["id"],
        moon,
        dict["v∞"],
        dict["v∞_dim"],
        dict["δmax"],
        dict["res_list"],
        dict["α_list"],
        [deserialize_keptransfer(tr, moon_by_name) for tr in dict["res_transfers"]],
        [deserialize_keptransfer(tr, moon_by_name) for tr in dict["nonres_transfers_I"]],
        [deserialize_keptransfer(tr, moon_by_name) for tr in dict["nonres_transfers_O"]],
        [deserialize_keptransfer(tr, moon_by_name) for tr in dict["vilts_outgoing"]],
        [deserialize_keptransfer(tr, moon_by_name) for tr in dict["vilts_incoming"]],
        dict["moon_change_set"]
    )
    return level
end