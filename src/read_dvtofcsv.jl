
function read_dvtof_csv(filepath::String)
    if !isfile(filepath)
        error("file not found")
        end
    df = CSV.read(filepath, DataFrame)
    dv_array = df[:,1]
    tof_array = df[:,2]
    return dv_array, tof_array
end