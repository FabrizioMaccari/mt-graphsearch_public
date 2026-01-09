"""
Module containting the base objects for the keplerian pathced conics model with circular planar moons orbits.
The main objects are:
- Moon: object containing the physical and orbital parameters of a moon
- KepArc: object containing the vinf parameters and the orbital parameters of a keplerian arc intersecting the orbit of a moon
- KepTransfer: object containing the parameters of a transfer between two encounters, specified in its two KepArc objects
"""
struct Moon 
    name ::String # string containing the name of the Moon
    id ::Integer

    # dimensional quantities
    μ_dim ::Float64  # gravitational parameter of the moon
    radius_dim ::Float64 # radius of the moon
    rga_dim ::Float64  # distance of the Moon from the primary
    vga_dim ::Float64  # velocity of the  Moon around the primary
    Tga_dim ::Float64  # orbital period of the Moon
    min_altfb_dim ::Float64  # minimum altitude of fly by
    min_rp_fb_dim ::Float64  # minimum pericenter radius of the fly by
    μP_dim ::Float64  # gravitational constant of the planet

    levels :: Vector{Any}  # empty list, to store the v_inf_level objects 
    vinf_values :: Vector{Float64}  # empty array to store the corresponding vinf values
    color ::String # Color for plotting

    # scaled quantities
    μ::Float64  # scaled gravitaional parameter
    radius ::Float64 # radius of the moon
    rga ::Float64  # distance of the Moon from the primary
    vga ::Float64  # velocity of the  Moon around the primary
    Tga ::Float64  # orbital period of the Moon
end

function init_moon(name::String, id::Integer, μ_dim::Float64, radius_dim::Float64, rga_dim::Float64, vga_dim::Float64, Tga_dim::Float64, min_altfb_dim::Float64, μP_dim::Float64, color::String)
    min_rpfb_dim = min_altfb_dim + radius_dim
    # non dimensional attributes
    VU = vga_dim
    TU = Tga_dim
    LU = rga_dim
    μ = (μ_dim/LU^3)*TU^2
    radius = radius_dim/LU
    rga = rga_dim/LU
    vga = vga_dim/VU
    Tga = Tga_dim/TU
    levels = Vector{Any}()  # empty list, to store the v_inf_level objects 
    vinf_values = Vector{Float64}()  # empty array to store the corresponding vinf values
    return Moon(name, id, μ_dim, radius_dim, rga_dim, vga_dim, Tga_dim, min_altfb_dim, min_rpfb_dim, μP_dim, levels, vinf_values, color, μ, radius, rga, vga, Tga)
end

struct KepArc 
    moon::Moon
    # encounter parameters
    v∞::Float64     # non-dim (vinf)
    α::Float64      # rad (alpha)
    κ::Float64      # rad (kappa)
    L::Float64      # rad
    enc_par::Vector{Float64}
    # orbital elements
    sma::Float64      # non-dim (sma)
    ecc::Float64      # non-dim (ecc)
    inc::Float64      # rad (inc)
    Ω::Float64      # rad
    ω::Float64      # rad
    f::Float64      # rad
    coe::Vector{Float64}
end

function keparc_vinf(moon::Moon, enc_par::Vector{Float64})
    # recover encounter parameter
    v∞ = enc_par[1]  # scaled, angles in radians
    α = enc_par[2]
    κ = enc_par[3]
    L = enc_par[4]
    # convert to coe and recover single elements
    coe = enc_to_coe(enc_par) # scaled sma, angles in radians
    sma = coe[1]  
    ecc = coe[2]
    inc = coe[3]
    Ω = coe[4]
    ω = coe[5]
    f = coe[6]
    return KepArc(moon, v∞, α, κ, L, enc_par, sma, ecc, inc, Ω, ω, f, coe)
end


# define the model object type
struct KepTransfer 
    id::Integer
    N::Integer 
    M::Integer
    K::Integer
    k_ei::Integer
    k_io1::Integer
    k_io2::Integer
    arc1::KepArc
    arc2::KepArc
    ΔV::Float64   # m/s
    ToF::Float64  # days
    res_ratio::Float64
    name::String
end

# define the constructors for the possible KepTransfer types
function vilt(id::Integer, N::Integer, M::Integer, K::Integer, k_ei::Integer, k_io1::Integer, k_io2::Integer, arc1::KepArc, arc2::KepArc, ΔV::Float64, ToF::Float64)
    res_ratio = N / M
    ei = k_ei == 1 ? "ext-" : "int-"
    dep_io = k_io1 == 1 ? "O" : "I"
    arr_io = k_io2 == 1 ? "O" : "I"
    name = string(arc1.moon.name, "-", ei, dep_io, arr_io, "-", string(N), ":", string(M), "(", string(K), ")")
    return KepTransfer(id, N, M, K, k_ei, k_io1, k_io2, arc1, arc2, ΔV, ToF, res_ratio, name)
end

function res(id::Integer, N::Integer, M::Integer, k_io::Integer, arc::KepArc)
    res_ratio = N / M
    k_ei = 0
    k_io1 = k_io
    k_io2 = k_io
    arc1 = arc
    arc2 = arc
    ΔV = 0.0
    ToF = N * arc1.moon.Tga_dim / 86400
    dep_io = k_io == 1 ? "O" : "I"
    arr_io = k_io == 1 ? "O" : "I"
    name = string(arc1.moon.name, "-", dep_io, arr_io, "-", string(N), ":", string(M))
    return KepTransfer(id, N, M, 0, k_ei, k_io1, k_io2, arc1, arc2, ΔV, ToF, res_ratio, name)
end

function nonres(id::Integer, N::Integer, M::Integer, k_io1::Integer, arc1::KepArc, arc2::KepArc, ToF::Float64)
    res_ratio = N / M
    k_ei = 0
    k_io2 = -k_io1
    ΔV = 0.0
    dep_io = k_io1 == 1 ? "O" : "I"
    arr_io = k_io1 == 1 ? "I" : "O"
    name = string(arc1.moon.name, "-", dep_io, arr_io, "-", string(N), ":", string(M))
    return KepTransfer(id, N, M, 0, k_ei, k_io1, k_io2, arc1, arc2, ΔV, ToF, res_ratio, name)
end

# define the constructors for the possible KepTransfer types
function moonchange(id::Integer, k_io1::Integer, k_io2::Integer, arc1::KepArc, arc2::KepArc, ToF::Float64)
    name = string(arc1.moon.name, " to ", arc2.moon.name)
    return KepTransfer(id, 1, 1, 0, 1, k_io1, k_io2, arc1, arc2, 0.0, ToF, 1.0, name)
end

# defines a model data structure for failed zero findings
function keptransfer_fail(moon::Moon)
    arc_fail = KepArc(moon, 0.0, 0.0, 0.0, 0.0, [0.0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, [0.0])
    return KepTransfer(0,1, 1, 0, 1, 1, 1, arc_fail, arc_fail, 100000.0, 0.0, 1.0, "Fail")
end