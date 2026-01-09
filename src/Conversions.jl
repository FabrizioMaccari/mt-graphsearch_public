"""
Module for conversion between different sets of orbital parameters and vinf related quantities.
"""


"""
    enc_to_coe(enc_par::Vector{Float64})

computes the primary-centred Keplerian classical orbital elements for an encounter with a minor body, given the
encounter paramters: v∞, α,κ, L.
The minor body is assumed to have a circular planar orbit
"""

function enc_to_coe(enc_par::Vector{Float64})
    v∞ = enc_par[1]
    α = enc_par[2]
    κ = enc_par[3]
    L = enc_par[4]

    rga = 1
    vga = 1

    # Keplerian parameters computation
    sma = rga / ( 2 - 1 / vga^2 * (v∞^2 + vga^2 + 2*v∞*vga*cos(α)) )  # semi major axis

    # inclination:
    inc = atan(sin(κ)*(sin(α)/(vga/v∞ + cos(α))))

    # correction for numerical error
    if abs(inc) < 1e-8
        inc = 0
    end

    # correction for quadrant
    if -v∞*cos(α)<= vga
        # prograde orbit
        inc = abs(inc)*sign(κ)
    else
        # retrograde
        inc = - abs(inc)*sign(κ)
    end
        # eccentricity
    ecc = sqrt( 1 - rga/sma * ( ( 3 - (v∞/vga)^2 - rga/sma )/( 2*cos(inc) ) )^2 )

    vsc = sqrt( v∞^2 + vga^2 + 2*v∞*vga*cos(α) )

    # true anomaly of the encounter
    p = sma*(1-ecc^2)
    cos_f = 1/ecc*(p -1)

    if abs(cos_f - 1) < 1e-8 # numerical correction for tangent encounters
        cos_f = 1
    end
    if abs(cos_f + 1) < 1e-8
        cos_f = - 1
    end
    f = acos(cos_f)

    
    # sign correction with κ
    if f!= 0 && f!= π 
        f = - f*sign(cos(κ))
    end

    # longitude of the ascending node - careful for κ = +- π
    sign_kappa = sign(κ)
    if κ == 0
        sign_kappa = 1
    end
    Ω = L + (1 - sign_kappa)*π/2
        

    # argument of the pericenter
    if sign(inc) >= 0
        sign_i = 1
    else
        sign_i = -1
    end

    cos_om = (sign_i/ecc) * ( sma/rga * (1-ecc^2) - 1 )
    if abs(cos_om - 1) < 1e-8
        cos_om = 1
    end
    if abs(cos_om  + 1) < 1e-8
        cos_om = - 1
    end

    ω = acos( cos_om ) 
    if inc >= 0
        ω = -f
    else
        ω = π - f
    end

    # k_io determination, from f
    if f>=0 && f<=π
        k_io = 1
    else
        k_io = -1    
    end 
     
    coe = [sma,ecc, inc, Ω, ω, f]
    return coe
end

# TODO:
# function coeset_to_enc(coe::Vector{Float})
#     # TO IMPLEMENT
# end

"""
    pumpcrank_to_vinfbody(v∞::Float64, α::Float64, κ::Float64)

computes the flyby-body centered vinf vector given its norm and the pump and crank angles. Flyby body orbit assumed to be planar and circular
"""
function pumpcrank_to_vinfbody(v∞::Float64, α::Float64, κ::Float64)
    γ_ga = 0  # circular moon orbits assumed
    
    v∞_vec = v∞*[-sin(α)*cos(κ)*cos(γ_ga) + cos(α)*sin(γ_ga), cos(α)*cos(γ_ga) - sin(α)*cos(κ)*sin(γ_ga), sin(α)*sin(κ)]
    
    return v∞_vec
end

"""
    Kep_To_Car(sma::Float64, ecc::Float64, inc::Float64, Ω::Float64, ω::Float64, f::Float64) 
Converts the orbital elements in cartesian r,v vectors in the primary centred inertial rf. 

"""
function Kep_To_Car(sma::Float64, ecc::Float64, inc::Float64, Ω::Float64, ω::Float64, f::Float64) 

    p = sma*(1-ecc^2) # semi latus rectum
    dist = p/(1+ecc*cos(f)) # distance from the primary
    # position and velocity in pqw frame
    r_pqw = [cos(f) , sin(f) , 0.0]*dist
    v_pqw = sqrt(1/p)*[ -sin(f) , ecc + cos(f) , 0.0]

    # rotation matrices
    R1 = [cos(Ω)  sin(Ω)  0;
        -sin(Ω) cos(Ω)  0;
        0       0      1]
        
    R2 = [1  0           0;
        0  cos(inc)  sin(inc);    
        0 -sin(inc)  cos(inc)]

    R3 = [cos(ω)  sin(ω)  0;
        -sin(ω) cos(ω)  0;
        0       0      1]

    R = R3*R2*R1

    r = R'*r_pqw
    v = R'*v_pqw

    return r,v
end


"""
    compute_pump_angle(N:Integer, M:Integer, vinf::Float64, vga::Float64)

computes the pump angle α for a resonant orbit with resonance N:M and given vinf and vga
"""

function compute_pump_angle(N::Integer, M::Integer, v∞::Float64, vga::Float64)

    cos_α =  ( vga^2*(2 - (M/N)^(2/3) ) - v∞^2 - vga^2) / ( 2*v∞*vga ) 
    if cos_α >= -1 && cos_α <= 1
        α = acos(cos_α)
    else
        α = NaN
    end
    return α
end

"""
    compute_pump_angle(v∞::Float64, sma:Float64)

computes the pump angle α given v∞ and semi-major axis (all scaled).
"""
function compute_pump_angle(v∞::Float64, sma::Float64)

    cos_α =  - 1/(2*v∞*sma) + 1/(2*v∞) - v∞/2

    if cos_α >= -1 && cos_α <= 1
        α = acos(cos_α)
    else
        α = NaN
    end
    return α
end


function v∞_vec_from_α_κ(v∞::Float64, α::Float64, κ::Float64)
    return v∞.*[-sin(α)*cos(κ), 
                cos(α)]
end

function v∞_vec_from_α_κ(arc::KepArc)
    return v∞_vec_from_α_κ(arc.v∞, arc.α, arc.κ)
end