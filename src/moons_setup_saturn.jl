"""
Saturn system constants
"""
μₚₛ = 3.7931e7 # Saturn gravitational parameter

# Titan
rga_Titan = 1221.87e3 # km - Titan orbital radius
μ_Titan = 8.9779e3 # Titan gravitational constant
radius_Titan = 2574.7 # km - Titan radius
vga_Titan = sqrt(μₚₛ/rga_Titan)
Tga_Titan = 2*π*sqrt(rga_Titan^3/μₚₛ)
min_altfb_Titan = 1600.0 # minimum altitude of flyby at Titan

Titan = init_moon("Titan", 1, μ_Titan, radius_Titan, rga_Titan, vga_Titan, Tga_Titan, min_altfb_Titan, μₚₛ, "#C18803")

# Rhea
rga_Rhea = 527.108e3 # km - Rhea orbital radius
μ_Rhea = 1.5394e2  # Rhea gravitational constant
radius_Rhea = 763.8 # km - Rhea radius
vga_Rhea = sqrt(μₚₛ/rga_Rhea)
Tga_Rhea = 2*π*sqrt(rga_Rhea^3/μₚₛ)
min_altfb_Rhea = 50.0 # minimum altitude of flyby at Rhea

Rhea = init_moon("Rhea", 2, μ_Rhea, radius_Rhea, rga_Rhea, vga_Rhea, Tga_Rhea, min_altfb_Rhea, μₚₛ, "#5045d8")

# Dione
rga_Dione = 377.396e3 # km - Dione orbital radius
μ_Dione = 73.110 # Dione gravitational constant
radius_Dione = 561.4 # km - Dione radius
vga_Dione = sqrt(μₚₛ/rga_Dione)
Tga_Dione = 2*π*sqrt(rga_Dione^3/μₚₛ)
min_altfb_Dione = 50.0 # minimum altitude of flyby at Dione

Dione = init_moon("Dione", 3, μ_Dione, radius_Dione, rga_Dione, vga_Dione, Tga_Dione, min_altfb_Dione, μₚₛ, "#873E23")

# Tethys
rga_Tethys = 294.619e3 # km - Tethys orbital radius
μ_Tethys = 41.209  # Tethys gravitational constant
radius_Tethys = 531.1 # km - Tethys radius
vga_Tethys = sqrt(μₚₛ/rga_Tethys)
Tga_Tethys = 2*π*sqrt(rga_Tethys^3/μₚₛ)
min_altfb_Tethys = 50.0 # minimum altitude of flyby at Tethys

Tethys = init_moon("Tethys", 4, μ_Tethys, radius_Tethys, rga_Tethys, vga_Tethys, Tga_Tethys, min_altfb_Tethys, μₚₛ, "#005B00")

# Enceladus
rga_Enceladus = 237.948e3 # km - Enceladus orbital radius
μ_Enceladus = 7.2094  # Enceladus gravitational constant
radius_Enceladus = 252.1 # km - Enceladus radius
vga_Enceladus = sqrt(μₚₛ/rga_Enceladus)
Tga_Enceladus = 2*π*sqrt(rga_Enceladus^3/μₚₛ)
min_altfb_Enceladus = 25.0 # minimum altitude of flyby at Enceladus

Enceladus = init_moon("Enceladus", 5, μ_Enceladus, radius_Enceladus, rga_Enceladus, vga_Enceladus, Tga_Enceladus, min_altfb_Enceladus, μₚₛ, "#15CED1")
