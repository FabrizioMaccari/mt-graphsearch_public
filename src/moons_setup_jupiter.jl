"""
Jupiter system constants
"""
μₚⱼ = 1.2669e8 # Jupiter gravitational parameter

# Callisto
rga_Callisto = 1.8827e6 # km - Callisto orbital radius
μ_Callisto = 7.17928936139727e3 # Callisto gravitational constant
radius_Callisto = 2410.3 # km - Callisto radius
vga_Callisto = sqrt(μₚⱼ/rga_Callisto)
Tga_Callisto = 2*π*sqrt(rga_Callisto^3/μₚⱼ)
min_altfb_Callisto = 100.0 # minimum altitude of flyby at Callisto

Callisto = init_moon("Callisto", 51, μ_Callisto, radius_Callisto, rga_Callisto, vga_Callisto, Tga_Callisto, min_altfb_Callisto, μₚⱼ, "#b29661")

# Ganymede
rga_Ganymede = 1.0704e6 # km - Ganymede orbital radius
μ_Ganymede = 9.887834453334144e3 # Ganymede gravitational constant
radius_Ganymede = 2634.1 # km - Ganymede radius
vga_Ganymede = sqrt(μₚⱼ/rga_Ganymede)
Tga_Ganymede = 2*π*sqrt(rga_Ganymede^3/μₚⱼ)
min_altfb_Ganymede = 100.0 # minimum altitude of flyby at Ganymede

Ganymede = init_moon("Ganymede", 52, μ_Ganymede, radius_Ganymede, rga_Ganymede, vga_Ganymede, Tga_Ganymede, min_altfb_Ganymede, μₚⱼ, "#915200")

# Europa
rga_Europa = 6.709e5 # km - Europa orbital radius
μ_Europa = 3.202738774922892e3 # Europa gravitational constant
radius_Europa = 1560.8  # km - Europa radius
vga_Europa = sqrt(μₚⱼ/rga_Europa)
Tga_Europa = 2*π*sqrt(rga_Europa^3/μₚⱼ)
min_altfb_Europa = 100.0 # minimum altitude of flyby at Europa

Europa = init_moon("Europa", 53, μ_Europa, radius_Europa, rga_Europa, vga_Europa, Tga_Europa, min_altfb_Europa, μₚⱼ, "#0fbdd2")

# Io
rga_Io = 4.217e5 # km - Io orbital radius
μ_Io = 5.959916033410404e3 # Io gravitational constant
radius_Io = 1821.6 # km - Io radius
vga_Io = sqrt(μₚⱼ/rga_Io)
Tga_Io = 2*π*sqrt(rga_Io^3/μₚⱼ)
min_altfb_Io = 100.0 # minimum altitude of flyby at Io

Io = init_moon("Io", 54, μ_Io, radius_Io, rga_Io, vga_Io, Tga_Io, min_altfb_Io, μₚⱼ, "#d6440c")

