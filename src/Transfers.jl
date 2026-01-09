"""
Module for the computation of vilts and nonresonant transfers.
"""

"""
    compute_nonres_transfer(moon::Moon,L_dep::Float64, v∞::Float64, M::Integer, Ne::Integer, k_io1::Integer, toll::Float64 = 1e-7, id::Integer = 0)

Given vinf, M,N, and inbound or outbound information computes the corresponding in plane resonant transfer. 
Resonant transfer solved as a zero finding problem. Non dimensional units used
Planar circular orbits of the moon assumed.
Computation in the non-dimensional framework (LU = rga, TU = Tga, VU = vga)
Output is a keplerian transfer object
"""

function compute_nonres_transfer(moon::Moon,L_dep::Float64, v∞::Float64, M::Integer, Ne::Integer, k_io1::Integer, toll::Float64 = 1e-7, id::Integer = 0)
    
    # Ma definition
    k_io1 == 1 ? Ma = M +1 : Ma = M
    
    # sma guess from resonance
    if Ma != 0 && Ne != 0
        sma_res = (Ne/M)^(2/3)
    else
        sma_res = 1 
    end

    # zero finding problem definition and solution
    prob = NonlinearProblem(time_error_nonres, sma_res, SA[Ne, Ma, v∞, k_io1])
    sol = solve(prob, alg = SimpleNewtonRaphson(), reltol=toll, abstol=toll, maxiters=200)
    
    if SciMLBase.successful_retcode(sol)
        sma = sol.u
        # Orbit characteristics recovery
        T = sma^(3/2)  # period
        ecc = sqrt( 1 - 1/sma*((3 - v∞^2 - 1/sma)/2)^2 )  # eccentricity

        f_dep = k_io1*acos( 1/ecc*( sma*(1 - ecc^2) -1 ) )  # true anomaly at departure

        # angular offset
        if k_io1 == -1
            y = abs(f_dep)/π
        else
            y = 1 - abs(f_dep)/π
        end

        # Longitude of the 2nd encounter
        L_arr = L_dep + 2*π*y

        # Time of flight
        ToF_scaled = y + Ne
        ToF = ToF_scaled*moon.Tga_dim/(86400)  # days
        # pump angle
        α = acos( ( 1 - (1/T)^(2/3) - v∞^2 ) / ( 2*v∞ ) )

        # krank angles assignment
        if k_io1 == -1  # inbound departure
            κ_1 = 0
            κ_2 = π
        else  # outbound departure
            κ_1 = π
            κ_2 = 0
        end

        # Transfer object build
        # arcs
        arc1 = keparc_vinf(moon, [v∞,α,κ_1,L_dep])
        arc2 = keparc_vinf(moon, [v∞,α,κ_2,L_arr])

        # transfer object build
        NonRes = nonres(id, Ne, M, k_io1, arc1, arc2, ToF)

    else
        NonRes = keptransfer_fail(moon)
    end

    return NonRes
end


function time_error_nonres(sma::Number, params::SVector{4,Float64})
    # unpacking
    Ne = params[1]
    Ma = params[2]
    v∞ = params[3]
    k_io1 = params[4]

    T = sma^(3/2)  # period
    ecc = sqrt( 1 - 1/sma*((3 - v∞^2 - 1/sma)/2)^2 )  # eccentricity

    cos_f = 1/ecc*( sma*(1 - ecc^2) -1 )
    if cos_f >= -1 && cos_f <=1
        f = k_io1*acos(cos_f)  # true anomaly at departure
    else
        f = NaN
    end

    cos_E = 1/ecc*( 1 - 1/sma )
    if cos_E >= -1 && cos_E <=1
        E = k_io1*acos(cos_E)  # eccentric anomaly at departure
    else
        E = NaN
    end
    τ = T/(2*π)*( E - ecc*sin(E) )  # Kepler's equation

    # offsets
    if k_io1 == - 1  # inbound departure
        y_io_ang = abs(f)/π
        y_io_time = Ma*T + 2*abs(τ) - Ne
    else  # outbound departure
        y_io_ang = 1 - abs(f)/π 
        y_io_time = Ma*T - 2*abs(τ) - Ne  
    end         
    err = y_io_time - y_io_ang

    return err
end


"""
    compute_vilt(moon::Moon, v∞1::Float64, v∞2::Float64, N::Int, M::Int, K::Int, k_io1::Float64, k_io2::Float64, k_ei::Float64, L1::Float64, toll::Float64 = 1e-7, id::Integer = 0)

    VILT computation, given the departure conditions, vinf2 and the VILT characteristics. 
Data to be provided in non-dimensional units
ASSUMPTIONS:
- moon in i=0 circular orbit
- spacecraft orbit in plane (out of plane only backflip vilts possible, more limited)

INPUTS:
moon - moon object (to be updated to a list of one or two moons, to include moon to moon vilts)
vinf1 = float/v_inf_level - vinf before the vilt
vinf2 = float/v_inf_level - vinf after the vilt 
N = int - number of full moon revs
M = int - number of full s/c revs
K = int - revolution of the burn
k_io1 - float:
    1  - for outbound departure
    -1  - for inbound departure
k_io2 - float:
    1  - for outbound arrival
    -1  - for inbound arrival
k_ei = float:
    1  - for exterior leveraging
    -1  - for interior leveraging
L1 = float - longitude of the first encounter on the planet centred inertial frame
toll = float - tolerance on the error 

OUTPUTS:
vilt - keplerian transfer object
"""
function compute_vilt(moon::Moon, v∞1::Float64, v∞2::Float64, N::Int, M::Int, K::Int, k_io1::Int, k_io2::Int, k_ei::Int, L1::Float64, toll::Float64 = 1e-7, id::Integer = 0)
    
    # adjust number of revolutions
    if k_io1 == 1 && k_io2 == -1
        Ma = M + 1  # sc revs
        Na = N + 1  # moon revs
    else 
        Ma = M  # sc revs
        Na = N  # moon revs
    end

    params = ( v∞1, v∞2, Na, Ma,K , k_ei, k_io1, k_io2)

    # sma guess
    if M != 0 && N != 0
        sma_res = (N/M)^(2/3)  # check wether to change to Na/Ma
    else
        sma_res = 1
    end
    
    # zero finding problem definition and solution
    prob = NonlinearProblem(time_error_vilt, sma_res, SVector{8}(params))
    sol = solve(prob, alg = SimpleNewtonRaphson(), reltol=toll, abstol=toll, maxiters=200)
    
    if SciMLBase.successful_retcode(sol)
        sma1 = sol.u
        # Orbit characteristics recovery
        r_la = sma1 + k_ei*sqrt( sma1^2 - 1/4*sma1*( 3 - v∞1^2 )^2 + 1/2*( 3 - v∞1^2 ) - 1/(4*sma1) )

        # a2 computation
        # polynomial coefficients
        A = (3 - v∞2^2)^2 - 8*r_la
        B = 4*r_la^2 - 2*(3 - v∞2^2)
        C = 1

        # find roots
        d = B^2 - 4*A*C
        a21 = (-B + sqrt(d)) / 2 / A
        a22 = (-B - sqrt(d)) / 2 / A

        if a21 >0 && a22<=0
            sma2 = a21
        elseif a21<=0 && a22>0
            sma2 = a22
        else
            # testing of both roots
            e21 = k_ei*( r_la/a21 - 1 )

            # check if the eccentric anomaly is computable with the first root
            # E21 = k_io2*acos( (a21 - 1)/(e21*a21) )
            cosE21 = (a21 - 1)/(e21*a21)
            if cosE21 >= -1 && cosE21 <= 1
                nanE21 = false
            else
                nanE21 = true
            end

            # check if the eccentric anomaly is computable with the second root
            # E22 = k_io2*acos( (a22 - 1)/(e22*a22) )
            e22 = k_ei*( r_la/a22 - 1 )
            cosE22 = (a22 - 1)/(e22*a22)
            if cosE22 >= -1 && cosE22 <= 1
                nanE22 = false
            else
                nanE22 = true
            end

            # choice of the root
            if e21<=0 || e21>=1 || nanE21 && e22>0 && e22<1 && !nanE22  # only a22 phisically acceptable
                sma2 = a22
            elseif e22<=0 || e22>=1 || nanE22 && e21>0 && e21<1 && !nanE21  # only a21 phisically acceptable
                sma2 = a21
            elseif e21 >0 && e21<1 && !nanE21 && e22>0 && e22<1 && !nanE22  # both phisically acceptable
                if abs(a21 - a1) < abs(a22 - a1)
                    sma2 = a21
                else
                    sma2 = a22
                end
            end
        end

        # eccentricities computation
        ecc1 = k_ei*( r_la/sma1 - 1 )
        ecc2 = k_ei*( r_la/sma2 - 1 )

        # DV cost computation
        ΔV = abs( sqrt( 2/r_la - 1/sma1 ) - sqrt( 2/r_la - 1/sma2 ) )*moon.vga_dim*1000  # m/s

        # Trajectory objects creation
        # sc velocities
        vsc1 = sqrt( 2 - 1/sma1 )
        vsc2 = sqrt( 2 - 1/sma2 )
        # pump angles, from velocities and vinfs
        α1 = acos( ( vsc1^2 - v∞1^2 - 1 )/( 2*v∞1 ) )
        α2 = acos( ( vsc2^2 - v∞2^2 - 1 )/( 2*v∞2 ) )


        # crank angles, from departure and arrival choice of the vilt, i = 0
        if k_io1 == 1
            κ1 = π
        else
            κ1 = 0
        end

        if k_io2 == 1
            κ2 = π
        else
            κ2 = 0
        end

        # true anomalies recovery, from leveraging apse and eccentricities
        f1 = k_io1*acos( 1/ecc1*( sma1*( 1 - ecc1^2 ) - 1 ) )
        f2 = k_io2*acos( 1/ecc2*( sma2*( 1 - ecc2^2 ) - 1 ) )

        # time of flight recovery, as the one of the moon
        # ToF computation
        ToF_scaled = Na + 1/(2π)*( f2 - f1 )
        ToF = ToF_scaled*moon.Tga_dim/(86400)  # days


        # longitude of the second encounter
        L2 = mod(L1 + ToF_scaled*2π,2π)

        # Trajectory objects creation
        arc1 = keparc_vinf(moon, [v∞1,α1,κ1,L1])
        arc2 = keparc_vinf(moon, [v∞2,α2,κ2,L2])

        # model object build
        VILT = vilt(id, N, M, K, k_ei, k_io1, k_io2, arc1, arc2, ΔV, ToF)
        return VILT
    else
        return keptransfer_fail(moon)
    end
end

function time_error_vilt(sma1::Number, params::SVector{8, Float64})
    # unpacking
    v∞1 = params[1]
    v∞2 = params[2]
    Na = params[3]
    Ma = params[4]
    K = params[5]
    k_ei = params[6]
    k_io1 = params[7]
    k_io2 = params[8]

    # leveraging apse computation
    FlagNan = 0
    # r_la = sma1 + real(k_ei*sqrt( Complex(sma1^2 - 1/4*sma1*( 3 - v∞1^2 )^2 + 1/2*( 3 - v∞1^2 ) - 1/(4*sma1)) ))
    sqterm = sma1^2 - 1/4*sma1*( 3 - v∞1^2 )^2 + 1/2*( 3 - v∞1^2 ) - 1/(4*sma1)
    
    if sqterm >= 0 && sma1 > 0 
        r_la = sma1 + k_ei*sqrt( sqterm )
    
        ecc1 = k_ei*( r_la/sma1 - 1 )
        cosE1 = (sma1 - 1)/(ecc1*sma1)

        if cosE1 < -1 || cosE1 > 1
            FlagNan = 1
        end
        # a2 computation
        # polynomial coefficients
        A = (3 - v∞2^2)^2 - 8*r_la
        B = 4*r_la^2 - 2*(3 - v∞2^2)
        C = 1

        # find roots
        d = B^2 - 4*A*C
        a21 = (-B + sqrt(d)) / 2 / A
        a22 = (-B - sqrt(d)) / 2 / A

        if a21 >0 && a22<=0
            sma2 = a21
            e2 = k_ei*( r_la/sma2 - 1 )
            cosE2 = (sma2 - 1)/(e2*sma2)
            if cosE2 < -1 || cosE2 > 1
                FlagNan = 1
            end
        elseif a21<=0 && a22>0
            sma2 = a22
            e2 = k_ei*( r_la/sma2 - 1 )
            cosE2 = (sma2 - 1)/(e2*sma2)
            if cosE2 < -1 || cosE2 > 1
                FlagNan = 1
            end
        else
            # testing of both roots
            e21 = k_ei*( r_la/a21 - 1 )

            # check if the eccentric anomaly is computable with the first root
            # E21 = k_io2*acos( (a21 - 1)/(e21*a21) )
            cosE21 = (a21 - 1)/(e21*a21)
            if cosE21 >= -1 && cosE21 <= 1
                nanE21 = false
            else
                nanE21 = true
            end

            # check if the eccentric anomaly is computable with the second root
            # E22 = k_io2*acos( (a22 - 1)/(e22*a22) )
            e22 = k_ei*( r_la/a22 - 1 )
            cosE22 = (a22 - 1)/(e22*a22)
            if cosE22 >= -1 && cosE22 <= 1
                nanE22 = false
            else
                nanE22 = true
            end
            
            # choice of the root 
            if (e21<=0 || e21>=1 || nanE21) && (e22>0 && e22<1 && !nanE22)  # only a22 phisically acceptable
                sma2 = a22
            elseif (e22<=0 || e22>=1 || nanE22) && (e21>0 && e21<1 && !nanE21)  # only a21 phisically acceptable
                sma2 = a21
            elseif (e21 >0 && e21<1 && !nanE21) && (e22>0 && e22<1 && !nanE22)  # both phisically acceptable
                if abs(a21 - a1) < abs(a22 - a1)
                    sma2 = a21
                else
                    sma2 = a22
                end
            else
                err = NaN
                FlagNan = 1
            end
        end

        if FlagNan == 0

            ecc2 = k_ei*( r_la/sma2 - 1 )
            cosE2 =  (sma2 - 1)/(ecc2*sma2)
            if (cosE1 < -1 || cosE1 > 1) || (cosE2 < -1 || cosE2 > 1)
                println("eccoci")
                println("sma1 " * string(sma1))
                println("cosE1 " * string(cosE1))
                println("cosE2 " * string(cosE2))
                println("cosE21 " * string(cosE21))
                println("cosE22 " * string(cosE22))
                println("e21 " * string(e21))
                println("e22 " * string(e22))
                println("a21 " * string(a21))
                println("a22 " * string(a22))
                println("sma2 " * string(sma2))
                println("ecc2 " * string(ecc2))
                println("v∞1 "   * string(params[1]))
                println("v∞2 "   * string(params[2]))
                println("Na  "   * string(params[3]))
                println("Ma  "   * string(params[4]))
                println("K"      * string(params[5]))
                println("k_ei "  * string(params[6]))
                println("k_io1 " * string(params[7]))
                println("k_io2 " * string(params[8]))
                # println(nanE21)
                # println(nanE22)
            end
            E1 = k_io1*acos( cosE1 )
            E2 = k_io2*acos( cosE2 )

            # periods
            T1 = sma1^(3/2)
            T2 = sma2^(3/2)

            # times from periapsis
            τ1 = T1/(2*pi)*( E1 - ecc1*sin(E1) )
            τ2 = T2/(2*pi)*( E2 - ecc2*sin(E2) )

            # s/c's time of flight
            # Number of sc revs
            # ToF computation - sc
            ToF_sc = τ2 - τ1 + T1*(K + ( 1 + k_ei )/4) + T2*(Ma - K - ( 1 + k_ei )/4)
            # print(ToF_sc)
            
            # true anomaly of the encounters
            cosf1 = 1/ecc1*( sma1*( 1 - ecc1^2 ) - 1 )
            cosf2 = 1/ecc2*( sma2*( 1 - ecc2^2 ) - 1 )
            if abs(cosf1 - 1) <= toll
                f1 = 0
            elseif abs(cosf1 + 1) <= toll
                f1 = π
            else
                f1 = k_io1*acos( cosf1 )
            end
            
            if abs(cosf2 - 1) <= toll
                f2 = 0
            elseif abs(cosf2 + 1) <= toll
                f2 = π
            else
                f2 = k_io2*acos( cosf2 )
            end

            # ToF computation - moon
            ToF_ga = Na + 1/(2*π)*( f2 - f1 )

            # error assessment
            err = ToF_ga - ToF_sc
        else
            err = NaN
        end 
    else
        err = NaN
    end
    return err
end

"""
    invert_vilt(vilt::KepTransfer, toll::Float64 = 1e-7, id::Integer = 0)
Given a vilts from v∞1 to v∞2, it returns the equivalent vilt (same Dv, same ToF) from v∞2 to v∞1, computed with toll tolerance.      
"""
function invert_vilt(vilt::KepTransfer, toll::Float64 = 1e-7, id::Integer = 0)
    # recover parameters of the Keplerian arcs
    moon = vilt.arc1.moon
    v∞₁_inv = vilt.arc2.v∞
    v∞₂_inv = vilt.arc1.v∞  
    
    N = vilt.N
    M = vilt.M
    # K = vilt.K
    K_inv = 0  # for now only N:M (0) vilts considered
    k_ei = vilt.k_ei
    k_io1 = vilt.k_io1
    k_io2 = vilt.k_io2
    
    # invert parameters
    if k_ei == 1
        if k_io1 == 1 && k_io2 == 1  # OO FW
            k_io1_inv = -1
            k_io2_inv = -1
            # kappa1_inv = 0
            # kappa2_inv = 0
            # K_inv = M - K - 1
        elseif k_io1 == 1 && k_io2 == -1  # OI FW
            k_io1_inv = 1
            k_io2_inv = -1
            # kappa1_inv = π
            # kappa2_inv = 0
            # K_inv = M - K
        elseif k_io1 == -1 && k_io2 == 1  # IO FW
            k_io1_inv = -1
            k_io2_inv = 1
            # kappa1_inv = 0
            # kappa2_inv = π
            # K_inv = M - K - 1
        else  # II FW
            k_io1_inv = 1
            k_io2_inv = 1
            # kappa1_inv = π
            # kappa2_inv = π
            # K_inv = M - K - 1
        end
    elseif k_ei == -1
        if k_io1 == 1 && k_io2 == 1  # OO FW
            k_io1_inv = -1
            k_io2_inv = -1
            # kappa1_inv = 0
            # kappa2_inv = 0
            # K_inv = M - K
        elseif k_io1 == 1 && k_io2 == -1  # OI FW
            k_io1_inv = 1
            k_io2_inv = -1
            # kappa1_inv = π
            # kappa2_inv = 0
            # K_inv = M - K + 1
        elseif k_io1 == -1 && k_io2 == 1  # IO FW
            k_io1_inv = -1
            k_io2_inv = 1
            # kappa1_inv = 0
            # kappa2_inv = π
            # K_inv = M - K
        else  # II FW
            k_io1_inv = 1
            k_io2_inv = 1
            # kappa1_inv = π
            # kappa2_inv = π
            # K_inv = M - K
        end
    end
    # if K_inv < 0
    #     K_inv = 0
    # end

    # compute the inverse vilt
    vilt_inv = compute_vilt(moon, v∞₁_inv, v∞₂_inv, N, M, K_inv, k_io1_inv, k_io2_inv, k_ei, 0.0, toll, id)
    return vilt_inv
end