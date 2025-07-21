using LinearAlgebra, Plots

include("./position_tuples.jl")


# Equações de momento
function momentum_conservation_equation_solver(
    αø_g::Vector{Float64},
    ρø_g::Vector{Float64},
    uø_g::Vector{Float64},
    αø_l::Vector{Float64},
    ρø_l::Vector{Float64},
    uø_l::Vector{Float64},
    α_g::Vector{Float64},
    ρ_g::Vector{Float64},
    u_g::Vector{Float64},
    α_l::Vector{Float64},
    ρ_l::Vector{Float64},
    u_l::Vector{Float64},
    Pø_::Vector{Float64}
    )

    ug_in = u_g[1]
    ul_in = u_l[1]  

    αg = face_maingrid(αø_g, N)
    ρg = face_maingrid(ρø_g, N)
    ug = face_subgrid(uø_g, N)

    αl = face_maingrid(αø_l, N)
    ρl = face_maingrid(ρø_l, N)
    ul = face_subgrid(uø_l, N)
    
    α0g = face_maingrid(α_g, N)
    ρ0g = face_maingrid(ρ_g, N)
    u0g = face_subgrid(u_g, N)
    
    α0l = face_maingrid(α_l, N)
    ρ0l = face_maingrid(ρ_l, N)
    u0l = face_subgrid(u_l, N)
    
    P = face_maingrid(Pø_, N)

    CATHARE = @. γ * (αg.e*αl.e*ρg.e*ρl.e)*(ug.e - ul.e)^2 / (αg.e*ρl.e + αl.e*ρg.e)

    uø_G = momentum_linear_system(αg, ρg, ug, α0g, ρ0g, u0g, P, CATHARE, ug_in)
    uø_L = momentum_linear_system(αl, ρl, ul, α0l, ρ0l, u0l, P, CATHARE, ul_in)
    
    return uø_G, uø_L
end

function momentum_linear_system(
    αk::facemaingrid,
    ρk::facemaingrid,
    uk::facesubgrid,
    α0k::facemaingrid,
    ρ0k::facemaingrid,
    u0k::facesubgrid,
    P::facemaingrid,
    CATHARE::Vector{Float64},
    uk_in::Float64
    )
    
    # Fluxos numéricos
    F_E = @. αk.E*ρk.E*uk.E
    F_P = @. αk.P*ρk.P*uk.P
    ΔF = F_E - F_P

    # Coeficientes pós-agrupamento
    a0_e = @. (α0k.e*ρ0k.e)*(Δx/Δt)
    a_ee = @. -min(F_E, 0)
    a_w = @. max(F_P, 0)
    a_e = a0_e + a_ee + a_w + ΔF

    # Matriz A
    ## Diagonal principal
    uk_D = a_e
    uk_D_in = 1.0
    uk_D_out = 1.0
    uk_D = vcat([uk_D_in], uk_D, [uk_D_out])
    ## Diagonal superior
    uk_DU = -a_ee
    uk_DU_in = 0.0
    uk_DU = vcat([uk_DU_in], uk_DU)
    ## Diagonal inferior
    uk_DL = -a_w
    uk_DL_out = -1.0
    uk_DL = vcat(uk_DL, [uk_DL_out])
    ## Construção da matriz A
    uk_A = Tridiagonal(uk_DL, uk_D, uk_DU)
    
    # Vetor b
    uk_b = @. (
        + a0_e*u0k.e
        + αk.e*ρk.e*g*sin(θ)*Δx
        + αk.e*(P.P - P.E)
        + CATHARE*(αk.E - αk.P)
    )
    uk_b_in = uk_in
    uk_b_out = 0.0
    uk_b = vcat([uk_b_in], uk_b, [uk_b_out])
    
    # Solução do sistema linear
    uk_x = uk_A \ uk_b

    return uk_x
end

# Equação de correção de pressão
function pressure_correction_equation_solver(
    αø_g::Vector{Float64},
    ρø_g::Vector{Float64},
    uø_g::Vector{Float64},
    αø_l::Vector{Float64},
    ρø_l::Vector{Float64},
    uø_l::Vector{Float64},
    α_g::Vector{Float64},
    ρ_g::Vector{Float64},
    α_l::Vector{Float64},
    ρ_l::Vector{Float64},
    Pø_::Vector{Float64}
    )
    
    ρg_ref = ρ_g[1]
    ρl_ref = ρ_l[1]

    A_ge, A_le, a_ge, a_le, Ag, Al, ag, al = momentum_coefficients_for_pressure_equation(   αø_g, ρø_g, uø_g, 
                                                                                            αø_l, ρø_l, uø_l,
                                                                                            α_g, ρ_g, α_l, ρ_l  )
    
    αg = center_maingrid(αø_g, N)
    ρg = center_maingrid(ρø_g, N)
    ug = center_subgrid(uø_g, N)
    αl = center_maingrid(αø_l, N)
    ρl = center_maingrid(ρø_l, N)
    ul = center_subgrid(uø_l, N)

    α0g = center_maingrid(α_g, N)
    ρ0g = center_maingrid(ρ_g, N)
    α0l = center_maingrid(α_l, N)
    ρ0l = center_maingrid(ρ_l, N)

    # Matriz A
    ## Diagonal principal
    δP_D = @. (
        + (αg.e*ρg.e)/ρg_ref * (Ag.e/ag.e)
        + (αg.w*ρg.w)/ρg_ref * (Ag.w/ag.w)
        + (Δx/Δt)*(αg.P)/ρg_ref * (1/cG^2)
        + (αl.e*ρl.e)/ρl_ref * (Al.e/al.e)
        + (αl.w*ρl.w)/ρl_ref * (Al.w/al.w)
        + (Δx/Δt)*(αl.P)/ρl_ref * (1/cL^2)
    )
    δP_D_in = 1.0
    δP_D_out = 1.0
    δP_D = vcat([δP_D_in], δP_D, [δP_D_out])
    ## Diagonal superior
    δP_DU = @. (
        - (αg.e*ρg.e)/ρg_ref * (Ag.e/ag.e)
        - (αl.e*ρl.e)/ρl_ref * (Al.e/al.e)
    )  
    δP_DU_in = -1.0
    δP_DU = vcat([δP_DU_in], δP_DU)
    ## Diagonal inferior
    δP_DL = @. (
        - (αg.w*ρg.w)/ρg_ref * (Ag.w/ag.w)
        - (αl.w*ρl.w)/ρl_ref * (Al.w/al.w)
    )
    δP_DL_out = 0.0
    δP_DL = vcat(δP_DL, [δP_DL_out])
    ## Construção da matriz A
    δP_A = Tridiagonal(δP_DL, δP_D, δP_DU)

    # Vetor b
    δP_b = @. (
        + ((αg.w*ρg.w*ug.w) - (αg.e*ρg.e*ug.e))/ρg_ref
        + ((αl.w*ρl.w*ul.w) - (αl.e*ρl.e*ul.e))/ρl_ref
        + (Δx/Δt)*(
            + ((α0g.P*ρ0g.P) - (αg.P*ρg.P))/ρg_ref
            + ((α0l.P*ρ0l.P) - (αl.P*ρl.P))/ρl_ref
        )
    )
    δP_b_in = 0.0
    δP_b_out = 0.0
    δP_b = vcat([δP_b_in], δP_b, [δP_b_out])


    #dfp = DataFrame(DL=δP_DL[1:N-1], D=δP_D[1:N-1], DU=δP_DU[1:N-1], B=δP_b[1:N-1])
    #CSV.write("exporteusimp.csv", dfp, delim=",")

    # Solução do sistema linear
    δP_x = δP_A \ δP_b

    # Correção dos valores
    ## Correção das massas específicas
    ρø_G = @. ρø_g + (1/cG^2)*δP_x
    ρø_G = @. ρø_l + (1/cL^2)*δP_x
    ## Correção das velocidades
    uø_G[2:N] = @. uø_g[2:N] + (A_ge/a_ge)*(δP_x[1:N-1] - δP_x[2:N])
    uø_G[N+1] = uø_g[N]
    uø_L[2:N] = @. uø_l[2:N] + (A_le/a_le)*(δP_x[1:N-1] - δP_x[2:N])
    uø_L[N+1] = uø_l[N]
    ## Correção da pressão
    Pø_ = @. Pø_ + δP_x
    
    return ρø_G, uø_G, ρø_L, uø_L, Pø_
end

function momentum_coefficients_for_pressure_equation(
    αø_g::Vector{Float64},
    ρø_g::Vector{Float64},
    uø_g::Vector{Float64},
    αø_l::Vector{Float64},
    ρø_l::Vector{Float64},
    uø_l::Vector{Float64},
    α0_g::Vector{Float64},
    ρ0_g::Vector{Float64},
    α0_l::Vector{Float64},
    ρ0_l::Vector{Float64}
    )

    αg = face_maingrid(αø_g, N)
    ρg = face_maingrid(ρø_g, N)
    ug = face_subgrid(uø_g, N)

    α0g = face_maingrid(α0_g, N)
    ρ0g = face_maingrid(ρ0_g, N)

    αl = face_maingrid(αø_l, N)
    ρl = face_maingrid(ρø_l, N)
    ul = face_subgrid(uø_l, N)

    α0l = face_maingrid(α0_l, N)
    ρ0l = face_maingrid(ρ0_l, N)

    # Fluxos numéricos
    ## Fase gasosa
    F_gE = @. αg.E*ρg.E*ug.E
    F_gP = @. αg.P*ρg.P*ug.P
    ΔFg = @. F_gE - F_gP
    ## Fase líquida
    F_lE = @. αl.E*ρl.E*ul.E
    F_lP = @. αl.P*ρl.P*ul.P
    ΔFl = @. F_lE - F_lP

    # Coeficientes pós-agrupamento
    ## Fase gasosa
    A_ge = @. αg.e              # Face Subgrid, 2:N
    a0_ge = @. (α0g.e*ρ0g.e)*(Δx/Δt)
    a_gee = @. -min(F_gE, 0)
    a_gw = @. max(F_gP, 0)
    a_ge = @. a0_ge + a_gee + a_gw + ΔFg   # Face Subgrid, 2:N
    ## Fase líquida
    A_le = @. αl.e              # Face Subgrid, 2:N
    a0_le = @. (α0l.e*ρ0l.e)*(Δx/Δt)
    a_lee = @. -min(F_lE, 0)
    a_lw = @. max(F_lP, 0)
    a_le = a0_le + a_lee + a_lw + ΔFl   # Face Subgrid, 2:N
    
    # Tuples de posição
    Ag = @views (
        w = A_ge[2:N-1],     # Center Subgrid, 3:N
        e = A_ge[1:N-2],     # Center Subgrid, 2:N-1
    )
    Al = @views (
        w = A_le[2:N-1],     # Center Subgrid, 3:N
        e = A_le[1:N-2],     # Center Subgrid, 2:N-1
    )
    ag = @views (
        w = a_ge[2:N-1],     # Center Subgrid, 3:N
        e = a_ge[1:N-2],     # Center Subgrid, 2:N-1
    )
    al = @views (
        w = a_le[2:N-1],     # Center Subgrid, 3:N
        e = a_le[1:N-2],     # Center Subgrid, 2:N-1
    )

    return A_ge, A_le, a_ge, a_le, Ag, Al, ag, al
end

# Equações de fração volumétrica
function void_fraction_equation_solver(
    ρø_g::Vector{Float64},
    uø_g::Vector{Float64},
    ρø_l::Vector{Float64},
    uø_l::Vector{Float64},
    α_g::Vector{Float64},
    ρ_g::Vector{Float64},
    α_l::Vector{Float64},
    ρ_l::Vector{Float64}
    )

    αg_in = α_g[1]
    αl_in = α_l[1] 

    ρg = center_maingrid(ρø_g, N)
    ug = center_subgrid(uø_g, N)
    ρl = center_maingrid(ρø_l, N)
    ul = center_subgrid(uø_l, N)

    α0g = center_maingrid(α_g, N)
    ρ0g = center_maingrid(ρ_g, N)
    α0l = center_maingrid(α_l, N)
    ρ0l = center_maingrid(ρ_l, N)

    αø_G = void_fraction_linear_system(ρg, ug, α0g, ρ0g, αg_in)
    αø_L = void_fraction_linear_system(ρl, ul, α0l, ρ0l, αl_in)

    return αø_G, αø_L
end

function void_fraction_linear_system(
    ρk::centermaingrid,
    uk::centersubgrid,
    α0k::centermaingrid,
    ρ0k::centermaingrid,
    αk_in::Float64
    )

    # Fluxos numéricos
    F_e = @. ρk.e*uk.e
    F_w = @. ρk.w*uk.w
    ΔF = F_e - F_w

    # Coeficientes pós-agrupamento
    a0_P = @. ρ0k.P*(Δx/Δt)
    a_E = @. -min(F_e, 0)
    a_W = @. max(F_w, 0)
    a_P = ρk.P*(Δx/Δt) - a_W - a_E - ΔF

    # Matriz A
    ## Diagonal principal
    αk_D = a_P
    αk_D_in = 1.0
    αk_D_out = 1.0
    αk_D = vcat([αk_D_in], αk_D, [αk_D_out])
    ## Diagonal superior
    αk_DU = a_E
    αk_DU_in = 0.0
    αk_DU = append!([αk_DU_in], αk_DU)
    ## Diagonal inferior
    αk_DL = a_W
    αk_DL_out = -1.0
    αk_DL = append!(αk_DL, [αk_DL_out])
    ## Construção da matriz A
    αk_A = Tridiagonal(αk_DL, αk_D, αk_DU)

    # Vetor b
    αk_b = @. a0_P*α0k.P
    αk_b_in = αk_in
    αk_b_out = 0.0
    αk_b = vcat([αk_b_in], αk_b, [αk_b_out])

    # Solução do sistema linear
    αk_x = αk_A \ αk_b

    return αk_x
end
