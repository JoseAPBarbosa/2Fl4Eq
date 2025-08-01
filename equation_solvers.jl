using LinearAlgebra, Plots

include("./position_tuples.jl")


# Equações de momento
function momentum_conservation_equation_solver(
    α0_G::Vector{Float64},
    α0_L::Vector{Float64},
    u0_G::Vector{Float64},
    u0_L::Vector{Float64},
    αg_G::Vector{Float64},
    αg_L::Vector{Float64},
    ug_G::Vector{Float64},
    ug_L::Vector{Float64},
    Pg_::Vector{Float64},
    ρ_G::Float64,
    ρ_L::Float64,
    N::Int64
    )

    u_G_in = ug_G[1]
    u_L_in = ug_L[1]  

    α0G = face_maingrid(α0_G, N)
    α0L = face_maingrid(α0_L, N)
    u0G = face_subgrid(u0_G, N)
    u0L = face_subgrid(u0_L, N)

    αgG = face_maingrid(αg_G, N)
    αgL = face_maingrid(αg_L, N)
    ugG = face_subgrid(ug_G, N)
    ugL = face_subgrid(ug_L, N)
    
    Pg = face_maingrid(Pg_, N)

    CATHARE = @. γ * (αgG.e*αgL.e*ρ_G*ρ_L)*(ugG.e - ugL.e)^2 / (αgG.e*ρ_L + αgL.e*ρ_G)

    uø_G = momentum_linear_system(α0G, u0G, αgG, ugG, Pg, CATHARE, ρ_G, u_G_in)
    uø_L = momentum_linear_system(α0L, u0L, αgL, ugL, Pg, CATHARE, ρ_L, u_L_in)
    
    return uø_G, uø_L
end

function momentum_linear_system(
    α0k::facemaingrid,
    u0k::facesubgrid,
    αgk::facemaingrid,
    ugk::facesubgrid,
    Pg::facemaingrid,
    CATHARE::Vector{Float64},
    ρ_k::Float64,
    u_k_in::Float64
    )
    
    # Fluxos numéricos
    F_E = @. ρ_k*(αgk.E*ugk.E)
    F_P = @. ρ_k*(αgk.P*ugk.P)
    ΔF = F_E - F_P

    # Coeficientes pós-agrupamento
    a0_e = @. ρ_k*(α0k.e)*(Δx/Δt)
    a_ee = @. max(0, -F_E)
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
        + ρ_k*(αgk.e*g*sin(θ))*Δx
        + αgk.e*(Pg.P - Pg.E)
        + CATHARE*(αgk.E - αgk.P)
    )
    uk_b_in = u_k_in
    uk_b_out = 0.0
    uk_b = vcat([uk_b_in], uk_b, [uk_b_out])
    
    # Solução do sistema linear
    uk_x = uk_A \ uk_b

    return uk_x
end

# Equação de correção de pressão
function pressure_correction_equation_solver(
    αg_G::Vector{Float64},
    αg_L::Vector{Float64},
    α0_G::Vector{Float64},
    α0_L::Vector{Float64},
    ug_G::Vector{Float64},
    ug_L::Vector{Float64},
    Pg_::Vector{Float64},
    ρ_G::Float64,
    ρ_L::Float64,
    N::Int64
    )
    ρ_G_ref = ρ_g
    ρ_L_ref = ρ_l
    
    A_ge, A_le, a_ge, a_le, Ag, Al, ag, al = momentum_coefficients_for_pressure_equation(   αø_g, αø_l, α_g, α_l,
                                                                                            uø_g, uø_l,
                                                                                            ρ_g, ρ_l    )
    
    αg = center_maingrid(αø_g, N)
    αl = center_maingrid(αø_l, N)
    α0g = center_maingrid(α_g, N)
    α0l = center_maingrid(α_l, N)
    
    ug = center_subgrid(uø_g, N)
    ul = center_subgrid(uø_l, N)

    # Matriz A
    ## Diagonal principal
    δP_D = @. (
        + ρ_g*(αg.e)/ρg_ref * (Ag.e/ag.e)
        + ρ_g*(αg.w)/ρg_ref * (Ag.w/ag.w)
        + ρ_l*(αl.e)/ρl_ref * (Al.e/al.e)
        + ρ_l*(αl.w)/ρl_ref * (Al.w/al.w)
    )
    δP_D_in = 1.0
    δP_D_out = 1.0
    δP_D = vcat([δP_D_in], δP_D, [δP_D_out])
    ## Diagonal superior
    δP_DU = @. (
        - ρ_g*(αg.e)/ρg_ref * (Ag.e/ag.e)
        - ρ_l*(αl.e)/ρl_ref * (Al.e/al.e)
    )  
    δP_DU_in = -1.0
    δP_DU = vcat([δP_DU_in], δP_DU)
    ## Diagonal inferior
    δP_DL = @. (
        - ρ_g*(αg.w)/ρg_ref * (Ag.w/ag.w)
        - ρ_l*(αl.w)/ρl_ref * (Al.w/al.w)
    )
    δP_DL_out = 0.0
    δP_DL = vcat(δP_DL, [δP_DL_out])
    ## Construção da matriz A
    δP_A = Tridiagonal(δP_DL, δP_D, δP_DU)

    # Vetor b
    δP_b = @. (
        + ρ_g*((αg.w*ug.w) - (αg.e*ug.e))/ρg_ref
        + ρ_l*((αl.w*ul.w) - (αl.e*ul.e))/ρl_ref
        + (Δx/Δt)*(
            + ρ_g*(α0g.P - αg.P)/ρg_ref
            + ρ_l*(α0l.P - αl.P)/ρl_ref
        )
    )
    δP_b_in = 0.0
    δP_b_out = 0.0
    δP_b = vcat([δP_b_in], δP_b, [δP_b_out])

    # Solução do sistema linear
    δP_x = δP_A \ δP_b

    # Correção dos valores
    ## Correção das velocidades
    uø_G[2:N] = @. uø_g[2:N] + (A_ge/a_ge)*(δP_x[1:N-1] - δP_x[2:N])
    uø_G[N+1] = uø_g[N]
    uø_L[2:N] = @. uø_l[2:N] + (A_le/a_le)*(δP_x[1:N-1] - δP_x[2:N])
    uø_L[N+1] = uø_l[N]
    ## Correção da pressão
    Pø_ = @. Pø_ + δP_x
    
    return uø_G, uø_L, Pø_
end

function momentum_coefficients_for_pressure_equation(
    αø_g::Vector{Float64},
    αø_l::Vector{Float64},
    α0_g::Vector{Float64},
    α0_l::Vector{Float64},
    uø_g::Vector{Float64},
    uø_l::Vector{Float64},
    ρ_g::Float64,
    ρ_l::Float64
    )

    αg = face_maingrid(αø_g, N)
    αl = face_maingrid(αø_l, N)
    α0g = face_maingrid(α0_g, N)
    α0l = face_maingrid(α0_l, N)

    ug = face_subgrid(uø_g, N)
    ul = face_subgrid(uø_l, N)

    # Fluxos numéricos
    ## Fase gasosa
    F_gE = @. ρ_g*(αg.E*ug.E)
    F_gP = @. ρ_g*(αg.P*ug.P)
    ΔFg = @. F_gE - F_gP
    ## Fase líquida
    F_lE = @. ρ_l*(αl.E*ul.E)
    F_lP = @. ρ_l*(αl.P*ul.P)
    ΔFl = @. F_lE - F_lP

    # Coeficientes pós-agrupamento
    ## Fase gasosa
    A_ge = @. αg.e                          # Face Subgrid, 2:N
    a0_ge = @. ρ_g*(α0g.e)*(Δx/Δt)
    a_gee = @. max(0, -F_gE)
    a_gw = @. max(F_gP, 0)
    a_ge = @. a0_ge + a_gee + a_gw + ΔFg    # Face Subgrid, 2:N
    ## Fase líquida
    A_le = @. αl.e                          # Face Subgrid, 2:N
    a0_le = @. ρ_l*(α0l.e)*(Δx/Δt)
    a_lee = @. max(0, -F_lE)
    a_lw = @. max(F_lP, 0)
    a_le = a0_le + a_lee + a_lw + ΔFl       # Face Subgrid, 2:N
    
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
    α_g::Vector{Float64},
    α_l::Vector{Float64},
    uø_g::Vector{Float64},
    uø_l::Vector{Float64},
    ρ_g::Float64,
    ρ_l::Float64
    )

    α_g_in = α_g[1]
    α_l_in = α_l[1] 

    α0g = center_maingrid(α_g, N)
    α0l = center_maingrid(α_l, N)

    ug = center_subgrid(uø_g, N)
    ul = center_subgrid(uø_l, N)

    αø_G = void_fraction_linear_system(α0g, ug, ρ_g, α_g_in)
    αø_L = void_fraction_linear_system(α0l, ul, ρ_l, α_l_in)

    return αø_G, αø_L
end

function void_fraction_linear_system(
    α0k::centermaingrid,
    uk::centersubgrid,
    ρ_k::Float64,
    αk_in::Float64
    )

    # Fluxos numéricos
    F_e = @. ρ_k*uk.e
    F_w = @. ρ_k*uk.w
    ΔF = F_e - F_w

    # Coeficientes pós-agrupamento
    a0_P = ρ_k*(Δx/Δt)
    a_E = @. max(0, -F_e)
    a_W = @. max(F_w, 0)
    a_P = @. a0_P - a_W - a_E - ΔF

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
