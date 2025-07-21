include("./position_tuples.jl")
using LinearAlgebra


# Equações de momento
function momentum_conservation_equation_solver(
    αØ_g::Vector{Float64},
    ρØ_g::Vector{Float64},
    uØ_g::Vector{Float64},
    αØ_l::Vector{Float64},
    ρØ_l::Vector{Float64},
    uØ_l::Vector{Float64},
    PØ_::Vector{Float64}
    )

    uØ_g_in = uØ_g[1]
    uØ_l_in = uØ_l[1]

    αg = face_maingrid(αØ_g, N)
    ρg = face_maingrid(ρØ_g, N)
    ug = face_subgrid(uØ_g, N)
    αl = face_maingrid(αØ_l, N)
    ρl = face_maingrid(ρØ_l, N)
    ul = face_subgrid(uØ_l, N)
    P = face_maingrid(PØ_, N)

    display(ul.P)

    CATHARE = @. γ * (αg.e*αl.e*ρg.e*ρl.e)*(ug.e - ul.e)^2 / (αg.e*ρl.e + αl.e*ρg.e)

    uØ_G = momentum_linear_system(αg, ρg, ug, P, CATHARE, uØ_g_in)
    uØ_L = momentum_linear_system(αl, ρl, ul, P, CATHARE, uØ_l_in)
    
    return uØ_G, uØ_L
end

function momentum_linear_system(
    αk::facemaingrid,
    ρk::facemaingrid,
    uk::facesubgrid,
    P::facemaingrid,
    CATHARE::Vector{Float64},
    uk_in::Float64
    )
    
    F_E = @. αk.E*ρk.E*uk.E
    F_W = @. αk.E*ρk.E*uk.E
    ΔF = @. F_E - F_W

    a0_e = @. (αk.e*ρk.e)/Δt
    a_E = @. max(-F_E, 0)/Δx
    a_W = @. max(F_W, 0)/Δx
    a_e = a0_e + a_E + a_W

    # Matriz A
    ## Diagonal principal
    uk_D = a_e
    uk_D_in = 1.0
    uk_D_out = 1.0
    uk_D = vcat([uk_D_in], uk_D, [uk_D_out])
    ## Diagonal superior
    uk_DU = a_E
    uk_DU_in = 0.0
    uk_DU = vcat([uk_DU_in], uk_DU)
    ## Diagonal inferior
    uk_DL = -a_W
    uk_DL_out = -1.0
    uk_DL = vcat(uk_DL, [uk_DL_out])
    ## Construção da matriz A
    uk_A = Tridiagonal(uk_DL, uk_D, uk_DU)
    
    # Vetor b
    uk_b = @. (
        + (αk.e*ρk.e*uk.e)/Δt
        + αk.e*ρk.e*g*sin(θ)
        + αk.e*(P.P - P.E)/Δx
        + CATHARE*(αk.E - αk.P)/Δx
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
    αØ_g::Vector{Float64},
    ρØ_g::Vector{Float64},
    uØ_g::Vector{Float64},
    αØ_l::Vector{Float64},
    ρØ_l::Vector{Float64},
    uØ_l::Vector{Float64},
    PØ_::Vector{Float64},
    α_g::Vector{Float64},  # valores em n
    ρ_g::Vector{Float64},  # valores em n
    α_l::Vector{Float64},  # valores em n
    ρ_l::Vector{Float64}   # valores em n
    )
    
    ρg_ref = ρØ_g[1]
    ρl_ref = ρØ_l[1]

    A_g, A_l, a_g, a_l, Ag, Al, ag, al = momentum_coefficients_for_pressure_equation(αØ_g, ρØ_g, uØ_g, αØ_l, ρØ_l, uØ_l)
    
    ωg = center_omega(uØ_g, N)
    αg = center_maingrid(αØ_g, ωg, N)
    ρg = center_maingrid(ρØ_g, ωg, N)
    ug = center_subgrid(uØ_g, N)
    ωl = center_omega(uØ_l, N)
    αl = center_maingrid(αØ_l, ωl, N)
    ρl = center_maingrid(ρØ_l, ωl, N)
    ul = center_subgrid(uØ_l, N)

    αg_n = center_maingrid(α_g, ωg, N)
    ρg_n = center_maingrid(ρ_g, ωg, N)
    αl_n = center_maingrid(α_l, ωl, N)
    ρl_n = center_maingrid(ρ_l, ωl, N)

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
            + ((αg_n.P*ρg_n.P) - (αg.P*ρg.P))/ρg_ref
            + ((αl_n.P*ρl_n.P) - (αl.P*ρl.P))/ρl_ref
        )
    )
    δP_b_in = 0.0
    δP_b_out = 0.0
    δP_b = vcat([δP_b_in], δP_b, [δP_b_out])

    # Solução do sistema linear
    δP_x = δP_A \ δP_b

    # Correção dos valores
    ## Correção das massas específicas
    ρØ_G = @. ρØ_g + (1/cG^2)*δP_x
    ρØ_G = @. ρØ_l + (1/cL^2)*δP_x
    ## Correção das velocidades
    uØ_G[2:N] = @. uØ_g[2:N] + (A_g/a_g)*(δP_x[1:N-1] - δP_x[2:N])
    uØ_G[N+1] = uØ_g[N]
    uØ_L[2:N] = @. uØ_l[2:N] + (A_l/a_l)*(δP_x[1:N-1] - δP_x[2:N])
    uØ_L[N+1] = uØ_l[N]
    ## Correção da pressão
    PØ_ = @. PØ_ + δP_x
    
    return ρØ_G, uØ_G, ρØ_L, uØ_L, PØ_
end

function momentum_coefficients_for_pressure_equation(
    αØ_g::Vector{Float64},
    ρØ_g::Vector{Float64},
    uØ_g::Vector{Float64},
    αØ_l::Vector{Float64},
    ρØ_l::Vector{Float64},
    uØ_l::Vector{Float64}
    )

    ωg = face_omega(uØ_g, N)
    αg = face_maingrid(αØ_g, ωg, N)
    ρg = face_maingrid(ρØ_g, ωg, N)
    ug = face_subgrid(uØ_g, N)
    ωl = face_omega(uØ_l, N)
    αl = face_maingrid(αØ_l, ωl, N)
    ρl = face_maingrid(ρØ_l, ωl, N)
    ul = face_subgrid(uØ_l, N)

    A_g = @. αg.e/Δx        # Face Subgrid, 2:N-1
    A_l = @. αl.e/Δx        # Face Subgrid, 2:N-1
    a_g = @. (αg.e*ρg.e)/Δt + (ωg.E/2 + ωg.P/2)*(αg.e*ρg.e*ug.e)/Δx # Face Subgrid, 2:N-1
    a_l = @. (αl.e*ρl.e)/Δt + (ωg.E/2 + ωl.P/2)*(αl.e*ρl.e*ul.e)/Δx # Face Subgrid, 2:N-1
    
    Ag = @views (
        w = A_g[2:N-1],     # Center Subgrid, 3:N
        e = A_g[1:N-2],     # Center Subgrid, 2:N-1
    )
    Al = @views (
        w = A_l[2:N-1],     # Center Subgrid, 3:N
        e = A_l[1:N-2],     # Center Subgrid, 2:N-1
    )
    ag = @views (
        w = a_g[2:N-1],     # Center Subgrid, 3:N
        e = a_g[1:N-2],     # Center Subgrid, 2:N-1
    )
    al = @views (
        w = a_l[2:N-1],     # Center Subgrid, 3:N
        e = a_l[1:N-2],     # Center Subgrid, 2:N-1
    )

    return A_g, A_l, a_g, a_l, Ag, Al, ag, al
end

# Equações de fração volumétrica
function void_fraction_equation_solver(
    αØ_g::Vector{Float64},
    ρØ_g::Vector{Float64},
    uØ_g::Vector{Float64},
    αØ_l::Vector{Float64},
    ρØ_l::Vector{Float64},
    uØ_l::Vector{Float64},
    α_g::Vector{Float64},  # valores em n
    ρ_g::Vector{Float64},  # valores em n
    α_l::Vector{Float64},  # valores em n
    ρ_l::Vector{Float64}   # valores em n
    )

    ωg = center_omega(uØ_g, N)
    ρg = center_maingrid(ρØ_g, ωg, N)
    ug = center_subgrid(uØ_g, N)
    ωl = center_omega(uØ_l, N)
    ρl = center_maingrid(ρØ_l, ωl, N)
    ul = center_subgrid(uØ_l, N)

    αg_n = center_maingrid(α_g, ωg, N)   # valores em n
    ρg_n = center_maingrid(ρ_g, ωg, N)   # valores em n
    αl_n = center_maingrid(α_l, ωl, N)   # valores em n
    ρl_n = center_maingrid(ρ_l, ωl, N)   # valores em n

    αØ_G = void_fraction_linear_system(ωg, ρg, ug, αg_n, ρg_n, αØ_g[1])
    αØ_L = void_fraction_linear_system(ωl, ρl, ul, αl_n, ρl_n, αØ_l[1])

    return αØ_G, αØ_L
end

function void_fraction_linear_system(
    ωk::centeromega,
    ρk::centermaingrid,
    uk::centersubgrid,
    αk_n::centermaingrid,
    ρk_n::centermaingrid,
    αk_in::Float64
    )

    # Matriz A
    ## Diagonal principal
    αk_D = @. (
        + ρk.P/Δt
        - (1/2 + ωk.e/2)*(ρk.e*uk.e)/Δx
        + (1/2 - ωk.w/2)*(ρk.w*uk.w)/Δx
    )
    αk_D_in = 1.0
    αk_D_out = 1.0
    αk_D = vcat([αk_D_in], αk_D, [αk_D_out])
    ## Diagonal superior
    αk_DU = @. - (1/2 + ωk.e/2)*(ρk.e*uk.e)/Δx
    αk_DU_in = 0.0
    αk_DU = append!([αk_DU_in], αk_DU)
    ## Diagonal inferior
    αk_DL = @. (1/2 + ωk.w/2)*(ρk.w*uk.w)/Δx
    αk_DL_out = -1.0
    αk_DL = append!(αk_DL, [αk_DL_out])
    ## Construção da matriz A
    αk_A = Tridiagonal(αk_DL, αk_D, αk_DU)

    # Vetor b
    αk_b = @. (αk_n.P*ρk_n.P)/Δt
    αk_b_in = αk_in
    αk_b_out = 0.0
    αk_b = vcat([αk_b_in], αk_b, [αk_b_out])

    # Solução do sistema linear
    αk_x = αk_A \ αk_b

    return αk_x
end
