include("./position_tuples.jl")


# Equações de momento
function momentum_conservation_equation_solver(
    α0_G::Vector{Float64},
    α0_L::Vector{Float64},
    u0_G::Vector{Float64},
    u0_L::Vector{Float64},
    αx_G::Vector{Float64},
    αx_L::Vector{Float64},
    ux_G::Vector{Float64},
    ux_L::Vector{Float64},
    Px_::Vector{Float64},
    ρ_G::Float64,
    ρ_L::Float64,
    N::Int64
    )

    uG_in = ux_G[1]
    uL_in = ux_L[1]  

    α0G = face_maingrid(α0_G, N)
    α0L = face_maingrid(α0_L, N)
    u0G = face_subgrid(u0_G, N)
    u0L = face_subgrid(u0_L, N)

    αxG = face_maingrid(αx_G, N)
    αxL = face_maingrid(αx_L, N)
    uxG = face_subgrid(ux_G, N)
    uxL = face_subgrid(ux_L, N)
    
    Px = face_maingrid(Px_, N)

    CATHARE = @. γ * (αxG.e*αxL.e*ρ_G*ρ_L)*(uxG.e - uxL.e)^2 / (αxG.e*ρ_L + αxL.e*ρ_G)

    ux_G = momentum_linear_system(α0G, u0G, αxG, uxG, Px, CATHARE, ρ_G, uG_in)
    ux_L = momentum_linear_system(α0L, u0L, αxL, uxL, Px, CATHARE, ρ_L, uL_in)
    
    return ux_G, ux_L
end

function momentum_linear_system(
    α0k::facemaingrid,
    u0k::facesubgrid,
    αxk::facemaingrid,
    uxk::facesubgrid,
    Px::facemaingrid,
    CATHARE::Vector{Float64},
    ρ_k::Float64,
    uk_in::Float64
    )
    
    # Fluxos numéricos
    F_E = @. ρ_k*αxk.E*uxk.E
    F_P = @. ρ_k*αxk.P*uxk.P
    ΔF = F_E - F_P
    
    # Coeficientes pós-agrupamento
    A_e = α0k.e
    a0_e = @. ρ_k*α0k.e*(Δx/Δt)
    a_w = @. max(F_P, 0)
    a_ee = @. max(0, -F_E)
    a_e = a0_e + a_w + a_ee + ΔF

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
        + A_e*(Px.P - Px.E)
        + a0_e*u0k.e
        + ρ_k*αxk.e*g*sin(θ)*Δx
        + CATHARE*(αxk.E - αxk.P)
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
    α0_G::Vector{Float64},
    α0_L::Vector{Float64},
    αx_G::Vector{Float64},
    αx_L::Vector{Float64},
    ux_G::Vector{Float64},
    ux_L::Vector{Float64},
    Px_::Vector{Float64},
    ρ_G::Float64,
    ρ_L::Float64,
    N::Int64
    )
    
    A_Ge, A_Le, a_Ge, a_Le, AG, AL, aG, aL = momentum_coefficients_for_pressure_equation(α0_G, α0_L, αx_G, αx_L, ux_G, ux_L, ρ_G, ρ_L, N)

    α0G = center_maingrid(α0_G, N)
    α0L = center_maingrid(α0_L, N)
    
    αxG = center_maingrid(αx_G, N)
    αxL = center_maingrid(αx_L, N)
    uxG = center_subgrid(ux_G, N)
    uxL = center_subgrid(ux_L, N)

    # Matriz A
    ## Diagonal principal
    δP_D = @. (
        + αxG.e * (AG.e/aG.e)
        + αxG.w * (AG.w/aG.w)
        + αxL.e * (AL.e/aL.e)
        + αxL.w * (AL.w/aL.w)
    )
    δP_D_in = 1.0
    δP_D_out = 1.0
    δP_D = vcat([δP_D_in], δP_D, [δP_D_out])
    ## Diagonal superior
    δP_DU = @. (
        - αxG.e * (AG.e/aG.e)
        - αxL.e * (AL.e/aL.e)
    )  
    δP_DU_in = -1.0
    δP_DU = vcat([δP_DU_in], δP_DU)
    ## Diagonal inferior
    δP_DL = @. (
        - αxG.w * (AG.w/aG.w)
        - αxL.w * (AL.w/aL.w)
    )
    δP_DL_out = 0.0
    δP_DL = vcat(δP_DL, [δP_DL_out])
    ## Construção da matriz A
    δP_A = Tridiagonal(δP_DL, δP_D, δP_DU)

    # Vetor b
    δP_b = @. (
        + αxG.w*uxG.w - αxG.e*uxG.e
        + αxL.w*uxL.w - αxL.e*uxL.e
        + (Δx/Δt)*(
            + α0G.P - αxG.P
            + α0L.P - αxL.P
        )
    )
    δP_b_in = 0.0
    δP_b_out = 0.0
    δP_b = vcat([δP_b_in], δP_b, [δP_b_out])

    # Solução do sistema linear
    δP_x = δP_A \ δP_b

    # Correção dos valores
    ## Correção das velocidades
    u_G[2:N] = @. ux_G[2:N] + (A_Ge/a_Ge)*(δP_x[1:N-1] - δP_x[2:N])
    u_G[N+1] = ux_G[N]
    u_L[2:N] = @. ux_L[2:N] + (A_Le/a_Le)*(δP_x[1:N-1] - δP_x[2:N])
    u_L[N+1] = ux_L[N]
    ## Correção da pressão
    P_ = @. Px_ + δP_x
    
    return u_G, u_L, P_
end

function momentum_coefficients_for_pressure_equation(
    α0_G::Vector{Float64},
    α0_L::Vector{Float64},
    αx_G::Vector{Float64},
    αx_L::Vector{Float64},
    ux_G::Vector{Float64},
    ux_L::Vector{Float64},
    ρ_G::Float64,
    ρ_L::Float64,
    N::Int64
    )

    α0G = face_maingrid(α0_G, N)
    α0L = face_maingrid(α0_L, N)
    
    αxG = face_maingrid(αx_G, N)
    αxL = face_maingrid(αx_L, N)
    uxG = face_subgrid(ux_G, N)
    uxL = face_subgrid(ux_L, N)

    # Fluxos numéricos
    ## Fase gasosa
    F_GE = @. αxG.E*uxG.E
    F_GP = @. αxG.P*uxG.P
    ΔFG = @. F_GE - F_GP
    ## Fase líquida
    F_LE = @. αxL.E*uxL.E
    F_LP = @. αxL.P*uxL.P
    ΔFL = @. F_LE - F_LP

    # Coeficientes pós-agrupamento
    ## Fase gasosa
    A_Ge = @. αxG.e/ρ_G                     # Face Subgrid, 2:N
    a0_Ge = @. α0G.e*(Δx/Δt)
    a_Gw = @. max(F_GP, 0)
    a_Gee = @. max(0, -F_GE)
    a_Ge = @. a0_Ge + a_Gw + a_Gee + ΔFG    # Face Subgrid, 2:N
    ## Fase líquida
    A_Le = @. αxL.e/ρ_L                     # Face Subgrid, 2:N
    a0_Le = @. α0L.e*(Δx/Δt)
    a_Lw = @. max(F_LP, 0)
    a_Lee = @. max(0, -F_LE)
    a_Le = a0_Le + a_Lw + a_Lee + ΔFL       # Face Subgrid, 2:N
    
    # Tuples de posição
    AG = @views (
        w = A_Ge[2:N-1],     # Center Subgrid, 3:N
        e = A_Ge[1:N-2],     # Center Subgrid, 2:N-1
    )
    AL = @views (
        w = A_Le[2:N-1],     # Center Subgrid, 3:N
        e = A_Le[1:N-2],     # Center Subgrid, 2:N-1
    )
    aG = @views (
        w = a_Ge[2:N-1],     # Center Subgrid, 3:N
        e = a_Ge[1:N-2],     # Center Subgrid, 2:N-1
    )
    aL = @views (
        w = a_Le[2:N-1],     # Center Subgrid, 3:N
        e = a_Le[1:N-2],     # Center Subgrid, 2:N-1
    )

    return A_Ge, A_Le, a_Ge, a_Le, AG, AL, aG, aL
end

# Equações de fração volumétrica
function void_fraction_equation_solver(
    α0_G::Vector{Float64},
    α0_L::Vector{Float64},
    u_G::Vector{Float64},
    u_L::Vector{Float64},
    ρ_G::Float64,
    ρ_L::Float64,
    N::Int64
    )

    αG_in = α0_G[1]
    αL_in = α0_L[1] 

    α0G = center_maingrid(α0_G, N)
    α0L = center_maingrid(α0_L, N)

    uG = center_subgrid(u_G, N)
    uL = center_subgrid(u_L, N)

    α_G = void_fraction_linear_system(α0G, uG, ρ_G, αG_in)
    α_L = void_fraction_linear_system(α0L, uL, ρ_L, αL_in)

    return α_G, α_L
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
    a_W = @. max(F_w, 0)
    a_E = @. max(0, -F_e)
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
