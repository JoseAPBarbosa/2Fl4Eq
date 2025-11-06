include("./position_structs.jl")


#=============================================================================#
# Equações de momento                                                         #
#=============================================================================#

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
    ρ_L::Float64
    )

    α0G = AlphaFace(α0_G, ux_G)
    α0L = AlphaFace(α0_L, ux_L)
    u0G = UVelFace(u0_G)
    u0L = UVelFace(u0_L)

    αxG = AlphaFace(αx_G, ux_G)
    αxL = AlphaFace(αx_L, ux_L)
    uxG = UVelFace(ux_G)
    uxL = UVelFace(ux_L)
    Px = PressFace(Px_)

    CATHARE = @. γ*((ρ_G*αxG.e*αxL.e*ρ_L)/(ρ_L*αxG.e + ρ_G*αxL.e))*(uxG.e - uxL.e)^2

    ux_G = momentum_linear_system(α0G, u0G, αxG, uxG, Px, ρ_G, CATHARE)
    ux_L = momentum_linear_system(α0L, u0L, αxL, uxL, Px, ρ_L, CATHARE)
    
    return ux_G, ux_L
end

function momentum_linear_system(
    α0k::AlphaFace,
    u0k::UVelFace,
    αxk::AlphaFace,
    uxk::UVelFace,
    Px::PressFace,
    ρ_k::Float64,
    CATHARE::Vector{Float64}
    )
    
    # Fluxos numéricos
    F_E = @. ρ_k*αxk.E*uxk.E/Δx
    F_P = @. ρ_k*αxk.P*uxk.P/Δx
    ΔF = F_E - F_P
    
    # Coeficientes pós-agrupamento
    A_e = αxk.e
    a0_e = @. ρ_k*α0k.e/Δt
    a_w = @. max(F_P, 0)
    a_ee = @. max(0, -F_E)
    a_e = a0_e + a_w + a_ee + ΔF

    # Matriz A
    ## Diagonal principal
    uk_D = zeros(N+1)
    uk_D[2:N] = a_e[:]
    uk_D[1] = 1.0
    uk_D[N+1] = 1.0
    ## Diagonal superior
    uk_DU = zeros(N)
    uk_DU[2:N] = -a_ee[:]
    uk_DU[1] = 0.0
    ## Diagonal inferior
    uk_DL = zeros(N)
    uk_DL[1:N-1] = -a_w[:]
    uk_DL[N] = -1.0
    ## Construção da matriz A
    uk_A = Tridiagonal(uk_DL, uk_D, uk_DU)
    
    # Vetor b
    uk_b = zeros(N+1)
    uk_b[2:N] = @. (
        + A_e*(Px.P - Px.E)
        + a0_e*u0k.e
        + CATHARE*(αxk.E - αxk.P)
        + ρ_k*αxk.e*g*sin(θ)
    )
    uk_b[1] = u0k.in
    uk_b[N+1] = 0.0

    # Solução do sistema linear
    uk_x = uk_A \ uk_b

    return uk_x
end


#=============================================================================#
# Equação de correção de pressão                                              #
#=============================================================================#

function pressure_correction_equation_solver(
    α0_G::Vector{Float64},
    α0_L::Vector{Float64},
    αx_G::Vector{Float64},
    αx_L::Vector{Float64},
    ux_G::Vector{Float64},
    ux_L::Vector{Float64},
    Px_::Vector{Float64},
    ρ_G::Float64,
    ρ_L::Float64
    )

    α0G = AlphaCenter(α0_G, ux_G)
    α0L = AlphaCenter(α0_L, ux_L)
    αxG = AlphaCenter(αx_G, ux_G)
    αxL = AlphaCenter(αx_L, ux_L)
    uxG = UVelCenter(ux_G)
    uxL = UVelCenter(ux_L)

    # Coeficientes da equação de momento
    D_G, DG = momentum_coefficients_for_pressure_equation(α0_G, αx_G, ux_G, ρ_G)
    D_L, DL = momentum_coefficients_for_pressure_equation(α0_L, αx_L, ux_L, ρ_L)

    # Matriz A
    ## Diagonal principal
    δP_D = zeros(N)
    δP_D[2:N-1] = @. (
        + αxG.e * DG.e
        + αxG.w * DG.w
        + αxL.e * DL.e
        + αxL.w * DL.w
    )
    δP_D[1] = 1.0
    δP_D[N] = 1.0
    ## Diagonal superior
    δP_DU = zeros(N-1)
    δP_DU[2:N-1] = @. (
        - αxG.e * DG.e
        - αxL.e * DL.e
    )  
    δP_DU[1] = -1.0
    ## Diagonal inferior
    δP_DL = zeros(N-1)
    δP_DL[1:N-2] = @. (
        - αxG.w * DG.w
        - αxL.w * DL.w
    )
    δP_DL[N-1] = 0.0
    ## Construção da matriz A
    δP_A = Tridiagonal(δP_DL, δP_D, δP_DU)
    
    # Vetor b
    δP_b = zeros(N)
    δP_b[2:N-1] = @. (
        + αxG.w*uxG.w - αxG.e*uxG.e
        + αxL.w*uxL.w - αxL.e*uxL.e
        + (Δx/Δt)*(
            + α0G.P - αxG.P
            + α0L.P - αxL.P
        )
    )
    δP_b[1] = 0.0
    δP_b[N] = 0.0

    # Solução do sistema linear
    δP_x = δP_A \ δP_b

    # Correção dos valores
    ## Correção das velocidades
    u_G[2:N] = @. ux_G[2:N] + D_G*(δP_x[1:N-1] - δP_x[2:N])
    u_G[N+1] = u_G[N]
    u_L[2:N] = @. ux_L[2:N] + D_L*(δP_x[1:N-1] - δP_x[2:N])
    u_L[N+1] = u_L[N]
    ## Correção da pressão
    P_ = @. Px_ + δP_x
    
    return u_G, u_L, P_
end

function momentum_coefficients_for_pressure_equation(
    α0k::Vector{Float64},
    αxk::Vector{Float64},
    uxk::Vector{Float64},
    ρ_k::Float64
    )

    α0k = AlphaFace(α0k, uxk)
    αxk = AlphaFace(αxk, uxk)
    uxk = UVelFace(uxk)

    # Fluxos numéricos
    ## Fase gasosa
    F_E = @. αxk.E*uxk.E/Δx
    F_P = @. αxk.P*uxk.P/Δx
    ΔF = @. F_E - F_P

    # Coeficientes pós-agrupamento
    ## Fase gasosa
    A_e = @. αxk.e/ρ_k
    a0_e = @. α0k.e/Δt
    a_w = @. max(F_P, 0)
    a_ee = @. max(0, -F_E)
    a_e = @. a0_e + a_w + a_ee + ΔF
    D_k = @. A_e/a_e

    # Tuple de posição
    Dk = @views (
        w = D_k[1:N-2],
        e = D_k[2:N-1],
    )

    return D_k, Dk
end

#=============================================================================#
# Equações de fração volumétrica                                              #
#=============================================================================#

function void_fraction_equation_solver(
    α0_G::Vector{Float64},
    α0_L::Vector{Float64},
    u_G::Vector{Float64},
    u_L::Vector{Float64},
    ρ_G::Float64,
    ρ_L::Float64
    )

    α0G = AlphaCenter(α0_G, u_G)
    α0L = AlphaCenter(α0_L, u_L)
    uG = UVelCenter(u_G)
    uL = UVelCenter(u_L)

    αx_G = void_fraction_linear_system(α0G, uG, ρ_G)
    αx_L = void_fraction_linear_system(α0L, uL, ρ_L)

    α_L[2:N] = @. αx_L[2:N] / (αx_L[2:N] + αx_G[2:N])
    α_G[2:N] = @. 1.0 - α_L[2:N]

    return α_G, α_L
end

function void_fraction_linear_system(
    α0k::AlphaCenter,
    uk::UVelCenter,
    ρ_k::Float64
    )

    # Fluxos numéricos
    F_e = @. ρ_k*uk.e/Δx
    F_w = @. ρ_k*uk.w/Δx
    ΔF = F_e - F_w

    # Coeficientes pós-agrupamento
    a0_P = ρ_k/Δt
    a_W = @. max(F_w, 0)
    a_E = @. max(0, -F_e)
    a_P = @. a0_P + a_W + a_E + ΔF

    # Matriz A
    ## Diagonal principal
    αk_D = zeros(N)
    αk_D[2:N-1] = a_P[:]
    αk_D[1] = 1.0
    αk_D[N] = 1.0
    ## Diagonal superior
    αk_DU = zeros(N-1)
    αk_DU[2:N-1] = -a_E[:]
    αk_DU[1] = 0.0
    ## Diagonal inferior
    αk_DL = zeros(N-1)
    αk_DL[1:N-2] = -a_W[:]
    αk_DL[N-1] = -1.0
    ## Construção da matriz A
    αk_A = Tridiagonal(αk_DL, αk_D, αk_DU)

    # Vetor b
    αk_b = zeros(N)
    αk_b[2:N-1] = @. a0_P*α0k.P
    αk_b[1] = α0k.in
    αk_b[N] = 0.0

    # Solução do sistema linear
    αk_x = αk_A \ αk_b

    return αk_x
end
