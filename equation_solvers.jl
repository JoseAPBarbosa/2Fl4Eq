include("./position_structs.jl")
include("./muscl.jl")


#=============================================================================#
# Equações de momento                                                         #
#=============================================================================#

function guessed_velocity_equation_solver(
    α0_G :: Vector{Float64},
    ρ0_G :: Vector{Float64},
    u0_G :: Vector{Float64},
    α0_L :: Vector{Float64},
    ρ0_L :: Vector{Float64},
    u0_L :: Vector{Float64},
    αx_G :: Vector{Float64},
    ρx_G :: Vector{Float64},
    ux_G :: Vector{Float64},
    αx_L :: Vector{Float64},
    ρx_L :: Vector{Float64},
    ux_L :: Vector{Float64},
    Px :: Vector{Float64},
    g :: Float64,
    θ :: Float64,
    γ :: Float64,
    N :: Int64,
    Δx :: Float64,
    Δt :: Float64
    )
    
    α0G = AlphaSubGrid(α0_G, u0_G, N)
    ρ0G = RhoSubGrid(ρ0_G, u0_G, N)
    u0G = UVelSubGrid(u0_G, N)

    α0L = AlphaSubGrid(α0_L, u0_L, N)
    ρ0L = RhoSubGrid(ρ0_L, u0_L, N)
    u0L = UVelSubGrid(u0_L, N)

    αxG = AlphaSubGrid(αx_G, ux_G, N)
    ρxG = RhoSubGrid(ρx_G, ux_G, N)
    uxG = UVelSubGrid(ux_G, N)

    αxL = AlphaSubGrid(αx_L, ux_L, N)
    ρxL = RhoSubGrid(ρx_L, ux_L, N)
    uxL = UVelSubGrid(ux_L, N)

    P  = PressSubGrid(Px, N)

    CATHARE = @. γ*((ρxG.e*αxG.e*αxL.e*ρxL.e)/(ρxL.e*αxG.e + ρxG.e*αxL.e))*(uxG.e - uxL.e)^2

    ux_G = momentum_muscl(ux_G, u0_G, α0G, ρ0G, u0G, P, CATHARE, g, θ, N, Δx, Δt)
    ux_L = momentum_muscl(ux_L, u0_L, α0L, ρ0L, u0L, P, CATHARE, g, θ, N, Δx, Δt)

    return ux_G, ux_L
end


#=============================================================================#
# Equação de correção de pressão                                              #
#=============================================================================#

function pressure_correction_equation_solver(
    α0_G :: Vector{Float64},
    ρ0_G :: Vector{Float64},
    α0_L :: Vector{Float64},
    ρ0_L :: Vector{Float64},
    αx_G :: Vector{Float64},
    ρx_G :: Vector{Float64},
    ux_G :: Vector{Float64},
    αx_L :: Vector{Float64},
    ρx_L :: Vector{Float64},
    ux_L :: Vector{Float64},
    Px :: Vector{Float64},
    CG :: Float64,
    CL :: Float64,
    N  :: Int64,
    Δx :: Float64,
    Δt :: Float64
    )

    ρG_ref = ρ0_G[1]
    ρL_ref = ρ0_L[1]

    α0G = AlphaMainGrid(α0_G, u0_G, N)
    ρ0G = AlphaMainGrid(ρ0_G, u0_G, N)

    α0L = AlphaMainGrid(α0_L, u0_L, N)
    ρ0L = AlphaMainGrid(ρ0_L, u0_L, N)

    αxG = AlphaMainGrid(αx_G, ux_G, N)
    ρxG = AlphaMainGrid(ρx_G, ux_G, N)
    uxG = UVelMainGrid(ux_G, N)

    αxL = AlphaMainGrid(αx_L, ux_L, N)
    ρxL = AlphaMainGrid(ρx_L, ux_L, N)
    uxL = UVelMainGrid(ux_L, N)

    # Coeficientes da equação de momento
    D_G, DG = momentum_coefficients_for_pressure_equation(α0_G, ρ0_G, u0_G, N)
    D_L, DL = momentum_coefficients_for_pressure_equation(α0_L, ρ0_L, u0_G, N)

    # Matriz A
    ## Diagonal principal
    δP_D = zeros(N)
    δP_D[2:N-1] = @. (
        + DG.e*(αxG.e*ρxG.e)/ρG_ref
        + DG.w*(αxG.w*ρxG.w)/ρG_ref
        + DL.e*(αxL.e*ρxL.e)/ρL_ref
        + DL.w*(αxL.w*ρxL.w)/ρL_ref
        + (Δx/Δt)*(
            + αxG.P/ρG_ref * (1/CG^2)
            + αxL.P/ρL_ref * (1/CL^2)
            )
        )
    δP_D[1] = 1.0
    δP_D[N] = 1.0
    ## Diagonal superior
    δP_DU = zeros(N-1)
    δP_DU[2:N-1] = @. (
        - DG.e*(αxG.e*ρxG.e)/ρG_ref
        - DL.e*(αxL.e*ρxL.e)/ρL_ref
        )  
    δP_DU[1] = -1.0
    ## Diagonal inferior
    δP_DL = zeros(N-1)
    δP_DL[1:N-2] = @. (
        - DG.w*(αxG.w*ρxG.w)/ρG_ref
        - DL.w*(αxL.w*ρxL.w)/ρL_ref
        )
    δP_DL[N-1] = 0.0
    ## Construção da matriz A
    δP_A = Tridiagonal(δP_DL, δP_D, δP_DU)
    
    # Vetor b
    δP_b = zeros(N)
    δP_b[2:N-1] = @. (
        + (αxG.w*ρxG.w*uxG.w - αxG.e*ρxG.e*uxG.e)/ρG_ref
        + (αxL.w*ρxL.w*uxL.w - αxL.e*ρxL.e*uxL.e)/ρL_ref
        + (Δx/Δt)*(
            + (α0G.P*ρ0G.P - αxG.P*ρxG.P)/ρG_ref
            + (α0L.P*ρ0L.P - αxL.P*ρxL.P)/ρL_ref
            )
    )
    δP_b[1] = 0.0
    δP_b[N] = 0.0

    # Solução do sistema linear
    δP_x = δP_A \ δP_b
    
    # Correção dos valores
    ## Correção da pressão
    P = @. Px + δP_x
    ## Correção das Densidades
    ### Fase gasosa
    ρ_G = @. ρx_G + (1/CG^2)*δP_x
    ### Fase líquida
    ρ_L = @. ρx_L + (1/CL^2)*δP_x
    ## Correção das velocidades
    u_G = zeros(N+1)
    u_L = zeros(N+1)
    ### Fase gasosa
    u_G[1] = uxG.in
    u_G[2:N] = @. ux_G[2:N] + D_G*(δP_x[1:N-1] - δP_x[2:N])
    u_G[N+1] = u_G[N]
    ### Fase líquida
    u_L[1] = uxL.in
    u_L[2:N] = @. ux_L[2:N] + D_L*(δP_x[1:N-1] - δP_x[2:N])
    u_L[N+1] = u_L[N]

    return ρ_G, ρ_L, u_G, u_L, P
end

function momentum_coefficients_for_pressure_equation(
    α0_k::Vector{Float64},
    ρ0_k::Vector{Float64},
    u0_k::Vector{Float64},
    N :: Int64
    )

    α0k = AlphaSubGrid(α0_k, u0_k, N)
    ρ0k = RhoSubGrid(ρ0_k, u0_k, N)
    u0k = UVelSubGrid(u0_k, N)

    # Coeficientes pós-agrupamento (upwind)
    A_e = @. α0k.e/Δx
    a_e = @. α0k.e/Δt + abs(u0k.e)*α0k.e*ρ0k.e/Δx
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
#=
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
    αk_b[2:N-1] = @. a0_P*α0k.P/Δx
    αk_b[1] = α0k.in
    αk_b[N] = 0.0

    # Solução do sistema linear
    αk_x = αk_A \ αk_b

    return αk_x
end
=#