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
    ρ_L::Float64
    )

    uG_in = u0_G[1]
    uL_in = u0_L[1]  

    α0G = alpha_face(α0_G, ux_G)
    α0L = alpha_face(α0_L, ux_L)
    u0G = uvel_face(u0_G)
    u0L = uvel_face(u0_L)

    αxG = alpha_face(αx_G, ux_G)
    αxL = alpha_face(αx_L, ux_L)
    uxG = uvel_face(ux_G)
    uxL = uvel_face(ux_L)
    Px = press_face(Px_)

    ux_G = momentum_linear_system(α0G, u0G, αxG, uxG, Px, ρ_G, uG_in)
    ux_L = momentum_linear_system(α0L, u0L, αxL, uxL, Px, ρ_L, uL_in)
    
    return ux_G, ux_L
end

function momentum_linear_system(
    α0k::alphaFace,
    u0k::uvelFace,
    αxk::alphaFace,
    uxk::uvelFace,
    Px::pressFace,
    ρ_k::Float64,
    uk_in::Float64
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
    uk_DL[1:N-1] = -a_w
    uk_DL[N] = -1.0
    ## Construção da matriz A
    uk_A = Tridiagonal(uk_DL, uk_D, uk_DU)
    
    # Vetor b
    uk_b = zeros(N+1)
    uk_b[2:N] = @. (
        + a0_e*u0k.e
        + A_e*(Px.P - Px.E)
        + ρ_k*αxk.e*g*sin(θ)
    )
    uk_b[1] = uk_in
    uk_b[N+1] = 0.0

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
    ρ_L::Float64
    )

    d_G = Δt/(ρ_G*Δx)
    d_L = Δt/(ρ_L*Δx)
    #A_Ge, a_Ge, AG, aG = momentum_coefficients_for_pressure_equation(α0_G, αx_G, ux_G, ρ_G)
    #A_Le, a_Le, AL, aL = momentum_coefficients_for_pressure_equation(α0_L, αx_L, ux_L, ρ_L)

    α0G = alpha_center(α0_G, ux_G)
    α0L = alpha_center(α0_L, ux_L)
    αxG = alpha_center(αx_G, ux_G)
    αxL = alpha_center(αx_L, ux_L)
    uxG = uvel_center(ux_G)
    uxL = uvel_center(ux_L)

    # Matriz A
    ## Diagonal principal
    δP_D = zeros(N)
    δP_D[2:N-1] = @. (
        + αxG.e * d_G
        + αxG.w * d_G
        + αxL.e * d_L
        + αxL.w * d_L
    )
    δP_D[1] = 1.0
    δP_D[N] = 1.0
    ## Diagonal superior
    δP_DU = zeros(N-1)
    δP_DU[2:N-1] = @. (
        - αxG.e * d_G
        - αxL.e * d_L
    )  
    δP_DU[1] = -1.0
    ## Diagonal inferior
    δP_DL = zeros(N-1)
    δP_DL[1:N-2] = @. (
        - αxG.w * d_G
        - αxL.w * d_L
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
    u_G[2:N] = @. ux_G[2:N] + d_G*(δP_x[1:N-1] - δP_x[2:N])
    u_G[N+1] = u_G[N]
    u_L[2:N] = @. ux_L[2:N] + d_L*(δP_x[1:N-1] - δP_x[2:N])
    u_L[N+1] = u_L[N]
    ## Correção da pressão
    P_ = @. Px_ + δP_x
    
    return u_G, u_L, P_
end

#=
function momentum_coefficients_for_pressure_equation(
    α0k::alphaface,
    αxk::alphaface,
    uxk::uvelFace,
    ρ_k::Float64
    )

    α0k = alpha_face(α0k, uxk)
    αxk = alpha_face(αxk, uxk)
    uxk = uvel_face(uxk)

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

    # Tuples de posição
    Ak = @views (
        w = A_e[2:N-1],
        e = A_e[1:N-2],
    )
    ak = @views (
        w = a_e[2:N-1],
        e = a_e[1:N-2],
    )

    return A_ke, a_ke, Ak, ak
end
=#
#=
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

    αx_G = void_fraction_linear_system(α0G, uG, ρ_G, αG_in)
    αx_L = void_fraction_linear_system(α0L, uL, ρ_L, αL_in)

    α_L[2:N] = @. αx_L[2:N] / (αx_L[2:N] + αx_G[2:N])
    α_G[2:N] = @. 1.0 - α_L[2:N]

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
=#