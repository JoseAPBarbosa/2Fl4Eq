include("./initialization.jl")
include("./position_structs.jl")
using Plots

function guessed_velocity_equation_solver(
    α0_G::Vector{Float64},
    ρ0_G::Vector{Float64},
    u0_G::Vector{Float64},
    α0_L::Vector{Float64},
    ρ0_L::Vector{Float64},
    u0_L::Vector{Float64},
    αx_G::Vector{Float64},
    ρx_G::Vector{Float64},
    ux_G::Vector{Float64},
    αx_L::Vector{Float64},
    ρx_L::Vector{Float64},
    ux_L::Vector{Float64},
    Px ::Vector{Float64},
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

    ux_G = MUSCL(ux_G, u0_G, α0G, ρ0G, u0G, P, CATHARE, g, θ, N, Δx, Δt)
    ux_L = MUSCL(ux_L, u0_L, α0L, ρ0L, u0L, P, CATHARE, g, θ, N, Δx, Δt)

    return ux_G, ux_L
end

function MUSCL(
    ux_k :: Vector{Float64},
    u0_k :: Vector{Float64},
    α0k :: AlphaSubGrid,
    ρ0k :: RhoSubGrid,
    u0k :: UVelSubGrid,
    P  :: PressSubGrid,
    CATHARE :: Vector{Float64},
    g  :: Float64,
    θ  :: Float64,
    N  :: Int64,
    Δx :: Float64,
    Δt :: Float64
    )
    
    F = MUSCL_flux(ux_k, u0_k, N, Δx, Δt)
    S = @. α0k.e*(P.P - P.E)/Δx + α0k.e*ρ0k.e*g*sin(θ) + CATHARE*(α0k.E - α0k.P)/Δx

    ux_k[1] = u0k.in
    ux_k[2:N] = @. u0k.e - (Δt/Δx)*(F.R - F.L) + (Δt/(α0k.e*ρ0k.e))*(S)
    ux_k[N+1] = ux_k[N]

    return ux_k
end

limiter(a, b) = a * b <= 0 ? 0 : (abs(a) < abs(b) ? a : b)      #minmod
function MUSCL_flux(
    u_k  :: Vector{Float64},
    u0_k :: Vector{Float64},
    N  :: Int64,
    Δx :: Float64,
    Δt :: Float64
    )
    
    # Reconstrução linear dos dados
    s = zeros(N+1)
    s_upwind   = @. (u_k[2:N]   - u_k[1:N-1])/Δx
    s_downwind = @. (u_k[3:N+1] - u_k[2:N]  )/Δx
    s[2:N] = @. limiter(s_upwind, s_downwind)
    s[1] = limiter(0                     , (u_k[2] - u_k[1]/Δx))
    s[N] = limiter((u_k[N+1] - u_k[N]/Δx), 0                   )

    uL_k = @. u_k - s*(Δx/2)
    uR_k = @. u_k + s*(Δx/2)
    
    # Evolucao dos valores (MUSCL-Hancock)
    UL_k = @. uL_k + (1/2)*(Δt/Δx)*u0_k*(uL_k - uR_k)
    UR_k = @. uR_k + (1/2)*(Δt/Δx)*u0_k*(uL_k - uR_k)
    
    # Fluxo de Rusanov
    F_L = @. (1/2)*(u0_k[2:N] + abs(u0_k[2:N]))*UL_k[1:N-1] + (1/2)*(u0_k[2:N] - abs(u0_k[2:N]))*UR_k[2:N]
    F_R = @. (1/2)*(u0_k[2:N] + abs(u0_k[2:N]))*UL_k[2:N]   + (1/2)*(u0_k[2:N] - abs(u0_k[2:N]))*UR_k[3:N+1]

    return (L=F_L, R=F_R,)
end

a, b = guessed_velocity_equation_solver(α0_G, ρ0_G, u0_G, α0_L, ρ0_L, u0_L, αx_G, ρx_G, ux_G, αx_L, ρx_L, ux_L, Px_, g, θ, γ, N, Δx, Δt)
