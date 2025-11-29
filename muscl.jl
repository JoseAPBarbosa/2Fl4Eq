function momentum_muscl(
    ux_k :: Vector{Float64},
    u0_k :: Vector{Float64},
    α0k :: AlphaSubGrid,
    ρ0k :: RhoSubGrid,
    u0k :: UVelSubGrid,
    P :: PressSubGrid,
    CATHARE :: Vector{Float64},
    g :: Float64,
    θ :: Float64,
    N :: Int64,
    Δx :: Float64,
    Δt :: Float64
    )
    
    F = momentum_flux_flux(ux_k, u0_k, N, Δx, Δt)
    S = @. α0k.e*(P.P - P.E)/Δx + α0k.e*ρ0k.e*g*sin(θ) + CATHARE*(α0k.E - α0k.P)/Δx

    ux_k[1] = u0k.in
    ux_k[2:N] = @. u0k.e - (Δt/Δx)*(F.R - F.L) + (Δt/(α0k.e*ρ0k.e))*S
    ux_k[N+1] = ux_k[N]

    return ux_k
end

limiter(a, b) = a * b <= 0 ? 0 : (abs(a) < abs(b) ? a : b)      #minmod
function momentum_flux(
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
    s[1] =   limiter(0                     , (u_k[2] - u_k[1]/Δx))
    s[N+1] = limiter((u_k[N+1] - u_k[N]/Δx), 0                 )

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

function void_fraction_muscl(
    αx_k :: Vector{Float64},
    u0_k :: Vector{Float64},
    αxk :: AlphaSubGrid,
    ρxk :: RhoSubGrid,
    uxk :: UVelSubGrid,
    N :: Int64,
    Δx :: Float64,
    Δt :: Float64
    )
    


    F = momentum_flux_flux(ux_k, u0_k, N, Δx, Δt)
    S = @. α0k.e*(P.P - P.E)/Δx + α0k.e*ρ0k.e*g*sin(θ) + CATHARE*(α0k.E - α0k.P)/Δx

    ux_k[1] = u0k.in
    ux_k[2:N] = @. u0k.e - (Δt/Δx)*(F.R - F.L) + (Δt/(α0k.e*ρ0k.e))*S
    ux_k[N+1] = ux_k[N]

    return ux_k
end

function void_fraction_flux(
    α_k  :: Vector{Float64},
    α0_k :: Vector{Float64},
    N  :: Int64,
    Δx :: Float64,
    Δt :: Float64
    )

    # Reconstrução linear dos dados
    s = zeros(N)
    s_upwind   = @. (α_k[2:N-1] - α_k[1:N-2])/Δx
    s_downwind = @. (α_k[3:N]   - α_k[2:N-1])/Δx
    s[2:N-1] = @. limiter(s_upwind, s_downwind)
    s[1] = limiter(0                     , (α_k[2] - α_k[1]/Δx))
    s[N] = limiter((α_k[N] - α_k[N-1]/Δx), 0                   )

    αL_k = @. α_k - s*(Δx/2)
    αR_k = @. α_k + s*(Δx/2)
    
    # Evolucao dos valores (MUSCL-Hancock)
    ux_kP = @. (ux_k[1:N] + ux_k[2:N+1])/2

    AL_k = @. αL_k + (1/2)*(Δt/Δx)*ρx_G*ux_kP*(αL_k - αR_k)
    AR_k = @. αR_k + (1/2)*(Δt/Δx)*ρx_G*ux_kP*(αL_k - αR_k)
    
    #= Fluxo de Rusanov
    F_L = @. (1/2)*(u0_k[2:N] + abs(u0_k[2:N]))*UL_k[1:N-1] + (1/2)*(u0_k[2:N] - abs(u0_k[2:N]))*UR_k[2:N]
    F_R = @. (1/2)*(u0_k[2:N] + abs(u0_k[2:N]))*UL_k[2:N]   + (1/2)*(u0_k[2:N] - abs(u0_k[2:N]))*UR_k[3:N+1]
    
    return (L=F_L, R=F_R,)=#
end

void_fraction_flux()