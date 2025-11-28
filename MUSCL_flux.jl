flux = u -> 1/2 * u^2       # F(u)
wave_speed = u -> u         # dF(u)/du
minmod(a, b) = a * b <= 0 ? 0 : (abs(a) < abs(b) ? a : b)

function MUSCL_fluxes(qi :: Vector{Float64}, Δx :: Float64, Δt :: Float64)
    N = length(qi)

    # reconstrução dos dados
    si = zeros(N)
    s_upwind = @. (qi[2:N-1] - qi[1:N-2])/Δx
    s_downwind = @. (qi[3:N] - qi[2:N-1])/Δx
    si[2:N-1] = @. minmod(s_upwind, s_downwind)
    si[1] = minmod(0, (qi[2] - qi[1]/Δx))
    si[N] = minmod((qi[N] - qi[N-1]/Δx), 0)

    qiL = @. qi - si*(Δx/2)
    qiR = @. qi + si*(Δx/2)
    
    # evolucao dos valores (MUSCL-Hancock)
    qiL_barra = @. qiL + (1/2)*(Δt/Δx)*(flux(qiL) - flux(qiR))
    qiR_barra = @. qiR + (1/2)*(Δt/Δx)*(flux(qiL) - flux(qiR))

    # fluxo de Rusanov
    s_max = @. max(wave_speed(qiL_barra[1:N-1]), wave_speed(qiR_barra[2:N]))
    F_L = @. (1/2)*(flux(qiL_barra[1:N-2]) + flux(qiR_barra[2:N-1])) - (1/2)*s_max[1:N-2]*(qiR_barra[2:N-1] - qiL_barra[1:N-2])
    F_R = @. (1/2)*(flux(qiL_barra[2:N-1]) + flux(qiR_barra[3:N]))   - (1/2)*s_max[2:N-1]*(qiR_barra[3:N]   - qiL_barra[2:N-1])
    
    return F_L, F_R
end
