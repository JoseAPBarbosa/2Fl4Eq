function upwind_interpolation(
    ϕk::Vector{Float64},
    ωk::SubArray{Float64, 1},
    N::Int64
    )
    # Função de interpolação upwind

    ϕk_e = @. (1/2 + ωk/2) * ϕk[1:N-1] + (1/2 - ωk/2) * ϕk[2:N]     # Sg,   2:N,    N-1
    ϕk_in = ϕk[1]                           # CC,   i=1,    1
    ϕk_out = 2 * ϕk_e[N-1] - ϕk_e[N-2]      # CC,   i=N+1,  1
    ϕk_e = vcat([ϕk_in], ϕk_e, [ϕk_out])    # Sg,   1:N+1,  N+1

    return ϕk_e     # Sg,   1:N+1,  N+1
end
