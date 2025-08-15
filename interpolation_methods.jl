function upwind_interpolation(
    ɸk::AbstractVector{Float64},
    ωk::SubArray{Float64, 1},
    N::Int64
    )

    ɸk_e = @. (1/2 + ωk/2) * ɸk[1:N-1] + (1/2 - ωk/2) * ɸk[2:N]     # Sg,   2:N,    N-1
    ɸk_in = ɸk[1]                           # CC,   i=1,    1
    ɸk_out = 2 * ɸk_e[N-1] - ɸk_e[N-2]      # CC,   i=N+1,  1
    ɸk_e = vcat([ɸk_in], ɸk_e, [ɸk_out])    # Sg,   1:N+1,  N+1

    return ɸk_e     # Sg,   1:N+1,  N+1
end
