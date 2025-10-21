function face_upwind_interpolation(
    ɸk::Vector{Float64},
    uk::Vector{Float64}
    )

    N = length(ɸk)
    ɸk_e[2:N] = @. (1+sign(uk[2:N]))/2 * ɸk[1:N-1] + (1-sign(uk[2:N]))/2 * ɸk[2:N]

    return ɸk_e
end
