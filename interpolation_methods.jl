function upwind_interpolation(
    ɸk::Vector{Float64},
    uk::Vector{Float64}
    )

    ɸk_e = @. (1+sign(uk))/2 * ɸk[1:end-1] + (1-sign(uk))/2 * ɸk[2:end]

    return ɸk_e
end
