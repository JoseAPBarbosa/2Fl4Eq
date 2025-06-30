include("./interpoladores.jl")

const omegaface = NamedTuple{(:P, :E, :e),Tuple{SubArray{Float64,1},SubArray{Float64,1},SubArray{Float64,1}}}
function omega_face(
    uk::Vector{Float64},
    N::Int64
)
    # Função de construção do Tuple dos vetores do sinal do
    # valor do fluxo

    ωk = @. sign(uk)
    ωk_interp = @. sign((uk[1:N] + uk[2:N+1]) / 2)
    P = @view ωk_interp[1:N-1]
    E = @view ωk_interp[2:N]
    e = @view ωk[2:N]

    return omegaface((P=P, E=E, e=e))
end

const mainface = NamedTuple{(:BC, :P, :E, :e),Tuple{Float64,SubArray{Float64,1},SubArray{Float64,1},SubArray{Float64,1}}}
function maingrid_face(
    ϕk::Vector{Float64},
    ωk::omegaface,
    N::Int64
)
    # Função de construção do Tuple dos vetores das posições
    # dos valores dos centros das células nas equações de
    # conservação do momento

    BC = ϕk[1]
    P = @view ϕk[1:N-1]
    E = @view ϕk[2:N]
    ϕk_interp = upwind_interpolation(ϕk, ωk.e, N)
    e = @view ϕk_interp[2:N]

    return mainface((BC=BC, P=P, E=E, e=e))
end

const subface = NamedTuple{(:BC, :w, :e, :ee),Tuple{Float64,SubArray{Float64,1},SubArray{Float64,1},SubArray{Float64,1}}}
function subgrid_face(
    ϕk::Vector{Float64},
    N::Int64
)
    # Função de construção do Tuple dos vetores das posições
    # dos valores das faces das células nas equações de
    # conservação do momento

    BC = ϕk[1]
    w = @view ϕk[1:N-1]
    e = @view ϕk[2:N]
    ee = @view ϕk[3:N+1]

    return subface((BC=BC, w=w, e=e, ee=ee))
end

const omegacenter = NamedTuple{(:W, :P, :E),Tuple{SubArray{Float64,1},SubArray{Float64,1},SubArray{Float64,1}}}
function omega_center(
    uk::Vector{Float64},
    N::Int64
)
    # Função de construção do Tuple dos vetores do sinal do
    # valor do fluxo

    ωk_interp = @. sign((uk[1:N] + uk[2:N+1]) / 2)
    W = @view ωk_interp[1:N-2]
    P = @view ωk_interp[2:N-1]
    E = @view ωk_interp[3:N]

    return omegacenter((W=W, P=P, E=E))
end

const maincenter = NamedTuple{(:BC, :P, :w, :e),Tuple{Float64,SubArray{Float64,1},SubArray{Float64,1},SubArray{Float64,1}}}
function maingrid_center(
    ϕk::Vector{Float64},
    ωk::omegacenter,
    N::Int64
)
    # Função de construção do Tuple dos vetores das posições
    # dos valores dos centros das células nas equações de 
    # continuidade e de correção de pressão

    BC = ϕk[1]
    P = @view ϕk[2:N-1]
    ϕk_interp = upwind_interpolation(ϕk, ωk.P, N - 1)
    e = @view ϕk_interp[3:N]
    w = @view ϕk_interp[2:N-1]

    return maincenter((BC=BC, P=P, w=w, e=e))
end

const subcenter = NamedTuple{(:BC, :w, :e),Tuple{Float64,SubArray{Float64,1},SubArray{Float64,1}}}
function subgrid_center(
    ϕk::Vector{Float64},
    N::Int64
)
    # Função de construção do Tuple dos vetores das posições
    # dos valores das faces das células nas equações de 
    # continuidade e de correção de pressão

    BC = ϕk[1]
    w = @view ϕk[2:N-1]
    e = @view ϕk[3:N]

    return subcenter((BC=BC, w=w, e=e))
end
