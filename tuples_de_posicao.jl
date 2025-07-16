include("./interpoladores.jl")

const faceomega = NamedTuple{(:P, :E, :e),Tuple{SubArray{Float64,1},SubArray{Float64,1},SubArray{Float64,1}}}
function face_omega(
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

    return faceomega((P=P, E=E, e=e))
end

const facemaingrid = NamedTuple{(:P, :E, :e),Tuple{SubArray{Float64,1},SubArray{Float64,1},SubArray{Float64,1}}}
function face_maingrid(
    ϕk::Vector{Float64},
    ωk::faceomega,
    N::Int64
)
    # Função de construção do Tuple dos vetores das posições
    # dos valores dos centros das células nas equações de
    # conservação do momento

    P = @view ϕk[1:N-1]
    E = @view ϕk[2:N]
    ϕk_interp = upwind_interpolation(ϕk, ωk.e, N)
    e = @view ϕk_interp[2:N]

    return facemaingrid((P=P, E=E, e=e))
end

const facesubgrid = NamedTuple{(:w, :e, :ee),Tuple{SubArray{Float64,1},SubArray{Float64,1},SubArray{Float64,1}}}
function face_subgrid(
    ϕk::Vector{Float64},
    N::Int64
)
    # Função de construção do Tuple dos vetores das posições
    # dos valores das faces das células nas equações de
    # conservação do momento

    w = @view ϕk[1:N-1]
    e = @view ϕk[2:N]
    ee = @view ϕk[3:N+1]

    return facesubgrid((w=w, e=e, ee=ee))
end

const centeromega = NamedTuple{(:W, :P, :E, :w, :e),Tuple{SubArray{Float64,1},SubArray{Float64,1},SubArray{Float64,1},SubArray{Float64,1},SubArray{Float64,1}}}
function center_omega(
    uk::Vector{Float64},
    N::Int64
)
    # Função de construção do Tuple dos vetores do sinal do
    # valor do fluxo

    ωk = @. sign(uk)
    ωk_interp = @. sign((uk[1:N] + uk[2:N+1]) / 2)
    W = @view ωk_interp[1:N-2]
    P = @view ωk_interp[2:N-1]
    E = @view ωk_interp[3:N]
    w = @view ωk[2:N-1]
    e = @view ωk[3:N]

    return centeromega((W=W, P=P, E=E, w=w, e=e))
end

const centermaingrid = NamedTuple{(:P, :w, :e),Tuple{SubArray{Float64,1},SubArray{Float64,1},SubArray{Float64,1}}}
function center_maingrid(
    ϕk::Vector{Float64},
    ωk::centeromega,
    N::Int64
)
    # Função de construção do Tuple dos vetores das posições
    # dos valores dos centros das células nas equações de 
    # continuidade e de correção de pressão

    P = @view ϕk[2:N-1]
    ϕk_interp = upwind_interpolation(ϕk, ωk.P, N - 1)
    w = @view ϕk_interp[2:N-1]
    e = @view ϕk_interp[3:N]

    return centermaingrid((P=P, w=w, e=e))
end

const centersubgrid = NamedTuple{(:w, :e),Tuple{SubArray{Float64,1},SubArray{Float64,1}}}
function center_subgrid(
    ϕk::Vector{Float64},
    N::Int64
)
    # Função de construção do Tuple dos vetores das posições
    # dos valores das faces das células nas equações de 
    # continuidade e de correção de pressão

    w = @view ϕk[2:N-1]
    e = @view ϕk[3:N]

    return centersubgrid((w=w, e=e))
end
