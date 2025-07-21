include("./interpolation_methods.jl")


const facemaingrid = NamedTuple{(:P, :E, :e),Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}
function face_maingrid(
    ϕk::Vector{Float64},
    N::Int64
    )
    # Função de construção do Tuple dos vetores das posições
    # dos valores dos centros das células nas equações de
    # conservação do momento

    P = ϕk[1:N-1]
    E = ϕk[2:N]
    e = @. (ϕk[1:N-1] + ϕk[2:N])/2

    return facemaingrid((P=P, E=E, e=e))
end

const facesubgrid = NamedTuple{(:P, :E, :e),Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}
function face_subgrid(
    ϕk::Vector{Float64},
    N::Int64
    )
    # Função de construção do Tuple dos vetores das posições
    # dos valores das faces das células nas equações de
    # conservação do momento

    P = @. (ϕk[1:N-1] + ϕk[2:N])/2
    E = @. (ϕk[2:N] + ϕk[3:N+1])/2
    e = ϕk[2:N]

    return facesubgrid((P=P, E=E, e=e))
end

const centermaingrid = NamedTuple{(:P, :w, :e),Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}
function center_maingrid(
    ϕk::Vector{Float64},
    N::Int64
    )
    # Função de construção do Tuple dos vetores das posições
    # dos valores dos centros das células nas equações de 
    # continuidade e de correção de pressão

    P = ϕk[2:N-1]
    w = @. (ϕk[1:N-2] + ϕk[2:N-1])/2
    e = @. (ϕk[2:N-1] + ϕk[3:N])/2

    return centermaingrid((P=P, w=w, e=e))
end

const centersubgrid = NamedTuple{(:w, :e),Tuple{Vector{Float64},Vector{Float64}}}
function center_subgrid(
    ϕk::Vector{Float64},
    N::Int64
    )
    # Função de construção do Tuple dos vetores das posições
    # dos valores das faces das células nas equações de 
    # continuidade e de correção de pressão

    w = ϕk[2:N-1]
    e = ϕk[3:N]

    return centersubgrid((w=w, e=e))
end
