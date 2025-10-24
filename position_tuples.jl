const alphaFace = NamedTuple{(:P, :E, :e),Tuple{SubArray{Float64, 1},SubArray{Float64, 1},Vector{Float64}}}
function alpha_face(αk::Vector{Float64}, uk::Vector{Float64})

    N = length(αk)
    P = @view αk[1:N-1]
    E = @view αk[2:N]
    #e = @. (1+sign(uk[2:N]))/2 * αk[1:N-1] + (1-sign(uk[2:N]))/2 * αk[2:N]
    e = @. αk[1:N-1]/2 + αk[2:N]/2

    return alphaFace((P=P, E=E, e=e))
end

const rhoFace = NamedTuple{(:P, :E, :e),Tuple{SubArray{Float64, 1},SubArray{Float64, 1},Vector{Float64}}}
function rho_face(ρk::Vector{Float64}, uk::Vector{Float64})

    N = length(ρk)
    P = @view ρk[1:N-1]
    E = @view ρk[2:N]
    #e = @. (1+sign(uk[2:N]))/2 * ρk[1:N-1] + (1-sign(uk[2:N]))/2 * ρk[2:N]
    e = @. ρk[1:N-1]/2 + ρk[2:N]/2

    return rhoFace((P=P, E=E, e=e))
end

const uvelFace = NamedTuple{(:P, :E, :e),Tuple{Vector{Float64},Vector{Float64},SubArray{Float64, 1}}}
function uvel_face(uk::Vector{Float64})

    N = length(uk)-1
    P = @. uk[1:N-1]/2 + uk[2:N]/2
    E = @. uk[2:N]/2 + uk[3:N+1]/2
    e = @view uk[2:N]
    
    return uvelFace((P=P, E=E, e=e))
end

const pressFace = NamedTuple{(:P, :E),Tuple{SubArray{Float64, 1},SubArray{Float64, 1}}}
function press_face(ρk::Vector{Float64})

    N = length(ρk)
    P = @view ρk[1:N-1]
    E = @view ρk[2:N]

    return pressFace((P=P, E=E))
end

const facemaingrid = NamedTuple{(:P, :E, :e),Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}
function face_maingrid(
    ɸk::Vector{Float64},
    uk::Vector{Float64}
    )
    # Função de construção do Tuple dos vetores das posições
    # dos valores dos centros das células nas equações de
    # conservação do momento

    N = length(ɸk)
    P = ɸk[1:N-1]
    E = ɸk[2:N]
    e = @. (1+sign(uk[2:N]))/2 * ɸk[1:N-1] + (1-sign(uk[2:N]))/2 * ɸk[2:N]

    return facemaingrid((P=P, E=E, e=e))
end

const facesubgrid = NamedTuple{(:P, :E, :e),Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}
function face_subgrid(
    ɸk::Vector{Float64},
    uk::Vector{Float64}
    )
    # Função de construção do Tuple dos vetores das posições
    # dos valores das faces das células nas equações de
    # conservação do momento

    N = length(ɸk)-1
    P = @. (1+sign((uk[1:N-1]+uk[2:N])/2))/2 * ɸk[1:N-1] + (1-sign((uk[1:N-1]+uk[2:N])/2))/2 * ɸk[2:N]
    E = @. (1+sign((uk[2:N]+uk[3:N+1])/2))/2 * ɸk[2:N] + (1-sign((uk[2:N]+uk[3:N+1])/2))/2 * ɸk[3:N+1]
    e = ɸk[2:N]

    return facesubgrid((P=P, E=E, e=e))
end
 
const centermaingrid = NamedTuple{(:P, :w, :e),Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}
function center_maingrid(
    ɸk::Vector{Float64},
    N::Int64
    )
    # Função de construção do Tuple dos vetores das posições
    # dos valores dos centros das células nas equações de 
    # continuidade e de correção de pressão

    P = ɸk[2:N-1]
    w = @. (ɸk[1:N-2] + ɸk[2:N-1])/2
    e = @. (ɸk[2:N-1] + ɸk[3:N])/2

    return centermaingrid((P=P, w=w, e=e))
end

const centersubgrid = NamedTuple{(:w, :e),Tuple{Vector{Float64},Vector{Float64}}}
function center_subgrid(
    ɸk::Vector{Float64},
    N::Int64
    )
    # Função de construção do Tuple dos vetores das posições
    # dos valores das faces das células nas equações de 
    # continuidade e de correção de pressão

    w = ɸk[2:N-1]
    e = ɸk[3:N]

    return centersubgrid((w=w, e=e))
end
