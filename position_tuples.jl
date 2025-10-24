#==========================#
# Face tuples
#==========================#
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
function press_face(p::Vector{Float64})

    N = length(p)
    P = @view p[1:N-1]
    E = @view p[2:N]

    return pressFace((P=P, E=E))
end

#==========================#
# Center tuples
#==========================#

const alphaCenter = NamedTuple{(:W, :P, :E, :w, :e),Tuple{SubArray{Float64, 1},SubArray{Float64, 1},SubArray{Float64, 1},Vector{Float64},Vector{Float64}}}
function alpha_center(
    αk::Vector{Float64},
    uk::Vector{Float64}
    )

    N = length(αk)
    W = @view αk[1:N-2]
    P = @view αk[2:N-1]
    E = @view αk[3:N]
    w = @. (1+sign(uk[2:N-1]))/2 * αk[1:N-2] + (1-sign(uk[2:N-1]))/2 * αk[2:N-1]
    e = @. (1+sign(uk[3:N]))/2  * αk[2:N-1]  + (1-sign(uk[3:N]))/2  * αk[3:N]

    return alphaCenter((W=W, P=P, E=E, w=w, e=e))
end

const rhoCenter = NamedTuple{(:W, :P, :E, :w, :e),Tuple{SubArray{Float64, 1},SubArray{Float64, 1},SubArray{Float64, 1},Vector{Float64},Vector{Float64}}}
function rho_center(
    ρk::Vector{Float64},
    uk::Vector{Float64}
    )

    N = length(ρk)
    W = @view ρk[1:N-2]
    P = @view ρk[2:N-1]
    E = @view ρk[3:N]
    w = @. (1+sign(uk[2:N-1]))/2 * ρk[1:N-2] + (1-sign(uk[2:N-1]))/2 * ρk[2:N-1]
    e = @. (1+sign(uk[3:N]))/2  * ρk[2:N-1] + (1-sign(uk[3:N]))/2  * ρk[3:N]

    return rhoCenter((W=W, P=P, E=E, w=w, e=e))
end

const uvelCenter = NamedTuple{(:w, :e),Tuple{SubArray{Float64, 1},SubArray{Float64, 1}}}
function uvel_center(uk::Vector{Float64})

    N = length(uk)-1
    w = @view uk[2:N-1]
    e = @view uk[3:N]
    
    return uvelCenter((w=w, e=e))
end
