include("./interpolation_methods.jl")


#=============================================================================#
# Tuple das Faces                                                             #
#=============================================================================#

const alphaFace = NamedTuple{(:in, :P, :E, :e,),Tuple{Float64,SubArray{Float64, 1},SubArray{Float64, 1},Vector{Float64}}}
function alpha_face(αk::Vector{Float64}, uk::Vector{Float64})

    N = length(αk)
    in = αk[1]
    P = @view αk[1:N-1]
    E = @view αk[2:N]
    e = upwind_interpolation(αk, uk[2:N])

    return alphaFace((in=in, P=P, E=E, e=e,))
end

const rhoFace = NamedTuple{(:in, :P, :E, :e),Tuple{Float64,SubArray{Float64, 1},SubArray{Float64, 1},Vector{Float64}}}
function rho_face(ρk::Vector{Float64}, uk::Vector{Float64})

    N = length(ρk)
    in = ρk[1]
    P = @view ρk[1:N-1]
    E = @view ρk[2:N]
    e = upwind_interpolation(ρk, uk[2:N])

    return rhoFace((in=in, P=P, E=E, e=e,))
end

const uvelFace = NamedTuple{(:in, :P, :E, :e),Tuple{Float64,Vector{Float64},Vector{Float64},SubArray{Float64, 1}}}
function uvel_face(uk::Vector{Float64})

    N = length(uk)-1
    in = uk[1]
    P = upwind_interpolation(uk[1:N], (uk[1:N-1]+uk[2:N])/2)
    E = upwind_interpolation(uk[2:N+1], (uk[2:N]+uk[3:N+1])/2)
    e = @view uk[2:N]

    return uvelFace((in=in, P=P, E=E, e=e,))
end

const pressFace = NamedTuple{(:out, :P, :E,),Tuple{Float64,SubArray{Float64, 1},SubArray{Float64, 1}}}
function press_face(P_::Vector{Float64})

    N = length(P_)
    out = P_[1]
    P = @view P_[1:N-1]
    E = @view P_[2:N]

    return pressFace((out=out, P=P, E=E,))
end


#=============================================================================#
# Tuple dos Centros                                                           #
#=============================================================================#


const alphaCenter = NamedTuple{(:in, :W, :P, :E, :w, :e,),Tuple{Float64,SubArray{Float64, 1},SubArray{Float64, 1},SubArray{Float64, 1},Vector{Float64},Vector{Float64}}}
function alpha_center(αk::Vector{Float64}, uk::Vector{Float64})

    N = length(αk)
    in = αk[1]
    W = @view αk[1:N-2]
    P = @view αk[2:N-1]
    E = @view αk[3:N]
    w = upwind_interpolation(αk[1:N-1], uk[2:N-1])
    e = upwind_interpolation(αk[2:N], uk[3:N])

    return alphaCenter((in=in, W=W, P=P, E=E, w=w, e=e,))
end

const rhoCenter = NamedTuple{(:in, :W, :P, :E, :w, :e,),Tuple{Float64,SubArray{Float64, 1},SubArray{Float64, 1},SubArray{Float64, 1},Vector{Float64},Vector{Float64}}}
function rho_center(ρk::Vector{Float64}, uk::Vector{Float64})

    N = length(ρk)
    in = ρk[1]
    W = @view ρk[1:N-2]
    P = @view ρk[2:N-1]
    E = @view ρk[3:N]
    w = upwind_interpolation(ρk[1:N-1], uk[2:N-1])
    e = upwind_interpolation(ρk[2:N], uk[3:N])

    return alphaCenter((in=in, W=W, P=P, E=E, w=w, e=e,))
end

const uvelCenter = NamedTuple{(:in, :P, :E, :e),Tuple{Float64,Vector{Float64},Vector{Float64},SubArray{Float64, 1}}}
function uvel_center(uk::Vector{Float64})

    N = length(uk)-1
    in = uk[1]
    e = @view uk[2:N-1]
    w = @view uk[3:N]

    return uvelCenter((in=in, w=w, e=e,))
end

const pressCenter = NamedTuple{(:out, :W, :P, :E,),Tuple{Float64,SubArray{Float64, 1},SubArray{Float64, 1},SubArray{Float64, 1}}}
function press_center(P_::Vector{Float64})

    N = length(P_)
    out = P_[1]
    W = @view αk[1:N-2]
    P = @view αk[2:N-1]
    E = @view αk[3:N]

    return pressCenter((out=out, W=W, P=P, E=E,))
end
