include("./interpolation_methods.jl")


#=============================================================================#
# Structs das faces                                                           #
#=============================================================================#

struct AlphaFace
    in::Float64
    P::SubArray{Float64, 1}
    E::SubArray{Float64, 1}
    e::Vector{Float64}

    function AlphaFace(αk::Vector{Float64}, uk::Vector{Float64})
        N = length(αk)
        in = αk[1]
        P = @view αk[1:N-1]
        E = @view αk[2:N]
        e = upwind_interpolation(αk, uk[2:N])

        return new(in, P, E, e)
    end
end

struct RhoFace
    in::Float64
    P::SubArray{Float64, 1}
    E::SubArray{Float64, 1}
    e::Vector{Float64}

    function RhoFace(ρk::Vector{Float64}, uk::Vector{Float64})
        N = length(ρk)
        in = ρk[1]
        P = @view ρk[1:N-1]
        E = @view ρk[2:N]
        e = upwind_interpolation(ρk, uk[2:N])

        return new(in, P, E, e)
    end
end

struct UVelFace
    in::Float64
    P::Vector{Float64}
    E::Vector{Float64}
    e::SubArray{Float64, 1}

    function UVelFace(uk::Vector{Float64})
        N = length(uk)-1
        in = uk[1]
        P = upwind_interpolation(uk[1:N], (uk[1:N-1]+uk[2:N])/2)
        E = upwind_interpolation(uk[2:N+1], (uk[2:N]+uk[3:N+1])/2)
        e = @view uk[2:N]

        return new(in, P, E, e)
    end
end

struct PressFace
    out::Float64
    P::Vector{Float64}
    E::Vector{Float64}

    function PressFace(P_::Vector{Float64})
        N = length(P_)
        out = P_[1]
        P = @view P_[1:N-1]
        E = @view P_[2:N]

        return new(out, P, E)
    end
end


#=============================================================================#
# Structs dos centros                                                         #
#=============================================================================#

struct AlphaCenter
    in::Float64
    W::SubArray{Float64, 1}
    P::SubArray{Float64, 1}
    E::SubArray{Float64, 1}
    w::Vector{Float64}
    e::Vector{Float64}

    function AlphaCenter(αk::Vector{Float64}, uk::Vector{Float64})
        N = length(αk)
        in = αk[1]
        W = @view αk[1:N-2]
        P = @view αk[2:N-1]
        E = @view αk[3:N]
        w = upwind_interpolation(αk[1:N-1], uk[2:N-1])
        e = upwind_interpolation(αk[2:N], uk[3:N])

        return new(in, W, P, E, w, e)
    end
end

struct RhoCenter
    in::Float64
    W::SubArray{Float64, 1}
    P::SubArray{Float64, 1}
    E::SubArray{Float64, 1}
    w::Vector{Float64}
    e::Vector{Float64}

    function RhoCenter(ρk::Vector{Float64}, uk::Vector{Float64})
        N = length(ρk)
        in = ρk[1]
        W = @view ρk[1:N-2]
        P = @view ρk[2:N-1]
        E = @view ρk[3:N]
        w = upwind_interpolation(ρk[1:N-1], uk[2:N-1])
        e = upwind_interpolation(ρk[2:N], uk[3:N])

        return new(in, W, P, E, w, e)
    end
end

struct UVelCenter
    in::Float64
    w::SubArray{Float64, 1}
    e::SubArray{Float64, 1}

    function UVelCenter(uk::Vector{Float64})
        N = length(uk)-1
        in = uk[1]
        e = @view uk[2:N-1]
        w = @view uk[3:N]

        return new(in, w, e)
    end
end

struct PressCenter
    out::Float64
    W::SubArray{Float64, 1}
    P::SubArray{Float64, 1}
    E::SubArray{Float64, 1}

    function PressCenter(P_::Vector{Float64})
        N = length(P_)
        out = P_[1]
        W = @view P_[1:N-2]
        P = @view P_[2:N-1]
        E = @view P_[3:N]

        return new(out, W, P, E)
    end
end
