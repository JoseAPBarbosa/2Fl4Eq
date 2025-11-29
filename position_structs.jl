include("./initialization.jl")
include("./interpolation_methods.jl")

#=============================================================================#
# Structs das equações da malha secundária                                    #
#=============================================================================#

struct AlphaSubGrid
    in :: Float64
    P :: SubArray{Float64, 1}
    E :: SubArray{Float64, 1}
    e :: Vector{Float64}

    function AlphaSubGrid(αk :: Vector{Float64}, uk :: Vector{Float64}, N :: Int64)
        in = αk[1]
        P = @view αk[1:N-1]
        E = @view αk[2:N]
        e = upwind_interpolation(αk, uk[2:N])
        new(in, P, E, e,)
    end
end

struct RhoSubGrid
    in :: Float64
    P :: SubArray{Float64, 1}
    E :: SubArray{Float64, 1}
    e :: Vector{Float64}

    function RhoSubGrid(ρk :: Vector{Float64}, uk :: Vector{Float64}, N :: Int64)
        in = ρk[1]
        P = @view ρk[1:N-1]
        E = @view ρk[2:N]
        e = upwind_interpolation(ρk, uk[2:N])
        new(in, P, E, e,)
    end
end

struct UVelSubGrid
    in :: Float64
    e :: SubArray{Float64, 1}

    function UVelSubGrid(uk :: Vector{Float64}, N :: Int64)
        in = uk[1]
        e = @view uk[2:N]
        new(in, e,)
    end
end

struct PressSubGrid
    out :: Float64
    P :: Vector{Float64}
    E :: Vector{Float64}

    function PressSubGrid(P_ :: Vector{Float64}, N :: Int64)
        out = P_[N]
        P = @view P_[1:N-1]
        E = @view P_[2:N]
        new(out, P, E)
    end
end


#=============================================================================#
# Structs das equações da malha primária                                      #
#=============================================================================#

struct AlphaCenter
    in :: Float64
    W :: SubArray{Float64, 1}
    P :: SubArray{Float64, 1}
    E :: SubArray{Float64, 1}
    w :: Vector{Float64}
    e :: Vector{Float64}

    function AlphaCenter(αk :: Vector{Float64}, uk :: Vector{Float64})
        N = length(αk)
        in = αk[1]
        W = @view αk[1:N-2]
        P = @view αk[2:N-1]
        E = @view αk[3:N]
        w = upwind_interpolation(αk[1:N-1], uk[2:N-1])
        e = upwind_interpolation(αk[2:N], uk[3:N])
        new(in, W, P, E, w, e)
    end
end

struct RhoCenter
    in :: Float64
    W :: SubArray{Float64, 1}
    P :: SubArray{Float64, 1}
    E :: SubArray{Float64, 1}
    w :: Vector{Float64}
    e :: Vector{Float64}

    function RhoCenter(ρk :: Vector{Float64}, uk :: Vector{Float64})
        N = length(ρk)
        in = ρk[1]
        W = @view ρk[1:N-2]
        P = @view ρk[2:N-1]
        E = @view ρk[3:N]
        w = upwind_interpolation(ρk[1:N-1], uk[2:N-1])
        e = upwind_interpolation(ρk[2:N], uk[3:N])
        new(in, W, P, E, w, e)
    end
end

struct UVelCenter
    in :: Float64
    w :: SubArray{Float64, 1}
    e :: SubArray{Float64, 1}

    function UVelCenter(uk :: Vector{Float64})
        N = length(uk)-1
        in = uk[1]
        e = @view uk[2:N-1]
        w = @view uk[3:N]
        new(in, w, e)
    end
end

struct PressCenter
    out :: Float64
    W :: SubArray{Float64, 1}
    P :: SubArray{Float64, 1}
    E :: SubArray{Float64, 1}

    function PressCenter(P_ :: Vector{Float64})
        N = length(P_)
        out = P_[N]
        W = @view P_[1:N-2]
        P = @view P_[2:N-1]
        E = @view P_[3:N]
        new(out, W, P, E)
    end
end
