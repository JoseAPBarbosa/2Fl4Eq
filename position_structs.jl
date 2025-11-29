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
        e = @. (1+sign(uk[2:N]))/2 * αk[1:N-1] + (1-sign(uk[2:N]))/2 * αk[2:N]
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
        e = @. (1+sign(uk[2:N]))/2 * ρk[1:N-1] + (1-sign(uk[2:N]))/2 * ρk[2:N]
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

struct AlphaMainGrid
    in :: Float64
    W :: SubArray{Float64, 1}
    P :: SubArray{Float64, 1}
    E :: SubArray{Float64, 1}
    w :: Vector{Float64}
    e :: Vector{Float64}

    function AlphaMainGrid(αk :: Vector{Float64}, uk :: Vector{Float64}, N :: Int64)
        in = αk[1]
        W = @view αk[1:N-2]
        P = @view αk[2:N-1]
        E = @view αk[3:N]
        w = @. (1+sign(uk[2:N-1]))/2 * αk[1:N-2] + (1-sign(uk[2:N-1]))/2 * αk[2:N-1]
        e = @. (1+sign(uk[3:N]))/2   * αk[2:N-1] + (1-sign(uk[3:N]))/2   * αk[3:N]
        new(in, W, P, E, w, e)
    end
end

struct RhoMainGrid
    in :: Float64
    W :: SubArray{Float64, 1}
    P :: SubArray{Float64, 1}
    E :: SubArray{Float64, 1}
    w :: Vector{Float64}
    e :: Vector{Float64}

    function RhoMainGrid(ρk :: Vector{Float64}, uk :: Vector{Float64}, N :: Int64)
        in = ρk[1]
        W = @view ρk[1:N-2]
        P = @view ρk[2:N-1]
        E = @view ρk[3:N]
        w = @. (1+sign(uk[2:N-1]))/2 * ρk[1:N-2] + (1-sign(uk[2:N-1]))/2 * ρk[2:N-1]
        e = @. (1+sign(uk[3:N]))/2   * ρk[2:N-1] + (1-sign(uk[3:N]))/2   * ρk[3:N]
        new(in, W, P, E, w, e)
    end
end

struct UVelMainGrid
    in :: Float64
    w :: SubArray{Float64, 1}
    e :: SubArray{Float64, 1}

    function UVelMainGrid(uk :: Vector{Float64}, N :: Int64)
        in = uk[1]
        e = @view uk[2:N-1]
        w = @view uk[3:N]
        new(in, w, e)
    end
end

struct PressMainGrid
    out :: Float64
    W :: SubArray{Float64, 1}
    P :: SubArray{Float64, 1}
    E :: SubArray{Float64, 1}

    function PressMainGrid(P_ :: Vector{Float64}, N :: Int64)
        out = P_[N]
        W = @view P_[1:N-2]
        P = @view P_[2:N-1]
        E = @view P_[3:N]
        new(out, W, P, E)
    end
end
