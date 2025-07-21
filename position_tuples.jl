include("./interpolation_methods.jl")


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

const facemaingrid = NamedTuple{(:P, :E, :e),Tuple{SubArray{Float64,1},SubArray{Float64,1},Vector{Float64}}}
function face_maingrid(
    ϕk::Vector{Float64},
    N::Int64
    )
    # Função de construção do Tuple dos vetores das posições
    # dos valores dos centros das células nas equações de
    # conservação do momento

    P = @view ϕk[1:N-1]
    E = @view ϕk[2:N]
    e = @. (ϕk[1:N-1] + ϕk[2:N])/2

    return facemaingrid((P=P, E=E, e=e))
end

const facesubgrid = NamedTuple{(:P, :E, :e),Tuple{Vector{Float64},Vector{Float64},SubArray{Float64,1}}}
function face_subgrid(
    ϕk::Vector{Float64},
    N::Int64
    )
    # Função de construção do Tuple dos vetores das posições
    # dos valores das faces das células nas equações de
    # conservação do momento

    P = @. (ϕk[1:N-1] + ϕk[2:N])/2
    E = @. (ϕk[2:N] + ϕk[3:N+1])/2
    e = @view ϕk[2:N]

    return facesubgrid((P=P, E=E, e=e))
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


#=

using JSON


# Entrada de dados
dados = JSON.parsefile("Entrada.json")
## Dados
Comp =  dados["Comprimento [m]"]
D =     dados["Diâmetro [m]"]
θ =     dados["Inclinação [rad]"]
tempo = dados["Tempo de simulação [s]"]
αG_in = dados["Fração de gás na entrada"]
αG_i =  dados["Fração de gás inicial no duto"]
ρG_i =  dados["Densidade do gás [kg/m3]"]
ρL_i =  dados["Densidade do líquido [kg/m3]"]
uG_in = dados["Velocidade de entrada de gás [m/s]"]
uL_in = dados["Velocidade de entrada de líquido [m/s]"]
P_out = dados["Pressão na saída do domínio [Pa]"]
cG =    dados["Velocidade do som do gás [m/s]"]
cL =    dados["Velocidade do som do líquido [m/s]"]
N =     dados["Número de segmentos"]
CFL =   dados["CFL"]
Tol_P = dados["Tolerância para a pressão"]
Tol_u = dados["Tolerância para a velocidades"]
Tol_α = dados["Tolerância para as frações"]

# Parâmetros
## Parâmetros físico
Area = (pi * D^2) / 4       # Área da seção do duto
αL_in = 1 - αG_in           # Fração de líquido na entrada
αL_i = 1 - αG_i             # Fração de líquido inicial
g = 9.81                    # Aceleração gravitacional
## Parâmetros de malha
Δx = Comp / N               # Passo espacial
Δt = (CFL * Δx) / uL_in     # Passo temporal
Nt = Int(tempo ÷ Δt)        # Número de passos temporais
γ = 1.2                     # Garantia de hiperbolicidade

# Variáveis
## Condições Iniciais
α_G = αG_i * ones(N)
α_L = αL_i * ones(N)
ρ_G = ρG_i * ones(N)
ρ_L = ρL_i * ones(N)
u_G = uG_in * ones(N + 1)
u_L = uL_in * ones(N + 1)
P_ = P_out * ones(N)
## Supostos
αØ_G = αG_i * ones(N)
αØ_L = αL_i * ones(N)
ρØ_G = ρG_i * ones(N)
ρØ_L = ρL_i * ones(N)
uØ_G = uG_in * ones(N + 1)
uØ_L = uL_in * ones(N + 1)
PØ_ = P_out * ones(N)
## Erros
e_uL = 1
e_uG = 1
e_αL = 1
e_αG = 1
e_P = 1


ωl = face_omega(uØ_L, N)
αl = face_maingrid(αØ_L, N)
ρl = face_maingrid(ρØ_L, N)
ul = face_subgrid(uØ_L, N)
=#