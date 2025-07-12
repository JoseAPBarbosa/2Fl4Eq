include("./tuples_de_posicao.jl")
using LinearAlgebra, JSON, Printf, Plots, CSV, DataFrames

# Entrada de dados
dados = JSON.parsefile("Entrada.json")

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
Area = (pi * D^2) / 4       # Área da seção do duto
αL_in = 1 - αG_in           # Fração de líquido na entrada
αL_i = 1 - αG_i             # Fração de líquido inicial
g = 9.81                    # aceleração gravitacional

Δx = Comp / N               # Passo espacial
Δt = (CFL * Δx) / uL_in     # Passo temporal
Nt = Int(tempo ÷ Δt)        # Número de passos temporais
γ = 1.2                     # coeficiente misterioso

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

function momentum_conservation_equation(
    αØ_g::Vector{Float64},
    ρØ_g::Vector{Float64},
    uØ_g::Vector{Float64},
    αØ_l::Vector{Float64},
    ρØ_l::Vector{Float64},
    uØ_l::Vector{Float64},
    PØ_::Vector{Float64}
)

    uØ_g_in = uØ_g[1]
    uØ_l_in = uØ_l[1]

    ωg = omega_face(uØ_g, N)
    αg = maingrid_face(αØ_g, ωg, N)
    ρg = maingrid_face(ρØ_g, ωg, N)
    ug = subgrid_face(uØ_g, N)
    ωl = omega_face(uØ_l, N)
    αl = maingrid_face(αØ_l, ωl, N)
    ρl = maingrid_face(ρØ_l, ωl, N)
    ul = subgrid_face(uØ_l, N)
    P = maingrid_face(PØ_, ωl, N)

    CATHARE = @. γ * (
        (αg.e*αl.e*ρg.e*ρl.e)*(ug.e - ul.e)^2 / (αg.e*ρl.e + αl.e*ρg.e)
    )

    uØ_g = momentum_linear_system(ωg, αg, ρg, ug, P, CATHARE, uØ_g_in)
    uØ_l = momentum_linear_system(ωl, αl, ρl, ul, P, CATHARE, uØ_l_in)
    
    return uØ_g, uØ_l
end

function momentum_linear_system(
    ωk::omegaface,
    αk::mainface,
    ρk::mainface,
    uk::subface,
    P::mainface,
    CATHARE::Vector{Float64},
    uk_in::Float64
    )
    
    # Matriz A
    ## Diagonal principal
    uk_D = @. (     # Sg,   2:N,    N-1
        (αk.e*ρk.e)/Δt + (ωk.E/2 + ωk.P/2)*(αk.e*ρk.e*uk.e)/Δx
        )
    uk_D_in = 1.0                                   # CC,   i=1,    1
    uk_D_out = 1.0                                  # CC,   i=N+1,  1
    uk_D = vcat([uk_D_in], uk_D, [uk_D_out])        # Sg,   1:N+1,  N+1    
    ## Diagonal superior
    uk_DU = @. (1/2 - ωk.E/2)*(αk.e*ρk.e*uk.e)/Δx   # Sg,   3:N+1,  N-1
    uk_DU_in = 0.0                                  # CC,   i=2,    1
    uk_DU = vcat([uk_DU_in], uk_DU)                 # Sg,   2:N+1,  N
    ## Diagonal inferior
    uk_DL = @. - (1/2 + ωk.P/2)*(αk.e*ρk.e*uk.e)/Δx     # Sg,   1:N-1   N-1
    uk_DL_out = -1.0                                    # CC,   i=N,    1
    uk_DL = vcat(uk_DL, [uk_DL_out])                    # Sg,   1:N,    N
    ## Construção da matriz A
    uk_A = Tridiagonal(uk_DL, uk_D, uk_DU)
    
    # Vetor b
    uk_b = @. (                         # Sg,   2:N,    N-1
        + (αk.e*ρk.e*uk.e)/Δt           # Termo u
        + αk.e*ρk.e*g*sin(θ)            # Termo g
        + αk.e*(P.P - P.E)/Δx           # Termo ΔP
        + CATHARE*(αk.E - αk.P)/Δx      # Termo ΔPi
    )
    uk_b_in = uk_in                                 # CC,   i=1,    1
    uk_b_out = 0.0                                  # CC,   i=N+1,  1
    uk_b = vcat([uk_b_in], uk_b, [uk_b_out])        # Sg,   1:N+1,  N+1
    
    # Solução do sistema linear
    uk_x = uk_A \ uk_b

    return uk_x
end

function pressure_correction_equation(
    ωØ_g::omegacenter,
    αØ_g::maincenter,
    ρØ_g::maincenter,
    uØ_g::subcenter,
    αØ_gn::maincenter,
    ρØ_gn::maincenter,
    ωØ_l::omegacenter,
    αØ_l::maincenter,
    ρØ_l::maincenter,
    uØ_l::subcenter,
    αØ_ln::maincenter,
    ρØ_ln::maincenter
    )
    
    AGe = @. αG.e/Δx
    ALe = @. αL.e/Δx
    aGe = @. (αG.e*ρG.e)/Δt + (ωG.E/2 + ωG.P/2)*(αG.e*ρG.e*uG.e)/Δx
    aLe = @. (αL.e*ρL.e)/Δt + (ωL.E/2 + ωL.P/2)*(αL.e*ρL.e*uL.e)/Δx

    AGw = @. αG.w/Δx
    ALw = @. αL.w/Δx
    aGw = @. (αG.w*ρG.w)/Δt + (ωG.P/2 + ωG.W/2)*(αG.w*ρG.w*uG.w)/Δx
    aLw = @. (αL.w*ρL.w)/Δt + (ωL.P/2 + ωL.W/2)*(αL.w*ρL.w*uL.w)/Δx

    # Matriz A
    ## Diagonal principal
    δP_D = @. (
        + (αG.e*ρG.e)/ρG.BC * (AGe/aGe)
        + (αG.w*ρG.w)/ρG.BC * (AGw/aGw)
        + (Δx/Δt)*(αG.P)/ρG.BC * (1/cG^2)
        + (αL.e*ρL.e)/ρL.BC * (ALe/aLe)
        + (αL.w*ρL.w)/ρL.BC * (ALw/aLw)
        + (Δx/Δt)*(αL.P)/ρL.BC * (1/cL^2)
    )
    δP_D_in = 1.0
    δP_D_out = 1.0
    δP_D = vcat([δP_D_in], δP_D, [δP_D_out])
    ## Diagonal superior
    δP_DU = @. (
        - (αG.e*ρG.e)/ρG.BC * (AGe/aGe)
        - (αL.e*ρL.e)/ρL.BC * (ALe/aLe)
    )  
    δP_DU_in = -1.0
    δP_DU = vcat([δP_DU_in], δP_DU)
    ## Diagonal inferior
    δP_DL = @. (
        - (αG.w*ρG.w)/ρG.BC * (AGw/aGw)
        - (αL.w*ρL.w)/ρL.BC * (ALw/aLw)
    )
    δP_DL_out = 0.0
    δP_DL = vcat(δP_DL, [δP_DL_out])
    ## Construção da matriz A
    δP_A = Tridiagonal(δP_DL, δP_D, δP_DU)

    # Vetor b
    δP_b = @. (
        + ((αG.w*ρG.w*uG.w) - (αG.e*ρG.e*uG.e))/ρG.BC
        + ((αL.w*ρL.w*uL.w) - (αL.e*ρL.e*uL.e))/ρL.BC
        + (Δx/Δt)*(
            + ((αG_n.P*ρG_n.P) - (αG.P*ρG.P))/ρG.BC
            + ((αL_n.P*ρL_n.P) - (αL.P*ρL.P))/ρL.BC
        )
    )
    δP_b_in = 0.0
    δP_b_out = 0.0
    δP_b = vcat([δP_b_in], δP_b, [δP_b_out])

    # Solução do sistema linear
    δP_x = δP_A \ δP_b

    # Correção dos valores
    ## Correção das massas específicas
    ### Fase gasosa
    ρØ_G = @. ρØ_G + (1/cG^2)*δP_x
    ### Fase líquida
    ρØ_L = @. ρØ_L + (1/cL^2)*δP_x

    ## Correção das velocidades
    ### Fase gasosa
    uØ_G[2:N] = @. uØ_G[2:N] + (A.G/a.G)*(δP_[1:N-1] - δP_[2:N]) #<--------
    uØ_G[N+1] = uØ_G[N]
    ### Fase líquida
    uØ_L[2:N] = @. uØ_L[2:N] + (A.G/a.G)*(δP_[1:N-1] - δP_[2:N]) #<----------
    uØ_L[N+1] = uØ_L[N]

    ## Correção da pressão
    PØ_ = @. PØ_ + δP_

    return δP_DL, δP_D, δP_DU, δP_B
end

function void_fraction_linear_system(
    ωk::omegacenter,
    αk_in::Float64,
    ρk::maincenter,
    uk::subcenter,
    αk_n::maincenter,
    ρk_n::maincenter
    )

    # Diagonal principal
    α_D = @. (
        + ρk.P/Δt
        - (1/2 + ωk.e/2)*(ρk.e*uk.e)/Δx
        + (1/2 - ωk.w/2)*(ρk.w*uk.w)/Δx
    )
    α_D_in = 1.0
    α_D_out = 1.0
    α_D = vcat([α_D_in], α_D, [α_D_out])

    # Diagonal superior
    α_DU = @. - (1/2 + ωk.e/2)*(ρk.e*uk.e)/Δx
    α_DU_in = 0.0
    α_DU = append!([α_DU_in], α_DU)

    # Diagonal inferior
    α_DL = @. (1/2 + ωk.w/2)*(ρk.w*uk.w)/Δx
    α_DL_out = -1.0
    α_DL = append!(α_DL, [α_DL_out])

    # Vetor B
    α_B = @. (αk_n.P*ρk_n.P)/Δt
    α_B_in = αk.BC
    α_B_out = 0.0
    α_B = append!([α_B_in], α_B, [α_B_out])

    return α_DL, α_D, α_DU, α_B
end

cont = 0
while (cont < 1) || (e_uL < Tol_u && e_uG < Tol_u && e_αL < Tol_α && e_αG < Tol_α && e_P < Tol_P)

    global α_G, α_L, ρ_G, ρ_L, u_G, u_L, P_
    global αØ_G, αØ_L, ρØ_G, ρØ_L, uØ_G, uØ_L, PØ_
    global cont, e_uL, e_uG, e_αL, e_αG, e_P
    global ρG_i, ρL_i, cG, cL
    local ωG, αØG, ρØG, uØG, ωL, αØL, ρØL, uØL, PØ

    # Equação do momento ----------------------------------------------------------------------------------------------------------------------------
    uØ_G, uØ_L = momentum_conservation_equation(αØ_G, ρØ_G, uØ_G, 
                                                αØ_L, ρØ_L, uØ_L, 
                                                PØ_)
    
    #= Equação de correção de pressão ----------------------------------------------------------------------------------------------------------------
    
    ## Atualização dos Tuples de posição
    ### Fase gasosa
    ωG = omega_center(uØ_G, N)
    αØG = maingrid_center(αØ_G, ωG, N)
    ρØG = maingrid_center(ρØ_G, ωG, N)
    uØG = subgrid_center(uØ_G, N)
    PØ = maingrid_center(PØ_, ωG, N)
    αØG_n = maingrid_center(α_G, ωG, N)   # valores em n
    ρØG_n = maingrid_center(ρ_G, ωG, N)   # valores em n
    ### Fase líquida
    ωL = omega_center(uØ_L, N)
    αØL = maingrid_center(αØ_L, ωL, N)
    ρØL = maingrid_center(ρØ_L, ωL, N)
    uØL = subgrid_center(uØ_L, N)
    PØ = maingrid_center(PØ_, ωL, N)
    αØL_n = maingrid_center(α_L, ωL, N)   # valores em n
    ρØL_n = maingrid_center(ρ_L, ωL, N)   # valores em n

    ## Sistema linear das equacoes de momento
    δP_DL, δP_D, δP_DU, δP_B = pressure_correction_linear_system(ωG, αØG, ρØG, uØG, αØG_n, ρØG_n, ωL, αØL, ρØL, uØL, αØL_n, ρØL_n)
    δP_A = Tridiagonal(δP_DL, δP_D, δP_DU)
    ### Solução do sistema linear
    δP_ = δP_A \ δP_B




    # Equação de conservação ------------------------------------------------------------------------------------------------------------------------
    
    ## Atualização dos Tuples de posição
    ### Fase gasosa
    ωG = omega_center(uØ_G, N)
    αØG = maingrid_center(αØ_G, ωG, N)
    ρØG = maingrid_center(ρØ_G, ωG, N)
    uØG = subgrid_center(uØ_G, N)
    PØ = maingrid_center(PØ_, ωG, N)
    αØG_n = maingrid_center(α_G, ωG, N)   # valores em n
    ρØG_n = maingrid_center(ρ_G, ωG, N)   # valores em n
    ### Fase líquida
    ωL = omega_center(uØ_L, N)
    αØL = maingrid_center(αØ_L, ωL, N)
    ρØL = maingrid_center(ρØ_L, ωL, N)
    uØL = subgrid_center(uØ_L, N)
    PØ = maingrid_center(PØ_, ωL, N)
    αØL_n = maingrid_center(α_L, ωL, N)   # valores em n
    ρØL_n = maingrid_center(ρ_L, ωL, N)   # valores em n
    
    ## Sistema linear das equacoes de continuidade
    ### Fase gasosa
    αG_DL, αG_D, αG_DU, αG_B = void_fraction_linear_system(ωG, αG_in, ρØG, uØG, αØG_n, ρØG_n)
    αG_A = Tridiagonal(αG_DL, αG_D, αG_DU)
    ### Fase líquida
    αL_DL, αL_D, αL_DU, αL_B = void_fraction_linear_system(ωL, αL_in, ρØL, uØL, αØL_n, ρØL_n)
    αL_A = Tridiagonal(αL_DL, αL_D, αL_DU)
    ### Solução do sistema linear
    αØ_G[:] = αG_A \ αG_B
    αØ_L[:] = αL_A \ αL_B
    =#

    # Atualização dos supostos ----------------------------------------------------------------------------------------------------------------------
#=
    e_uL = norm(u_L - uØ_L) / norm(uØ_L)
    e_uG = norm(u_G - uØ_G) / norm(uØ_G)
    e_αL = norm(α_L - αØ_L) / norm(αØ_L)
    e_αG = norm(α_G - αØ_G) / norm(αØ_G)
    e_P = norm(P_ - PØ_) / norm(PØ_)
=#
    cont = cont + 1

    #display(e_P)
    #display(cont)

end
