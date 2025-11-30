using LinearAlgebra, JSON
using Printf, Plots, CSV, DataFrames

# Entrada de dados
dados = JSON.parsefile("Entrada.json")
## Dados
global Comp    = dados["Comprimento [m]"]
global D       = dados["Diâmetro [m]"]
global θ       = dados["Inclinação [rad]"]
global tempo   = dados["Tempo de simulação [s]"]
global α_G_in  = dados["Fração de gás na entrada"]
global α_G_i   = dados["Fração de gás inicial no duto"]
global ρ_G     = dados["Densidade do gás [kg/m3]"]
global ρ_L     = dados["Densidade do líquido [kg/m3]"]
global u_G_in  = dados["Velocidade de entrada de gás [m/s]"]
global u_L_in  = dados["Velocidade de entrada de líquido [m/s]"]
global P_out   = dados["Pressão na saída do domínio [Pa]"]
global CG      = dados["Velocidade do som do gás [m/s]"]
global CL      = dados["Velocidade do som do líquido [m/s]"]
global N       = dados["Número de segmentos"]
global CFL     = dados["CFL"]
global Tol_P   = dados["Tolerância para a pressão"]
global Tol_u   = dados["Tolerância para a velocidades"]
global Tol_α   = dados["Tolerância para as frações"]

# Parâmetros
## Parâmetros físico
Area = (pi * D^2) / 4       # Área da seção do duto
α_L_in = 1 - α_G_in         # Fração de líquido na entrada
α_L_i = 1 - α_G_i           # Fração de líquido inicial
global g = 9.81             # Aceleração gravitacional
## Parâmetros de malha
global Δx = Comp / N                # Passo espacial
global Δt = (CFL * Δx) / u_L_in     # Passo temporal
Nt = Int(tempo ÷ Δt)                # Número de passos temporais
global γ = 1.2                      # Garantia de hiperbolicidade

# Domínio das variáveis
## Variáveis em t=n
α0_G = fill(α_G_i, N)
ρ0_G = fill(ρ_G, N)
u0_G = fill(u_G_in, N+1)
α0_L = fill(α_L_i, N)
ρ0_L = fill(ρ_L, N)
u0_L = fill(u_L_in, N+1)
P0 = fill(P_out, N)
## Variáveis em t=κ
αx_G = fill(α_G_i, N)
ρx_G = fill(ρ_G, N)
ux_G = fill(u_G_in, N+1)
αx_L = fill(α_L_i, N)
ρx_L = fill(ρ_L, N)
ux_L = fill(u_L_in, N+1)
Px = fill(P_out, N)
## Variáveis em t=κ+1
α_G = fill(α_G_i, N)
ρ_G = fill(ρ_G, N)
u_G = fill(u_G_in, N+1)
α_L = fill(α_L_i, N)
ρ_L = fill(ρ_L, N)
u_L = fill(u_L_in, N+1)
P = fill(P_out, N)
## Erros
αGaux = zeros(N)
αLaux = zeros(N)
uGaux = zeros(N+1)
uLaux = zeros(N+1)
Paux = zeros(N)
e_uL = 1
e_uG = 1
e_αL = 1
e_αG = 1
e_P = 1
