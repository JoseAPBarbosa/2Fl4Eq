using JSON


# Entrada de dados
dados = JSON.parsefile("Entrada.json")
## Dados
Comp    = dados["Comprimento [m]"]
D       = dados["Diâmetro [m]"]
θ       = dados["Inclinação [rad]"]
tempo   = dados["Tempo de simulação [s]"]
α_G_in  = dados["Fração de gás na entrada"]
α_G_i   = dados["Fração de gás inicial no duto"]
ρ_G     = dados["Densidade do gás [kg/m3]"]
ρ_L     = dados["Densidade do líquido [kg/m3]"]
u_G_in  = dados["Velocidade de entrada de gás [m/s]"]
u_L_in  = dados["Velocidade de entrada de líquido [m/s]"]
P_out   = dados["Pressão na saída do domínio [Pa]"]
N       = dados["Número de segmentos"]
CFL     = dados["CFL"]
Tol_P   = dados["Tolerância para a pressão"]
Tol_u   = dados["Tolerância para a velocidades"]
Tol_α   = dados["Tolerância para as frações"]

# Parâmetros
## Parâmetros físico
Area = (pi * D^2) / 4       # Área da seção do duto
α_L_in = 1 - α_G_in         # Fração de líquido na entrada
α_L_i = 1 - α_G_i           # Fração de líquido inicial
g = 9.81                    # Aceleração gravitacional
## Parâmetros de malha
Δx = Comp / N               # Passo espacial
Δt = (CFL * Δx) / u_L_in     # Passo temporal
Nt = Int(tempo ÷ Δt)        # Número de passos temporais
γ = 1.2                     # Garantia de hiperbolicidade

# Domínio das variáveis
## Variáveis em t=n
α0_G = α_G_i * ones(N)
α0_L = α_L_i * ones(N)
u0_G = u_G_in * ones(N + 1)
u0_L = u_L_in * ones(N + 1)
P0_ = P_out * ones(N)
## Variáveis em t=κ
αg_G = α_G_i * ones(N)
αg_L = α_L_i * ones(N)
ug_G = u_G_in * ones(N + 1)
ug_L = u_L_in * ones(N + 1)
Pg_ = P_out * ones(N)
## Variáveis em t=κ+1
α_G = α_G_i * ones(N)
α_L = α_L_i * ones(N)
u_G = u_G_in * ones(N + 1)
u_L = u_L_in * ones(N + 1)
P_ = P_out * ones(N)
## Erros
e_uL = 1
e_uG = 1
e_αL = 1
e_αG = 1
e_P = 1
