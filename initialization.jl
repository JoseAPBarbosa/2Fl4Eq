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
αø_G = αG_i * ones(N)
αø_L = αL_i * ones(N)
ρø_G = ρG_i * ones(N)
ρø_L = ρL_i * ones(N)
uø_G = uG_in * ones(N + 1)
uø_L = uL_in * ones(N + 1)
Pø_ = P_out * ones(N)
## Erros
e_uL = 1
e_uG = 1
e_αL = 1
e_αG = 1
e_P = 1
