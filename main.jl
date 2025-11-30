### Inicialização ###
include("./initialization.jl")
using Plots

### Processamento ####
include("./equation_solvers.jl")

cont = 0
for t = 1:Nt
    global α0_G, ρ0_G, u0_G, α0_L, ρ0_L, u0_L, P0
    global αx_G, ρx_G, ux_G, αx_L, ρx_L, ux_L, Px
    global α_G, ρ_G, u_G, α_L, ρ_L, u_L, P
    global αGaux, αLaux, uGaux, uLaux, Paux, cont
    local e_αG, e_αL, e_uG, e_uL, e_P

    e_αG = 1
    e_αL = 1
    e_uG = 1
    e_uL = 1
    e_P = 1

    while e_αG > Tol_α || e_αL > Tol_α || e_uG > Tol_u || e_uL > Tol_u || e_P > Tol_P
    
        αGaux[:] = αx_G[:]
        αLaux[:] = αx_L[:]
        uGaux[:] = ux_G[:]
        uLaux[:] = ux_L[:]
        Paux[:] = Px[:]

        # Equação do momento
        ux_G, ux_L = guessed_velocity_equation_solver(α0_G, ρ0_G, u0_G, α0_L, ρ0_L, u0_L, ux_G, ux_L, Px, g, θ, γ, N, Δx, Δt)

        # Equação de correção de pressão
        ρ_G, ρ_L, u_G, u_L, P = pressure_correction_equation_solver(α0_G, ρ0_G, α0_L, ρ0_L, αx_G, ρx_G, ux_G, αx_L, ρx_L, ux_L, Px, CG, CL, N, Δx, Δt)
        
        # Equação de conservação
        α_G, α_L = continuity_equation_solver(αx_G, ρx_G, ux_G, αx_L, ρx_L, ux_L, ρ_G, ρ_L, N, Δx, Δt)

        # Cálculo dos erros
        e_αG = norm(α_G - αx_G) / norm(αGaux)
        e_αL = norm(α_L - αx_L) / norm(αLaux)
        e_uG = norm(u_G - ux_G) / norm(uGaux)
        e_uL = norm(u_L - ux_L) / norm(uLaux)
        e_P = norm(P - Px) / norm(Paux)
        
        # Atualização dos supostos
        αx_G[:] = α_G[:]
        αx_L[:] = α_L[:]
        ux_G[:] = u_G[:]
        ux_L[:] = u_L[:]
        Px[:] = P[:]
        cont = cont + 1
        print(cont)
        print(" | ")
        print(norm(α_G - αx_G))
        print(" | ")
        print(norm(α_L - αx_L))
        print(" | ")
        print(norm(u_G - ux_G))
        print(" | ")
        println(norm(u_L - ux_L))        

    end

    # Atualização Temporal -------------------------------------------------
    α0_G[:] = α_G[:]
    α0_L[:] = α_L[:]
    u0_G[:] = u_G[:]
    u0_L[:] = u_L[:]
    P0[:] = P[:]

end
println("ok")