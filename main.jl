### Inicialização ###
include("./initialization.jl")


### Processamento ####
include("./equation_solvers.jl")

t = 0
while t < 1
    global t
    #cont = 0
    while e_uL > Tol_u || e_uG > Tol_u || e_αL > Tol_α || e_αG > Tol_α || e_P > Tol_P
        global α0_G, α0_L, u0_G, u0_L, P0_
        global αx_G, αx_L, ux_G, ux_L, Px_
        global α_G, α_L, u_G, u_L, P_
        global N, ρ_G, ρ_L, e_uL, e_uG, e_αL, e_αG, e_P, cont

        # Equação do momento
        ux_G, ux_L = momentum_conservation_equation_solver( α0_G, α0_L, u0_G, u0_L, 
                                                            αx_G, αx_L, ux_G, ux_L, 
                                                            Px_, ρ_G, ρ_L   )

        # Equação de correção de pressão
        u_G, u_L, P_ = pressure_correction_equation_solver( α0_G, α0_L, αx_G, αx_L,
                                                            ux_G, ux_L, Px_, 
                                                            ρ_G, ρ_L,       )

        #= Equação de conservação
        α_G, α_L = void_fraction_equation_solver(   α0_G, α0_L, 
                                                    u_G, u_L, 
                                                    ρ_G, ρ_L, N )

        # Cálculo dos erros
        e_αL = norm(α_L - αx_L) / norm(αx_L)
        e_αG = norm(α_G - αx_G) / norm(αx_G)
        e_uL = norm(u_L - ux_L) / norm(ux_L)
        e_uG = norm(u_G - ux_G) / norm(ux_G)
        e_P = norm(P_ - Px_) / norm(Px_)
        
        # Atualização dos supostos
        αx_G[:] = α_G[:]
        αx_L[:] = α_L[:]
        ux_G[:] = u_G[:]
        ux_L[:] = u_L[:]
        Px_[:] = P_[:]
        =#
        break
    end
    break
    #= Atualização Temporal -------------------------------------------------
    α0_G[:] = α_G[:]
    α0_L[:] = α_L[:]
    u0_G[:] = u_G[:]
    u0_L[:] = u_L[:]
    P0_[:] = P_[:]
    t = t+1
    =#
end
