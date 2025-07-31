using Printf, Plots, CSV, DataFrames

include("./initialization.jl")
include("./equation_solvers.jl")


cont = 0
while (e_uL < Tol_u && e_uG < Tol_u && e_αL < Tol_α && e_αG < Tol_α && e_P < Tol_P)  ||  (cont < 1)
    global α_G, α_L, u_G, u_L, P_
    global αø_G, αø_L, uø_G, uø_L, Pø_
    global ρ_G, ρ_L, e_uL, e_uG, e_αL, e_αG, e_P, cont

    # Equação do momento
    uø_G, uø_L = momentum_conservation_equation_solver( αø_G, αø_L, α_G, α_L,
                                                        uø_G, uø_L, u_G, u_L,
                                                        Pø_, ρ_G, ρ_L   )

    # Equação de correção de pressão
    uø_G, uø_L, Pø_ = pressure_correction_equation_solver(  αø_G, αø_L, α_G, α_L,
                                                            uø_G, uø_L, Pø_,
                                                            ρ_G, ρ_L    )

    # Equação de conservação
    αø_G, αø_L = void_fraction_equation_solver( α_G, α_G,
                                                uø_L, uø_L, 
                                                ρ_G, ρ_L    )

    # Atualização dos supostos ----------------------------------------------------------------------------------------------------------------------
    α_G = αø_G
    α_L = αø_L
    u_G = uø_G
    u_L = uø_L
    P_ = Pø_    
    
    e_αL = norm(α_L - αø_L) / norm(αø_L)
    e_αG = norm(α_G - αø_G) / norm(αø_G)
    e_uL = norm(u_L - uø_L) / norm(uø_L)
    e_uG = norm(u_G - uø_G) / norm(uø_G)
    e_P = norm(P_ - Pø_) / norm(Pø_)

    cont = cont + 1

    display(cont)

end

display("ok")
