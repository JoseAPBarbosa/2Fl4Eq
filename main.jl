using Printf, Plots, CSV, DataFrames

include("./initialization.jl")
include("./equation_solvers.jl")


cont = 0
while (e_uL < Tol_u && e_uG < Tol_u && e_αL < Tol_α && e_αG < Tol_α && e_P < Tol_P)  ||  (cont < 1)
    global α_G, α_L, ρ_G, ρ_L, u_G, u_L, P_
    global αø_G, αø_L, ρø_G, ρø_L, uø_G, uø_L, Pø_
    global cG, cL, e_uL, e_uG, e_αL, e_αG, e_P, cont

    # Equação do momento
    uø_G, uø_L = momentum_conservation_equation_solver( αø_G, ρø_G, uø_G, 
                                                        αø_L, ρø_L, uø_L, 
                                                        α_G, ρ_G, u_G, 
                                                        α_L, ρ_L, u_L, 
                                                        Pø_ )

    #= Equação de correção de pressão
    ρØ_G, uØ_G, ρØ_L, uØ_L, PØ_ = pressure_correction_equation_solver(αØ_G, ρØ_G, uØ_G, αØ_L, ρØ_L, uØ_L, PØ_, α_G, ρ_G, α_L, ρ_L)

    # Equação de conservação
    αØ_G, αØ_L = void_fraction_equation_solver(αØ_G, ρØ_G, uØ_G, αØ_L, ρØ_L, uØ_L, α_G, ρ_G, α_L, ρ_L)

    # Atualização dos supostos ----------------------------------------------------------------------------------------------------------------------
    e_uL = norm(u_L - uØ_L) / norm(uØ_L)
    e_uG = norm(u_G - uØ_G) / norm(uØ_G)
    e_αL = norm(α_L - αØ_L) / norm(αØ_L)
    e_αG = norm(α_G - αØ_G) / norm(αØ_G)
    e_P = norm(P_ - PØ_) / norm(PØ_)
=#
    cont = cont + 1

    #display(cont)

end

display("ok")
