# Equação de momento

## Equação diferencial:

$\left(\alpha_k \rho_k \right) \dfrac{\partial u_k}{\partial t} + \left(\alpha_k \rho_k u_k \right) \dfrac{\partial u_k}{\partial x} = -\alpha_k \dfrac{\partial P}{\partial x} - \left( p_k - p^i_k \right) \dfrac{\partial  \alpha_k}{\partial x} + \alpha_k \rho_kg \sin\theta$

Discretizando ( em $e$ ):

$\left(\alpha_k \rho_k \right)^n_e \dfrac{\left( u_k \right)^{n+1}_e - \left( u_k \right)^{n}_e}{\Delta t} + \left(\alpha_k \rho_k u_k\right)^n_e \dfrac{\left( u_k \right)^{n+1}_E - \left( u_k \right)^{n+1}_P}{\Delta x}= -\left( \alpha_k \right)^{n}_e \dfrac{P^{n+1}_E - P^{n+1}_P}{\Delta x} + \left( \alpha_k \rho_k \right)^n_e g \sin\theta + \Delta P^n_e\dfrac{\alpha^n_{kE} - \alpha^n_{kP}}{\Delta x}$

A equação pode ser reescrita como:

$a_e \left( u_{ke} \right) = a_w \left( u_{kw} \right) + a_{ee} \left( u_{kee} \right) + A_e \left( P_P - P_E \right) + b$

Utilizado o esquema Upwind de primeira ordem:

$a_e = \left[ \dfrac{\left( \alpha_k \rho_k \right)^n_e}{\Delta t} + \left( \dfrac{\omega_{kE}}{2} + \dfrac{\omega_{kP}}{2} \right) \dfrac{\left( \alpha_k \rho_k u_k \right)^n_e}{\Delta x} \right]$

$a_w = \left( \dfrac{1}{2} + \dfrac{\omega_P}{2} \right) \dfrac{\left( \alpha_k \rho_k u_k \right)^n_e}{\Delta x}$

$a_{ee} = -\left( \dfrac{1}{2} - \dfrac{\omega_E}{2} \right) \dfrac{\left( \alpha_k \rho_k u_k \right)^n_e}{\Delta x}$

$A_e = \dfrac{\alpha^n_{ke}}{\Delta x}$

$b = \dfrac{\left( \alpha_k \rho_k u_k \right)^n_{e}}{\Delta t} + \left( \alpha_k \rho_k \right)^n_e g \sin\theta + \Delta P^n_e\dfrac{\alpha^n_{kE} - \alpha^n_{kP}}{\Delta x}$

Onde $\Delta P^n_e$ é dado pelo modelo Cathare:

$\Delta P^n_e = \gamma \dfrac{\left( \alpha_g \alpha_l \rho_g \rho_l \right)^n_e}{\left( \alpha_g \rho_l \right)^n_k + \left( \alpha_l \rho_g \right)^n_e}\, , \quad \gamma = 1.2$

# Equação de correção de pressão

## Forma discretizada ( em $P$ ):

$\dfrac{1}{\rho_{g,ref}} \left[ \left( \alpha_g \rho_g \right)^\kappa_e u^{\kappa+1}_{ge} - \left( \alpha_g \rho_g \right)^\kappa_w u^{\kappa+1}_{gw} \right] + \dfrac{1}{\rho_{l,ref}} \left[ \left( \alpha_l \rho_l \right)^\kappa_e u^{\kappa+1}_{le} - \left( \alpha_l \rho_l \right)^\kappa_w u^{\kappa+1}_{lw} \right] + \dfrac{\Delta x}{\Delta t} \left\{ \dfrac{\left( \alpha^\kappa_g \rho^{\kappa+1}_g \right)_P - \left( \alpha_g \rho_g \right)^n_P}{\rho_{g,ref}} + \dfrac{\left( \alpha^\kappa_l \rho^{\kappa+1}_l \right)_P - \left( \alpha^n_l \rho^n_l \right)_P}{\rho_{l,ref}} + \right\} = 0$

considerando:

$u^{\kappa+1}_{ke} = u^\kappa_{ke} + \delta u_{ke}\, , \quad \delta u_{ke} = \left( \dfrac{A_e}{a_e} \right)_k \left( \delta P_P - \delta P_E \right)$

$A_{ke} = \dfrac{\alpha^n_{ke}}{\Delta x}\, , \quad a_{ke} = \left[ \dfrac{\left( \alpha_k \rho_k \right)^n_e}{\Delta t} + \left( \dfrac{\omega_{kE}}{2} + \dfrac{\omega_{kP}}{2} \right) \dfrac{\left( \alpha_k \rho_k u_k \right)^n_e}{\Delta x} \right]$

$u^{\kappa+1}_{kw} = u^\kappa_{kw} + \delta u_{kw}\, , \quad \delta u_{kw} = \left( \dfrac{A_w}{a_w} \right)_k \left( \delta P_W - \delta P_P \right)$

$A_{kw} = \dfrac{\alpha^n_{kw}}{\Delta x}\, , \quad a_{kw} = \left[ \dfrac{\left( \alpha_k \rho_k \right)^n_w}{\Delta t} + \left( \dfrac{\omega_{kP}}{2} + \dfrac{\omega_{kW}}{2} \right) \dfrac{\left( \alpha_k \rho_k u_k \right)^n_w}{\Delta x} \right]$

$\left(\rho_k \right)^{\kappa+1}_P \approxeq \left( \rho_k \right)^\kappa_P + \left( \dfrac{1}{c^2_k} \right) \delta P$

A equação pode ser reescrita como:

$\left[ \dfrac{\left( \alpha_g \rho_g \right)^\kappa_e}{\rho_{g,ref}} \left( \dfrac{A_e}{a_e} \right)_g +  \dfrac{\left( \alpha_g \rho_g \right)^\kappa_w}{\rho_{g,ref}} \left( \dfrac{A_w}{a_w} \right)_g + \dfrac{\Delta x}{\Delta t} \dfrac{\alpha^\kappa_{gP}}{\rho_{g,ref}} \dfrac{1}{c^2_g} + \dfrac{\left(\alpha_l \rho_l\right)^\kappa_e}{\rho_{l,ref}} \left( \dfrac{A_e}{a_e} \right)_l + \dfrac{\left(\alpha_l \rho_l\right)^\kappa_w}{\rho_{l,ref}} \left( \dfrac{A_w}{a_w} \right)_l + \dfrac{\Delta x}{\Delta t} \dfrac{\alpha^\kappa_{lP}}{\rho_{l,ref}} \dfrac{1}{c^2_l} \right] \delta P_P$
$+ \left[ - \dfrac{\left(\alpha_g \rho_g\right)^\kappa_e}{\rho_{g,ref}} \left( \dfrac{A_e}{a_e} \right)_g - \dfrac{(\alpha_l \rho_l)^\kappa_e}{\rho_{l,ref}} \left( \dfrac{A_e}{a_e} \right)_l \right] \delta P_E$
$+ \left[ - \dfrac{\left(\alpha_g \rho_g\right)^\kappa_w}{\rho_{g,ref}} \left( \dfrac{A_w}{a_w} \right)_g - \dfrac{(\alpha_l \rho_l)^\kappa_w}{\rho_{l,ref}} \left( \dfrac{A_w}{a_w} \right)_l \right] \delta P_W$

$= \dfrac{\left(\alpha_g \rho_g u_g \right)^\kappa_w - \left(\alpha_g \rho_g u_g \right)^\kappa_e}{\rho_{g,ref}} + \dfrac{\left(\alpha_l \rho_l u_l \right)^\kappa_w - \left(\alpha_l \rho_l u_l \right)^\kappa_e}{\rho_{l,ref}} + \dfrac{\Delta x}{\Delta t} \left\{ \dfrac{\left(\alpha_g \rho_g \right)^n_P - \left(\alpha_g \rho_g \right)^\kappa_P}{\rho_{g,ref}} + \dfrac{\left(\alpha_l \rho_l \right)^n_P - \left(\alpha_l \rho_l \right)^\kappa_P}{\rho_{l,ref}}\right\}$

# Equação de conservação da massa

## Equação diferencial:

$\dfrac{\partial}{\partial t} \left( \alpha_k \rho_k \right) + \dfrac{\partial}{\partial x} \left( \alpha_k \rho_k u_k \right) = 0$

Discretizando ( em $P$ ):

$\left( \alpha_k \rho_k \right)_P^{\kappa+1} = \left( \alpha_k \rho_k \right)_P^n + \dfrac{\Delta t}{\Delta x} \left[ \left( \alpha_k \rho_k u_k \right)_e^{\kappa+1} - \left( \alpha_k \rho_k u_k \right)_w^{\kappa+1} \right]$

Onde:

$\left( \alpha \right)^{\kappa +1}_e = \left( \dfrac{1}{2} + \dfrac{\omega_e}{2} \right) \alpha^{\kappa +1}_P + \left( \dfrac{1}{2} - \dfrac{\omega_e}{2} \right) \alpha^{\kappa +1}_E$

$\left( \alpha \right)^{\kappa +1}_w = \left( \dfrac{1}{2} + \dfrac{\omega_w}{2} \right) \alpha^{\kappa +1}_W + \left( \dfrac{1}{2} - \dfrac{\omega_w}{2} \right) \alpha^{\kappa +1}_P$

Rearranjando:

$\left[ \dfrac{\left( \rho_k \right)^{\kappa +1}_P}{\Delta t} - \left( \dfrac{1}{2} + \dfrac{\omega_e}{2} \right) \dfrac{\left( \rho_k u_k \right)^{\kappa +1}_e}{\Delta x} + \left( \dfrac{1}{2} - \dfrac{\omega_w}{2} \right) \dfrac{\left( \rho_k u_k \right)^{\kappa +1}_w}{\Delta x} \right] \left( \alpha_k \right)^{\kappa +1}_P$
$+ \left[ - \left( \dfrac{1}{2} - \dfrac{\omega_e}{2} \right) \dfrac{\left( \rho_k u_k \right)^{\kappa +1}_e}{\Delta x} \right] \left( \alpha_k \right)^{\kappa +1}_E$
$+ \left[ \left( \dfrac{1}{2} + \dfrac{\omega_w}{2} \right) \dfrac{\left( \rho_k u_k \right)^{\kappa +1}_w}{\Delta x} \right] \left( \alpha_k \right)^{\kappa +1}_W$
$= \dfrac{\left( \alpha_k \rho_k \right)^n_P}{\Delta t}$
