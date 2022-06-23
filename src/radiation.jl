# This file is part of the ATMOStool.jl module and
# contains functions related to radiation calculations
# including:
# 
# See LICENCE
# ***************************************************

"""
Function to estimate the surface temperature give the LW downwelling
and upwelling radiation fluxes.

> Tₛ = T_surf(F_dw, F_up)

with Tₛ in K.

Optional parameter is ϵₛ the surface emissivity (default value = 0.98)
"""
function T_surf(Flw_dw, Flw_up; ϵₛ = 0.98)
    #where σ = 5.670374419E−8 # [W⋅m−2⋅K−4] Stephan-Boltzman constant
    # and ϵₛ the surface emissivity:
    return ((Flw_up - (1-ϵₛ)*Flw_dw)/ϵₛ/σ)^0.25
end


# end of file
