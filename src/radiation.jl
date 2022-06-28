# This file is part of the ATMOStool.jl module and
# contains functions related to radiation calculations
# including:
# 
# See LICENCE
# ***************************************************

"""
Function to estimate the surface temperature give the LW downwelling
and upwelling radiation fluxes.

julia> Tₛ = T_surf(F_dw, F_up)

or, giving a specific emissivity value of 0.96:

julia> Tₛ = T_surf(F_dw, F_up, εₛ=0.96)

with Tₛ in K.

Optional parameter is εₛ the surface emissivity (default value = 0.985)
"""
function T_surf(Flw_dw, Flw_up; εₛ = 0.985)
    #where σ = 5.670374419E−8 # [W⋅m−2⋅K−4] Stephan-Boltzman constant
    # and ϵₛ the surface emissivity:
    return ((Flw_up - (1-εₛ)*Flw_dw)/εₛ/σ)^0.25
end

"""
Uncertainty function for T_surf due to emissivity

julia> δTₛ = δT_surf(Flw_dw, Flw_up, δε)

with δε the error in percentage e.g. (Δε/ε).
The optional emissivity variable has a default value of εₛ = 0.985

"""
function δT_surf(Flw_dw, Flw_up, δε; εₛ = 0.985)
    faktor0 = let Tₛ= T_surf(Flw_dw, Flw_up, εₛ = εₛ)
        1/4/(Tₛ)^3
    end
    faktor1 = (Flw_up - Flw_dw)/εₛ/σ
    
    return faktor0*faktor1*δε
end
# end of file
