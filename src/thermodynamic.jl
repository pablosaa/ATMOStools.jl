# This file contains thermodynamic functions
# Part of ATMOStools.jl
# See LICENCE

"""
Virtual Temperature:
input:
* T: Temperature [K]
* qv: water vapour mixing ratio [g/g]
output:
* T_v: Virtual temperature [K]
"""
function T_v(T, qv)
    # Function to compute the virtual Temperature
    ϵ = 0.622
    T_v = T*(1.0 + qv/ϵ)/(1.0 + qv)
    return T_v
end 

"""
Virtual Potential Temperature
Input:
* T : Temperature [K]
* P : Pressure [hPa]
* qv: water vapour mixing ratio [g/g]
Output:
* theta_v: virtual potential temperature [K]
"""
function theta_v(T, P, qv)
    # Function to compute the virtual Potential Temperature:
    κ = 0.286
    T_virtual = T_v(T,qv)
    θ_v =  T_virtual*(100/P)^κ
    return θ_v  
end 

# end of file.
