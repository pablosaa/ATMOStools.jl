# This file is part of the ATMOStool.jl module and
# contains atmospheric thermodynamic functions
# 
# See LICENCE

# ********************************************************************
# Defining global Atmospheric constatnts:
#
# Boltzman gas constant R (dry air)
const Rd = 287.058   # [J/kg/K]

# Rv gas constant water vapour:
const Rv = 461.5  # [J/kg/K]

# Specific Heat capacity (isobaric)
const c_p = 1e3  # [J/kg/K]

# Ratio R/c_p
const κ = Rd/c_p  # ~ 0.286

# Reference Pressure P0 [hPa]
const P_0 = 1000.0

# Gravity acceleration g [m/s²]
const gravity_0 = 9.81

const ϵ = Rd/Rv   # ~ 0.622

# **********************************************************************
# Function potential temperature
"""
    θ = θ(T, P)

INPUT:
* T -> Temperature [K]
* P -> Pressure [hPa]

OUTPUT:
* θ -> Potential temperature [K]
"""
function θ(T, P)
    return T*(P_0/P)^κ
end
function θ(T::Matrix, P::Matrix)::Matrix
    return θ.(T, P)
end
# ----/
# ********************************************************************
# Function Virtual Temperature
"""
! ____________________________________________________________________
! Subroutine to calculate the Virtual Temperature given T, MIXR
!
! --------------------------------------------------------------------
! Function convert Temperature to Virtual Temp.
! -> T    : Temperature [K]
! -> MIXR : Vapour Mixing ratio [kg/kg]
! <- Tv   : Virtual Temperature [K]
! ---
! (c) 2020, Pablo Saavedra G.
! Geophysical Institute, University of Bergen
! See LICENSE
! ---
"""
function VirtualTemperature(T, MIXR)
    #real(kind=8), intent(in) :: T, MIXR
    #real(kind=8) :: TV
    #! ** local variables
    TV = T*(1.0 + MIXR/ϵ )/(1.0 + MIXR)

    return TV
end #function VirtualTemperature
function VirtualTemperature(T::Matrix, MIXR::Matrix)
    return VirtualTemperature.(T, MIXR)
end
# ----/

# **********************************************************************
# Alternative Function for Virtual Temperature
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
function T_v(T::Matrix, qv::Matrix)::Matrix
    return T_v.(T, qv)
end
# ----/

# *********************************************************************
# Function Virtual Potential Temperature
"""
Virtual Potential Temperature

θ_v = Theta_virtual(T, P, qv)

Input:
* T : Temperature [K]
* P : Pressure [hPa]
* qv: water vapour mixing ratio [g/g]
Output:
* θ_v: virtual potential temperature [K]
"""
function Theta_virtual(T, P, qv)
    # Function to compute the virtual Potential Temperature:
    
    T_virtual = VirtualTemperature(T, qv)
    return θ(T_virtual, P)  
end
function Theta_virtual(T::Matrix, P::Matrix, qv::Matrix)::Matrix
    return Theta_virtual.(T, P, qv)
end
# ----/

"""
! -------------------------------------------------------------------
! Function to convert Specific [kg/kg] to mixing rations [kg/kg].
! Input, Q_x  : Specific content in [kg/kg]
! Output, mixr : Mixing ratio in [kg/kg]
! ---
"""
function qx_to_mixr(Q_x)
  
  return Q_x/(1.0 - Q_x)
end #function qx_to_mixr
function qx_to_mixr(Q_x::Matrix)::Matrix
    return qx_to_mixr.(Q_x)
end
# ----/

"""
! -------------------------------------------------------------
! Function to calculate Partial pressure water vapour
-> T  : Temperature [K]
-> P  : Pressure [hPa]
<- Es : Partial pressure [hPa]
! ---
"""
function PWS(T)

    COEFF = 2.16679 # [g K J^-1]
    Tc = 647.096  # critical temperature [K]
    Pc = 220640   # critical pressure [hPa]
    B  = 0.6219907 # constant for air [kg/kg]
    CC = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
    EE = (1.0, 1.5, 3.0, 3.5, 4.0, 7.5)

    etha = 1.0 - T/Tc
    A = 0.0

    map(1:6) do i
        A = A + CC[i]*(etha^EE[i])
    end
    
    E_s = Pc*exp(A*Tc/T)
    return E_s
end
# ----/

"""
! ---------------------------------------------------------------
! ELEMENTAL FUNCTION Convert mixing ratio to Relative Humidity.
! HUMIDITY CONVERSION FORMULAS
! Calculation formulas for humidity (B210973EN-F)
! By (c) VAISALA 2003
! https://www.hatchability.com/Vaisala.pdf
!
! ---
! (c) 2019 Pablo Saavedra G.
! Geophysical Institute, University of Bergen
! See LICENSE
!
! ---
! -> MIXR: Vapour mixing ratio [kg/kg]
! -> P   : Pressure [hPa]
! -> T   : Temperature [K]
! <- RH  : Relative Humidity [%]
"""
function mixr_to_rh(MIXR, P, T)
    B  = 0.6219907 # constant for air [kg/kg]

    #RH = 100*PW/PWS
    RH = 100*MIXR*P/(MIXR + B)/PWS(T)
    return RH
end ##function mixr_to_rh
function mixr_to_rh(MIXR::Matrix, P::Matrix, T::Matrix)::Matrix
    return mixr_to_rh.(MIXR, P, T)
end
# ----/

"""
! _________________________________________________________________
! ELEMENTAL FuNCTION Specific Humidity to Relative Humidity
! -> Qv : Specific Humidity [kg/kg]
! -> P  : Pressure [hPa]
! -> T  : Temperature [K]
! <- RH : Relative Humidity [%]
! ---
"""
function qv_to_rh(QV, P, T)
  
  # Converting Specific Humidity to Mixing Ratio:
  MIXR = qx_to_mixr(QV)
  
  RH = mixr_to_rh(MIXR, P, T)

  return RH
end #function qv_to_rh
# ----/

"""
! --------------------------------------------------------------------
! Function to convert kg/kg to kg/m^3
! -> Qx  : Specific quantity e.g. humidity [kg/kg]
! -> T   : Temperature [K]
! -> P   : Pressure [hPa]
! -> MIXR: vapour mixing ration [kg/kg]
! <- RHOx: Specific quantity [kg/m^3]
! ---
"""
function MassRatio2MassVolume(Q_x, T, P ,MIXR)
  #real(kind=8), intent(in) :: Q_x, T, P, MIXR
  #real(kind=8) :: RHO_x
  # Local variables:
  # real(kind=8) :: Tv, RHO_air
#  const Rd = 287  # [J/kg/K]
#  const Rv = 461.5  # [J/kg/K]
  
  Tv = VirtualTemperature(T, MIXR)
  RHO_air = (1.0E2*P)/Tv/Rd  # [kg/m^3]
  RHO_x   = Q_x*RHO_air      # [kg/m^3]
  
  return RHO_x
end #function MassRatio2MassVolume
# ----/

"""
! ----------------------------------------------------------------------
! Function to convert dew point to relative humidity
! -> T  : Temperature [K]
! -> Td : Dew point []
! <- RH : Relative Humidity [%]
! ---
"""
function dewpoint_to_rh(T, Td)

    E_0 = 6.11    # [hPa]
    T_0 = 273.15  # [K]
    L_Rv = 5423.0 # [K]

    # according to an approximation of the Clausius-Clapeyron Equation (CCE):
    CCE(Tx) = E_0*exp(L_Rv*(1/T_0 - 1/Tx))
    
    RH = 100.0*CCE(Td)/CCE(T)
    return RH
end
function dewpoint_to_rh(T::Matrix, Td::Matrix)::Matrix
    return dewpoint_to_rh.(T, Td)
end
# ----/

"""
! -----------------------------------------------------------------------
! Function to convert from relative humidity to mixing ratio
! -> RH  : relative humidity [%]
! -> T   : temperature [K]
! -> P   : pressure [hPa]
! <- mixr: mixing ratio [g/g]
"""
function rh_to_mixr(RH, T, P)

    B  = 0.6219907 # constant for air [kg/kg]

    Es = PWS(T)
    mixr = RH*B*Es/(100.0*P - RH*Es)
    return mixr
end
function rh_to_mixr(RH::Matrix, T::Matrix, P::Matrix)::Matrix
    return rh_to_mixr.(RH, T, P)
end
# ----/

# ***************************************************************************
# end of file.
