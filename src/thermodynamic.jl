# This file contains thermodynamic functions
# Part of ATMOStools.jl
# See LICENCE

"""
! _______________________________________________________________________
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
  #! Rd gas constant dry air, Rv gas constant water vapour:
  const Rd = 287  ! [J/kg/K]
  const Rv = 461.5  ! [J/kg/K]
  const ϵ = Rd/Rv           ! ~ 0.622

  TV = T*(1.0 + MIXR/ϵ )/(1.0 + MIXR)
  return TV
end #function VirtualTemperature
# ----/


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

"""
! -------------------------------------------------------------------
! Function to convert Specific [kg/kg] to mixing rations [kg/kg].
! Input, Q_x  : Specific content in [kg/kg]
! Output, mixr : Mixing ratio in [kg/kg]
! ---
"""
function qx_to_mixr(Q_x)
  
  MIXR = Q_x/(1.0 - Q_x)
  return MIXR
end #function qx_to_mixr
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
  real :: PWS, etha, A
  const COEFF = 2.16679 # [g K J^-1]
  const Tc = 647.096  # critical temperature [K]
  const Pc = 220640   # critical pressure [hPa]
  const B  = 0.6219907 # constant for air [kg/kg]
  const CC = [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502]
  const EE = [1.0, 1.5, 3.0, 3.5, 4.0, 7.5]

  etha = 1.0 - T/Tc
  A = 0.0
  #do i=1, 6
  #   A = A + CC(i)*(etha**EE(i))
  #end do
  map(1:6) do i
       A = A + CC[i]*(etha**EE[i])
  end
  
  PWS = Pc*exp(A*Tc/T)
  #RH = 1E-3*MIXR*T/PWS/COEFF
  #RH = 100*PW/PWS
  RH = 100*MIXR*P/(MIXR + B)/PWS
  return RH
end ##function mixr_to_rh
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
  #implicit none
  #real(kind=8), intent(in) :: QV, P, T
  #real(kind=8) :: RH

  # Local variables:
  #real(kind=8) :: MIXR
  
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
  const Rd = 287  # [J/kg/K]
  const Rv = 461.5  # [J/kg/K]
  
  Tv = VirtualTemperature(T, MIXR)
  RHO_air = (1.0E2*P)/Tv/Rd  # [kg/m^3]
  RHO_x   = Q_x*RHO_air      # [kg/m^3]
  
  return RHO_x
end #function MassRatio2MassVolume
# ----/

# end of file.
