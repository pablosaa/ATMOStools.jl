# subroutine to implement the Liebe MPM93 model
#=
MPM93 - subroutines adapted by Jeff Haferman (NASA/GSFC 5/97)
  from Liebe's MPM93 model.  His comments are included below.
  I've based this adaptation on Frank Evans' MPM92 extraction.

---------------------------------------------------------------------
       ATMOSPHERIC ATTENUATION AND DELAY RATES UP TO 1000 GHz
       June 1993
       Hans J. Liebe     (303-497-3310)    
       George A. Hufford (       -3457)
       Michael G. Cotton (       -7346)
       Institute for Telecommunication Sciences
       NTIA/ITS.S3 
       325 BROADWAY
       Boulder, CO  80303,  USA

       FAX   :  (303) 497-5993 (ITS), 497-3680 (ITS.S2)
       E-Mail:  HLIEBE@NTIA.ITS.BLDRDOC.GOV

 COMMENTS:
 
   The Millimeter-wave Propagation Model (MPM85) was reported in Ref.
 [1]. Molecular absorption by O2, H2O, and N2 is considered, as well as
 dielectric loss for haze and fog/cloud conditions (Rayleigh absorption
 approximation), and dielectric plus scatter losses (aR**b -
 approximation to Mie's theory) under rain conditions. The complex
 atmospheric refractivity N (or path-specific rates of attenuation A and
 delay B) were continued to be upgraded as discussed in [2] - [7].
 
   Features of the current version, MPM93, are:
 
 - Haze model to predict the water droplet density for 
       U = 80 to 99.95%RH , when a hygroscopic aerosol reference density
       wa(80%RH) and a climatic code ('A, B, C, or D') are provided [2],[3]   
 
 - Improved model for the dielectric properties of liquid water to
       calculate RAYLEIGH absorption and delay by suspended water droplets
       for haze, fog, and cloud conditions [6],[7]
 
 - Rain attenuation model for Laws & Parsons drop-sizes by Olsen et al. 
       [11], and associated dispersive delay, approximated from results 
       reported by Zuffery [12]
 
 - New temperature-dependent linewidth data (b3 to b6) for the water
       vapor lines below 1 THz, and a 5 percent increase in the 
       strength b1 of the 22-GHz and 183-GHz lines [9]
 
 - New set of line mixing coefficients (a5, a6) for dry air, and 
       their improved fit to the extensive 60-GHz lab. data [8],[9]
 
 - Revised water vapor saturation pressure equation [10] 
 
 - Approximation for Zeeman (O2) [4] and Doppler (H2O) line-broadening
       to cover heights up to 100 km.
 
 - New pseudo-line water vapor continuum formulation [9]   
 
 - Detailed treatment of the anisotropic, mesospheric Zeeman effect
   of O2 microwave lines [5]. The ZPM  code [9].
 
 
                                 REFERENCES
 
  [1] H. Liebe, "An updated model for millimeter-wave propagation in
       moist air", Radio Science, vol. 20, no. 5, pp. 1069-1089, 1985.
 
  [2] H. Liebe,"A contribution to modeling atmospheric mm-wave properties",
       FREQUENZ, vol.41, no. 1/2, pp. 31-36, 1987.
 
  [3] H. Liebe and D. Layton, "MM-wave Properties of the Atmosphere:
       Laboratory Studies and Propagation Modeling",
       NTIA Report 87-224, 80p., Oct. 1987 (NTIS Order No. PB88-164215/AF).
       
  [4] H. Liebe,"MPM89 - An atmospheric mm-wave propagation model",
       Int. J. IR & MM Waves, vol.10, no.6, pp. 631-650, June 1989.
 
  [5] G. Hufford and H. Liebe, "MM-Wave Propagation in the Mesosphere",
       NTIA Report 89-249, 67p., Sept. 1989 (NTIS Order No. PB90-119868/AS).
 
  [6] H. Liebe, T. Manabe, and G. Hufford, "Mm-wave attenuation and delay
       rates due to fog/cloud conditions", IEEE Trans. Ant. Prop.,
       vol. 37, no. 12, pp. 1617-1623, Dec. 1989.
 
  [7] H. Liebe, G. Hufford (ice), and T. Manabe, "A model for the complex
       refractivity of water (ice) at frequencies below 1 THz",
       Int. J. IR & MM Waves, vol. 12, no. 7, 659-682, 1991.
   
  [8] H. Liebe, P. Rosenkranz, and G. Hufford, "Atmospheric 60-GHz   
       oxygen spectrum: New laboratory measurements and line parameters", 
       J. Quant. Spectr. Rad. Transf., vol. 48, no. 5/6, pp. 629-643, 1992.
 
  [9] H. Liebe, G. Hufford, and M. Cotton, "Propagation modeling of moist air 
       and suspended water/ice particles at frequencies below 1000 GHz", 
       Proc. AGARD Conf. Paper 3/1-10, Palma De Mallorca, Spain, May 1993.
  
 [10] W. Boegel, "Neue Naeherungsgleichungen fuer den Saettigungsdruck des
       Wasserdampfes, DFVLR Bericht DLR-FB 77-52, 1977.
 
 [11] R.L. Olsen, D.V. Rogers, and D.B. Hodge, "The aRb relation in the
       calculation of rain attenuation",
       IEEE Trans. Ant. Prop., vol. AP-26, no. 2, pp. 318-329, 1978.
 
 [12] C.H. Zuffery, "A study of rain effects on EM waves in the
       1 to 600 GHz range", MS-THesis, Dept. Electrical Eng.,
       University of Colorado, Boulder,  CO 80309, Feb., 1972.
-----------------------------------------------------------------------
=# 


"""
 Computes volume absorption coefficient for an atmospheric
 layer given the meteorological properties. The allowed frequency
 range is from 1 to 1000 GHz.  This routine is hacked from Liebe's
 GAS1 subroutine in his MPM93 model, taking out rain and dispersion
 computations.
    Included is dry air attenuation, oxygen and "psuedo"
 water vapor line-continuum absorption, and Rayleigh cloud droplet
 absorption.

USAGE:
```julia-repl
julia> αₜₒₜ = mpm93(ν, Pₐ, Tc, RH)
julia> αₜₒₜ = mpm93(ν, Pₐ, Tc, RH, qₗ)
```
WHERE:
```julia-repl
* ν  frequency (GHz)
* Pₐ Total pressure (hPa)
* Tc Temperature (C)
* RH Relative Humidity [0…1]
* qₗ  (optional) cloud liquid water content (g m⁻³), default nothing.
```
OUTPUT:
```julia-repl
* αₜₒₜ::Dict  Absorption coefficient (dB km⁻¹) for:
-- :gas   => sum of oxygen, dry air, and water vapour [dB],
-- :O₂    => oxygem [dB],
-- :dry   => dry air [dB],
-- :vapor => water vapour [dB],
-- :cloud => cloud liquid water [dB]
```

(c) 2023 P. Saavedra Garfias

"""
function mpm93(ν, Pₐ, Tc, RH, qₗ)
    
    oxytab = [# frequency  # A1     # A2   # A3   # A4   # A5     # A6   
              50.474239    0.094    9.694  0.890  0.000  0.240    0.790  
          50.987747    0.246    8.694  0.910  0.000  0.220    0.780  
          51.503349    0.608    7.744  0.940  0.000  0.197    0.774  
          52.021412    1.414    6.844  0.970  0.000  0.166    0.764  
          52.542393    3.102    6.004  0.990  0.000  0.136    0.751  
          53.066906    6.410    5.224  1.020  0.000  0.131    0.714  
          53.595749    12.470   4.484  1.050  0.000  0.230    0.584  
          54.130001    22.800   3.814  1.070  0.000  0.335    0.431  
          54.671158    39.180   3.194  1.100  0.000  0.374    0.305  
          55.221367    63.160   2.624  1.130  0.000  0.258    0.339  
          55.783802    95.350   2.119  1.170  0.000  -0.166   0.705  
          56.264774    54.890   0.015  1.730  0.000  0.390    -0.113 
          56.363388    134.400  1.660  1.200  0.000  -0.297   0.753  
          56.968204    176.300  1.260  1.240  0.000  -0.416   0.742  
          57.612484    214.100  0.915  1.280  0.000  -0.613   0.697  
          58.323875    238.600  0.626  1.330  0.000  -0.205   0.051  
          58.446590    145.700  0.084  1.520  0.000  0.748    -0.146 
          59.164207    240.400  0.391  1.390  0.000  -0.722   0.266  
          59.590984    211.200  0.212  1.430  0.000  0.765    -0.090 
          60.306061    212.400  0.212  1.450  0.000  -0.705   0.081  
          60.434776    246.100  0.391  1.360  0.000  0.697    -0.324 
          61.150558    250.400  0.626  1.310  0.000  0.104    -0.067 
          61.800156    229.800  0.915  1.270  0.000  0.570    -0.761 
          62.411217    193.300  1.260  1.230  0.000  0.360    -0.777 
          62.486259    151.700  0.083  1.540  0.000  -0.498   0.097  
          62.997978    150.300  1.665  1.200  0.000  0.239    -0.768 
          63.568520    108.700  2.115  1.170  0.000  0.108    -0.706 
          64.127769    73.350   2.620  1.130  0.000  -0.311   -0.332 
          64.678902    46.350   3.195  1.100  0.000  -0.421   -0.298 
          65.224068    27.480   3.815  1.070  0.000  -0.375   -0.423 
          65.764771    15.300   4.485  1.050  0.000  -0.267   -0.575 
          66.302094    8.009    5.225  1.020  0.000  -0.168   -0.700 
          66.836830    3.946    6.005  0.990  0.000  -0.169   -0.735 
          67.369598    1.832    6.845  0.970  0.000  -0.200   -0.744 
          67.900864    0.801    7.745  0.940  0.000  -0.228   -0.753 
          68.431007    0.330    8.695  0.920  0.000  -0.240   -0.760 
          68.960312    0.128    9.695  0.900  0.000  -0.250   -0.765 
          118.750343   94.500   0.009  1.630  0.000  -0.036   0.009  
          368.498352   6.790    0.049  1.920  0.600  0.000    0.000  
          424.763123   63.800   0.044  1.930  0.600  0.000    0.000  
          487.249359   23.500   0.049  1.920  0.600  0.000    0.000  
          715.393127   9.960    0.145  1.810  0.600  0.000    0.000  
          773.839661   67.100   0.130  1.820  0.600  0.000    0.000  
          834.145325   18.000   0.147  1.810  0.600  0.000    0.000
          ];
          

    watvaptab = [# frequency    # B1       # B2     # B3     # B4    # B5   # B6
                 22.235081      0.01130      2.143    2.811    4.80    0.69   1.00
             67.803963      0.00012      8.735    2.858    4.93    0.69   0.82
             119.995941     0.00008      8.356    2.948    4.78    0.70   0.79
             183.310089     0.24200      0.668    3.050    5.30    0.64   0.85
             321.225647     0.00483      6.181    2.303    4.69    0.67   0.54
             325.152924     0.14990      1.540    2.783    4.85    0.68   0.74
             336.222595     0.00011      9.829    2.693    4.74    0.69   0.61
             380.197357     1.15200      1.048    2.873    5.38    0.54   0.89
             390.134521     0.00046      7.350    2.152    4.81    0.63   0.55
             437.346680     0.00650      5.050    1.845    4.23    0.60   0.48
             439.150818     0.09218      3.596    2.100    4.29    0.63   0.52
             443.018280     0.01976      5.050    1.860    4.23    0.60   0.50
             448.001068     1.03200      1.405    2.632    4.84    0.66   0.67
             470.888947     0.03297      3.599    2.152    4.57    0.66   0.65
             474.689117     0.12620      2.381    2.355    4.65    0.65   0.64
             488.491119     0.02520      2.853    2.602    5.04    0.69   0.72
             503.568542     0.00390      6.733    1.612    3.98    0.61   0.43
             504.482697     0.00130      6.733    1.612    4.01    0.61   0.45
             547.676453     0.97010      0.114    2.600    4.50    0.70   1.00
             552.020935     1.47700      0.114    2.600    4.50    0.70   1.00
             556.935974     48.74000     0.159    3.210    4.11    0.69   1.00
             620.700806     0.50120      2.200    2.438    4.68    0.71   0.68
             645.866150     0.00713      8.580    1.800    4.00    0.60   0.50
             658.005310     0.03022      7.820    3.210    4.14    0.69   1.00
             752.033203     23.96000     0.396    3.060    4.09    0.68   0.84
             841.053955     0.00140      8.180    1.590    5.76    0.33   0.45
             859.962341     0.01472      7.989    3.060    4.09    0.68   0.84
             899.306702     0.00605      7.917    2.985    4.53    0.68   0.90
             902.616150     0.00426      8.432    2.865    5.10    0.70   0.95
             906.207336     0.01876      5.111    2.408    4.70    0.70   0.53
             916.171570     0.83400      1.442    2.670    4.78    0.70   0.78
             923.118408     0.00869      10.220   2.900    5.00    0.70   0.80
             970.315002     0.89720      1.920    2.550    4.94    0.64   0.67
             987.926758     13.21000     0.258    2.985    4.55    0.68   0.90
             1780.000000    2230.00000   0.952    17.620   30.50   2.00   5.00
            ];
                                     
                        
    # *********************************************************************
    # DEFINING VARIABLES AND FUNCTIONS:

    Tₖ = Tc + 273.15    # Temperature in K
    # * Saturated water vapour pressure [hPa]
    eₛ = ATMOStools.PWS(Tₖ)

    # * Pᵥ partial water vapour pressure [hPa]:
    Pᵥ = eₛ*RH
    
    # Dry air Pressure: Atmospheric P minus partial vapour pressure [hPa]:
    𝑃 = Pₐ - Pᵥ
    # Defining the relative inverse temperature [K]:
    𝑇 = 300/Tₖ

    𝐹(νᵢ, δ, γ) = ν/νᵢ*(complex(1., -δ)/complex(νᵢ-ν, -γ) - complex(1., δ)/complex(νᵢ+ν, γ))
    # Defining Attenuation function:
    α(𝑁𝑓, ν) = 0.182ν*imag(𝑁𝑓)
    
    # *********************************************************************
    # STARTING TO COMPUTE ABSORPTION COEFFICIENTS:
    
    #=
    OXYGEN CONTINUUM
    =#

    𝑁𝑓ₒₓ = mapreduce(+, eachrow(oxytab)) do (νᵢ, a₁, a₂, a₃, a₄, a₅, a₆)
        Sᵢ = 1e-6a₁*exp((1 - 𝑇)a₂)*𝑃*𝑇^3
        γᵢ = let γₜ = 1e-3a₃*(𝑃*𝑇^(0.8 - a₄) + 1.1𝑇*Pᵥ)
            √(γₜ^2 + (25*0.6e-4)^2 )
        end

        δᵢ = 1e-3(a₅ + a₆*𝑇)*(𝑃 + Pᵥ)*𝑇^0.8
        # return ΣSᵢFᵢ
        Sᵢ*𝐹(νᵢ, δᵢ, γᵢ)
    end

    # Oxygen line absorption:  AT1=.182*F*AIMAG(ZN)
    αₒₓ = α(𝑁𝑓ₒₓ, ν)

    #=
    DRY AIR CONTINUUM
    =#

    𝑁𝑓ₐᵢᵣ = let γₒ = 0.56e-3(𝑃 + Pᵥ)*𝑇^0.8
        Sₒ = 6.14e-5𝑃*𝑇^2
        𝑁𝑓ₒ = -ν/complex(ν, γₒ)
        Sₙ = 1.40e-12𝑃*𝑃*𝑇^3.5
        𝑁𝑓ₙ = complex(0.0, ν/(1.93e-5ν^1.5 + 1.0) )

        # return (Sₒ𝑁𝑓ₒ + Sₙ𝑁𝑓ₙ)
        Sₒ*𝑁𝑓ₒ + Sₙ*𝑁𝑓ₙ
    end

    # NONRESONAT DRY AIR ABSORPTION: 
    αₐᵢᵣ =  α(𝑁𝑓ₐᵢᵣ, ν)

    #=
    WATER VAPOUR
    =#
    𝑁𝑓ᵥ = mapreduce(+, eachrow(watvaptab)) do (νₕ₂ₒ, b₁, b₂, b₃, b₄, b₅, b₆)
        S = b₁*exp((1 - 𝑇)b₂)*Pᵥ*𝑇^3.5
        # Doppler approximation
        
        γd2 = 1e-12/𝑇*(1.46*νₕ₂ₒ)^2
        γₕ = 1e-3b₃*(𝑃*𝑇^b₅ + b₄*Pᵥ*𝑇^b₆)
        γₕ = 0.535γₕ + √(0.217γₕ^2 + γd2 )
        δₕ = 0.0
        𝑍𝑓 = 𝐹(νₕ₂ₒ, δₕ, γₕ)
        # return ΣS*𝑍𝑓
        S*𝑍𝑓
    end

    #= WATER VAPOR LINE ABSORPTION 
    SEE LIEBE'S COMMENT REGARDING "PSUEDO-LINE WATER VAPOR CONTINUUM" - JLH
    =#
    αᵥ = α(𝑁𝑓ᵥ, ν)

    #= ********************************************************************* 
     LIQUID WATER PERMITTIVITY [8]
     Use exponential form for gamma for T<0 extrapolation (a la Frank Evans)
    =#

    if !isnothing(qₗ) && qₗ ≠ 0.0
        fD = 20.1*exp(7.88*(1-𝑇)) #20.20-146.4*(𝑇 - 1) + 316*(𝑇 - 1)^2 # 
        fS=39.8fD
        ϵ = 103.3*(𝑇-1)+77.66
        ϵ_∞ = 0.0671ϵ

        Eopt = 3.52
        # Complex Permittivity of water (double-Debye model)
        ZEp = ϵ - ν*((ϵ - ϵ_∞)/complex(ν, fD) + (ϵ_∞ - Eopt)/complex(ν, fS))
    
        #=C
        C ICE PERMITTIVITY [8]
        (PSG: not implemented yet) =#

        # SUSPENDED PARTICLE RAYLEIGH APPROXIMATION [6]
        ZNw = 1.5qₗ*((ZEp - 1)/(ZEp + 2) - 1 + 3/(ϵ + 2))

        # SUSPENDED WATER DROPLET EXTINCTION 
        αₗ = α(ZNw, ν)
    else
        αₗ = 0.0
    end

    αₜₒₜₐₗ = (O₂ = αₒₓ,
              dry = αₐᵢᵣ,
              vapor = αᵥ,
              gas = αₒₓ + αₐᵢᵣ + αᵥ,
              cloud = αₗ)

    return αₜₒₜₐₗ

end
# -- OR --
function mpm93(ν::AbstractFloat, Pₐ::Matrix, Tc::Matrix, RH::Matrix, qₗ::Matrix)

    αₐₜₜ = mpm93.(ν, Pₐ, Tc, RH, qₗ)

    αₜₒₜₐₗ = Dict(k=>getproperty.(αₐₜₜ, k) for k ∈ keys(αₐₜₜ[1]))

    return αₜₒₜₐₗ
end

# end of file

