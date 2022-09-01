"""

ACM(lat, doy, t_max, t_min, lai, rad, ca, ce)

This is the Aggregated Canopy Model (ACM),
originally described in Williams et al. (1997),
and coupled to DALEC models in Williams et al., (2005).
Code implementation based on code and materials 
from Fox et al. (2009) "REFLEX" experiment,
and subsequently adapted for Bloom & Williams 2015.

## Arguments
- `lat::Real`: latitude
- `doy::Real`: day-of-year, range: 1-365
- `t_max::Real`: maximum temperature in °C
- `t_min::Real`: minimum temperature in °C
- `lai::Real`: leaf area index
- `rad::Real`: solar irradiance [MJ/m2/day]
- `ca::Real`: atmospheric CO₂ concentration
- `ce::Real`: canopy efficiency
- `FT::Type`: Float type, either Float32 or Float64

## Return
- `gpp::Real`: gross primary production in gC/m^2/day
"""
function ACM(;lat, doy, t_max, t_min, lai, rad, ca, ce, FT)

    # these parameters are directly from the ACM C files,
    # which differ from the original
    # implementation in Williams et al. (1997)
    
    d1 = FT(0.0156935) #day length constant
    θ = FT(4.22273) # 
    k = FT(208.868)
    d2 = FT(0.0453194)
    b2 = FT(0.37836)
    c1 = FT(7.19298)
    a2 = FT(0.011136)
    c2 = FT(2.1001)
    eb1T = FT(0.789798)
    ψ_d = FT(-2)
    H = FT(1)
    
    # compute daily canopy conductance, gc
    gc = (abs(ψ_d)^eb1T) / (b2 * H + FT(0.5) * (t_max - t_min))
    
    # compute p parameter needed for ci
    p = lai * FT(1) * ce * exp(a2 * t_max) / gc
    
    # compute the q parameter needed for ci
    q = θ - k
    
    # compute the internal CO2 concentration, ci
    ci = FT(0.5) * (ca + q - p + sqrt((ca + q - p)^FT(2) - FT(4) * (ca * q - p * θ)))
    #println("ACM ci: ", ci)
    
    # compute canopy-level quantum yield, e0
    e0 = c1 * (lai ^ FT(2)) / (c2 + lai ^ FT(2))
    
    # compute the day length dayl
    dec = FT(-23.4) * cos((FT(360) * (doy + FT(10)) / FT(365)) * FT(π) / FT(180)) * FT(π) / FT(180)
    mult = tan(lat * FT(π) / FT(180)) * tan(dec)
    if mult >= FT(1)
        dayl = FT(24)  
    elseif mult <= -1 
        dayl = FT(0)
    else
      dayl = FT(24) * FT(acos(-mult)) / FT(π)
    end
    
    # compute co2 rate of diffusion to the site of fixation, pd
    pd = gc * (ca - ci)
    
    # compute light limitation pi
    pi = e0 * rad * pd / (e0 * rad + pd)
    
    # compute gpp
    gpp = pi * (d1 * dayl + d2)
    
    return gpp
end;