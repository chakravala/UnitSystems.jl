
#   This file is part of UnitSystems.jl. It is licensed under the MIT license
#   UnitSystems Copyright (C) 2020 Michael Reed

@doc """
    T/kelvin == (9/5)*T

Converts temperature `T` from Kelvin to degrees Rankine (°R).
""" kelvin

@doc """
    T/rankine == (5/9)*T

Converts temperature `T` from degrees Rankine to Kelvin (K).
""" rankine

"""
    moles(N::Real,U::UnitSystem=Metric) = N/avogadro(U)

Converts the number of molecules `N` to number of moles (mol).
"""
@pure moles(N::Real,U::UnitSystem=Metric) = N/avogadro(U)

"""
    molecules(n::Real,U::UnitSystem=Metric) = n*avogadro(U)

Converts the number of moles `n` to number of molecules (dimensionless).
"""
@pure molecules(n::Real,U::UnitSystem=Metric) = n*avogadro(U)

# thermodynamics

"""
$(convertext(:temperature,"mass(U,S)*speed(U,S)^2/boltzmann(U,S)"))

Measurement scale for thermodynamic energy or `temperature` (K), unit conversion factor.

```Julia
julia> temperature(Metric,SI2019) # K⋅K⁻¹
$(temperature(Metric,SI2019))

julia> temperature(English,SI2019) # K⋅°R⁻¹
$(temperature(English,SI2019))

julia> temperature(English,Metric) # K⋅°R⁻¹
$(temperature(English,Metric))

julia> temperature(PlanckGauss,Metric) # K⋅TP⁻¹
$(temperature(PlanckGauss,Metric))
```
"""
@pure temperature(U::UnitSystem,S::UnitSystem) = unit((boltzmann(U)*electronmass(S)*lightspeed(S)^2)/(boltzmann(S)*electronmass(U)*lightspeed(U)^2))

"""
$(convertext(:entropy,"energy(U,S)/temperature(U,S)"))

Heat capacity or `energy` per `temperature` or `entropy` (J⋅K⁻¹), unit conversion factor.

```Julia
julia> entropy(Metric,SI2019) # K⋅K⁻¹
$(entropy(Metric,SI2019))

julia> entropy(CGS,Metric) # J⋅erg⁻¹
$(entropy(CGS,Metric))

julia> entropy(English,SI2019) # J⋅°R⋅K⁻¹⋅ft⁻¹⋅lb⁻¹
$(entropy(English,SI2019))

julia> entropy(EnglishUS,English) # ftUS²⋅°R⋅°ft⁻²⋅°R⁻¹
$(entropy(EnglishUS,English))
```
"""
@pure entropy(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)/temperature(U,S))

"""
$(convertext(:specificentropy,"specificenergy(U,S)/temperature(U,S)"))

Specific heat capacity or `specificentropy` (J⋅K⁻¹⋅kg⁻¹), unit conversion factor.

```Julia
julia> specificentropy(Metric,SI2019) # m²⋅K⋅K⁻¹⋅cm⁻²
$(specificentropy(Metric,SI2019))

julia> specificentropy(CGS,Metric) # m²⋅cm⁻²
$(specificentropy(CGS,Metric))

julia> specificentropy(English,Metric) # m²⋅°R⋅K⁻¹⋅ft⁻²
$(specificentropy(English,Metric))

julia> specificentropy(EnglishUS,English) # ft²⋅°R⋅ftUS⁻²⋅°R⁻¹
$(specificentropy(EnglishUS,English))
```
"""
@pure specificentropy(U::UnitSystem,S::UnitSystem) = unit(specificenergy(U,S)/temperature(U,S))

"""
$(convertext(:volumeheatcapacity,"entropy(U,S)/volume(U,S)"))

The `entropy` per `volume` or `volumeheatcapacity` (J⋅K⁻¹⋅m⁻³), unit conversion factor.

```Julia
julia> volumeheatcapacity(Metric,SI2019) # K⋅K⁻¹
$(volumeheatcapacity(Metric,SI2019))

julia> volumeheatcapacity(CGS,Metric) # J⋅cm³⋅erg⁻¹⋅m⁻³
$(volumeheatcapacity(CGS,Metric))

julia> volumeheatcapacity(English,SI2019) # J⋅ft²⋅°R⋅K⁻¹⋅lb⁻¹⋅m⁻³
$(volumeheatcapacity(English,SI2019))

julia> volumeheatcapacity(EnglishUS,English) # ftUS⁵°R⋅°ft⁻⁵⋅°R⁻¹
$(volumeheatcapacity(EnglishUS,English))
```
"""
@pure volumeheatcapacity(U::UnitSystem,S::UnitSystem) = unit(entropy(U,S)/volume(U,S))

"""
$(convertext(:thermalconductivity,"force(U,S)/time(U,S)/temperature(U,S)"))

Heat conductivity or `thermalconductivity` (W⋅m⁻¹⋅K⁻¹), unit conversion factor.

```Julia
julia> thermalconductivity(Metric,SI2019) # K⋅K⁻¹
$(thermalconductivity(Metric,SI2019))

julia> thermalconductivity(CGS,Metric) # N⋅dyn⁻¹
$(thermalconductivity(CGS,Metric))

julia> thermalconductivity(English,Metric) # N⋅°R⋅K⁻¹⋅ft⁻¹⋅lb⁻¹
$(thermalconductivity(English,Metric))
```
"""
@pure thermalconductivity(U::UnitSystem,S::UnitSystem) = unit(force(U,S)/time(U,S)/temperature(U,S))

"""
$(convertext(:thermalconductance,"thermalconductivity(U,S)*length(U,S)"))

Reciprocal of `thermalresistance` (W⋅K⁻¹), unit conversion factor.

```Julia
julia> thermalconductance(Metric,SI2019) # K⋅K⁻¹
$(thermalconductance(Metric,SI2019))

julia> thermalconductance(CGS,Metric) # W⋅s⋅erg⁻¹
$(thermalconductance(CGS,Metric))

julia> thermalconductance(English,Metric) # J⋅°R⋅K⁻¹⋅ft⁻¹⋅lb⁻¹
$(thermalconductance(English,Metric))
```
"""
@pure thermalconductance(U::UnitSystem,S::UnitSystem) = unit(thermalconductivity(U,S)*length(U,S))

"""
$(convertext(:thermalresistance,"1/thermalconductivity(U,S)/length(U,S)"))

Resistance to heat flow or `thermalresistance` (K⋅W⁻¹), unit conversion factor.

```Julia
julia> thermalresistance(Metric,SI2019) # K⋅K⁻¹
$(thermalresistance(Metric,SI2019))

julia> thermalresistance(CGS,Metric) # erg⋅s⁻¹⋅W⁻¹
$(thermalresistance(CGS,Metric))

julia> thermalresistance(English,Metric) # ft⋅lb⋅K⋅°R⁻¹⋅J⁻¹
$(thermalresistance(English,Metric))
```
"""
@pure thermalresistance(U::UnitSystem,S::UnitSystem) = thermalconductance(S,U)

"""
$(convertext(:thermalexpansion,"1/temperature(U,S)"))

Measurement scale for coefficient of `thermalexpansion` (K⁻¹), unit conversion factor.

```Julia
julia> thermalexpansion(Metric,SI2019) # K⋅K⁻¹
$(thermalexpansion(Metric,SI2019))

julia> thermalexpansion(English,SI2019) # °R⋅K⁻¹
$(thermalexpansion(English,SI2019))

julia> thermalexpansion(English,Metric) # °R⋅K⁻¹
$(thermalexpansion(English,Metric))
```
"""
@pure thermalexpansion(U::UnitSystem,S::UnitSystem) = temperature(S,U)

"""
$(convertext(:lapserate,"temperature(U,S)/length(U,S)"))

Temperature gradient over `length` or `lapserate` (K⋅m⁻¹), unit conversion factor.

```Julia
julia> lapserate(Metric,SI2019) # K⋅K⁻¹
$(lapserate(Metric,SI2019))

julia> lapserate(English,SI2019) # K⋅ft⋅°R⁻¹⋅m⁻¹
$(lapserate(English,SI2019))

julia> lapserate(English,Metric) # K⋅ft⋅°R⁻¹⋅m⁻¹
$(lapserate(English,Metric))

julia> lapserate(EnglishUS,English) # °R⋅ftUS⋅°R⁻¹⋅ft⁻¹
$(lapserate(EnglishUS,English))
```
"""
@pure lapserate(U::UnitSystem,S::UnitSystem) = unit(temperature(U,S)/length(U,S))

# molar

@doc """
$(convertext(:molarmass,"molarmass(S)/molarmass(U)"))

Molar mass or `mass` per `mole` (kg⋅mol⁻¹), unit conversion factor.

```Julia
julia> molarmass(CGS,Metric) # mol⋅mol⁻¹
$(molarmass(CGS,Metric))

julia> molarmass(Metric,SI2019) # mol⋅mol⁻¹
$(molarmass(Metric,SI2019))
```
""" molarmass(::UnitSystem,::UnitSystem)

"""
$(convertext(:molality,"molarmass(U)/molarmass(S)"))

Molality or `mole` per `mass` (mol⋅kg⁻¹), unit conversion factor.

```Julia
julia> molality(CGS,Metric) # mol⋅mol⁻¹
$(molality(CGS,Metric))

julia> molality(Metric,SI2019) # mol⋅mol⁻¹
$(molality(Metric,SI2019))
```
"""
@pure molality(U::UnitSystem,S::UnitSystem) = molarmass(S,U)

"""
$(convertext(:mole,"mass(U,S)*molality(U,S)"))

Amount of molecular substance or `mole` (mol), unit conversion factor.

```Julia
julia> mole(SI2019,Metric) # mol⋅mol⁻¹
$(mole(SI2019,Metric))

julia> mole(English,SI2019) # mol⋅slug-mol⁻¹
$(mole(English,SI2019))

julia> mole(English,SI2019)/lbm # mol⋅lb-mol⁻¹
$(mole(English,SI2019)/lbm)
```
"""
@pure mole(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)*molality(U,S))

"""
$(convertext(:molarity,"mole(U,S)/volume(U,S)"))

Molar concentration or `mole` per `volume` or `molarity` (mol⋅m⁻³), unit conversion factor.

```Julia
julia> molarity(CGS,Metric) # cm³⋅m⁻³
$(molarity(CGS,Metric))

julia> molarity(English,SI2019) # ft³⋅m⁻³
$(molarity(English,SI2019))
```
"""
@pure molarity(U::UnitSystem,S::UnitSystem) = unit(mole(U,S)/volume(U,S))

"""
$(convertext(:molarvolume,"volume(U,S)/mole(U,S)"))

Occupied `volume` per `mole` or `molarvolume` (m³⋅mol⁻¹), unit conversion factor.

```Julia
julia> molarvolume(CGS,Metric) # m³⋅cm⁻³
$(molarvolume(CGS,Metric))

julia> molarvolume(English,SI2019) # m³⋅ft⁻³
$(molarvolume(English,SI2019))
```
"""
@pure molarvolume(U::UnitSystem,S::UnitSystem) = unit(volume(U,S)/mole(U,S))

"""
$(convertext(:molarentropy,"entropy(U,S)/mole(U,S)"))

Molar heat capacity or `entropy` per `mole` (J⋅K⁻¹⋅mol⁻¹), unit conversion factor.

```Julia
julia> molarentropy(CGS,Metric) # J⋅erg⁻¹
$(molarentropy(CGS,Metric))

julia> molarentropy(English,SI2019) # J⋅°R⋅slug-mol⋅ft⁻¹⋅lb⁻¹⋅K⁻¹⋅mol⁻¹
$(molarentropy(English,SI2019))
```
"""
@pure molarentropy(U::UnitSystem,S::UnitSystem) = unit(entropy(U,S)/mole(U,S))

"""
$(convertext(:molarenergy,"energy(U,S)/mole(U,S)"))

Gibbs free `energy` per `mole` or `molarenergy` (J⋅mol⁻¹), unit conversion factor.

```Julia
julia> molarenergy(CGS,Metric) # J⋅erg⁻¹
$(molarenergy(CGS,Metric))

julia> molarenergy(English,SI2019) # J⋅slug-mol⋅ft⁻¹⋅lb⁻¹⋅mol⁻¹
$(molarenergy(English,SI2019))
```
"""
@pure molarenergy(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)/mole(U,S))

"""
$(convertext(:molarconductivity,"conductivity(U,S)*area(U,S)/mole(U,S)"))

Conductivity per `molarvolume` or `molarconductivity` (S⋅m²⋅mol⁻¹), unit conversion factor.

```Julia
julia> molarconductivity(EMU,Metric) # S⋅m²⋅abΩ⋅cm⁻²
$(molarconductivity(EMU,Metric))

julia> molarconductivity(ESU,Metric) # S⋅m²⋅statΩ⋅cm⁻²
$(molarconductivity(ESU,Metric))
```
"""
@pure molarconductivity(U::UnitSystem,S::UnitSystem) = unit(conductivity(U,S)*area(U,S)/mole(U,S))

"""
$(convertext(:molarsusceptibility,"specificsusceptibility(U,S)*molarmass(U,S)"))

Magnetic/electric molar mass `susceptibility` (m³⋅mol⁻¹), unit conversion factor.

```Julia
julia> molarsusceptibility(CGS,Metric) # m³⋅cm⁻³
$(molarsusceptibility(CGS,Metric))

julia> molarsusceptibility(Metric,SI2019) # m³⋅mol⋅mol⁻¹⋅cm⁻³
$(molarsusceptibility(Metric,SI2019))
```
"""
@pure molarsusceptibility(U::UnitSystem,S::UnitSystem) = unit(specificsusceptibility(U,S)*molarmass(U,S))

"""
$(convertext(:catalysis,"mole(U,S)/time(U,S)"))

Catalytic activity or `mole` per `time` or `catalysis` (kat, mol⋅s⁻¹), unit conversion factor.

```Julia
julia> catalysis(English,Metric) # kat⋅s⋅slug-mol⁻¹
$(catalysis(English,Metric))
```
"""
@pure catalysis(U::UnitSystem,S::UnitSystem) = unit(mole(U,S)/time(U,S))

"""
$(convertext(:specificity,"volume(U,S)/mole(U,S)/time(U,S)"))

Catalytic efficiency or `volume` per `mole` per `time` (m³⋅mol⁻¹⋅s⁻¹), unit conversion factor.

```Julia
julia> specificity(CGS,Metric) # m³⋅cm⁻³
$(specificity(CGS,Metric))

julia> specificity(English,Metric) # m³⋅slug-mol⋅mol⁻¹⋅ft⁻³
$(specificity(English,Metric))
```
"""
@pure specificity(U::UnitSystem,S::UnitSystem) = unit(volume(U,S)/mole(U,S)/time(U,S))

# photometrics

"""
$(convertext(:luminousflux,"luminousefficacy(U,S)*planck(U,S)/time(U,S)^2"))

Perceived power of light or `luminousflux` (lm, cd⋅sr), unit conversion factor.
"""
@pure luminousflux(U::UnitSystem,S::UnitSystem) = unit(frequency(U,S)^2*(luminousefficacy(S)*planckreduced(S))/(luminousefficacy(U)*planckreduced(U)))

"""
$(convertext(:luminance,"luminousflux(U,S)/area(U,S)"))

Luminous intensity per `area` or `luminance` (lx, lm⋅m⁻², cd⋅m⁻²⋅sr), unit conversion factor.

```Julia
julia> luminance(CGS,Metric) # lx⋅ph⁻¹
$(luminance(CGS,Metric))

julia> luminance(IAU,Metric) # lx⋅au²⋅lm⁻¹
$(luminance(IAU,Metric))

julia> luminance(English,Metric) # ft²⋅m⁻²
$(luminance(English,Metric))

julia> 1/10.76 # lx⋅fc⁻¹
$(1/10.76)
```
"""
@pure luminance(U::UnitSystem,S::UnitSystem) = unit(luminousflux(U,S)/area(U,S))

"""
$(convertext(:luminousenergy,"luminousflux(U,S)*time(U,S)"))

Perceived quantity of light or `luminousenergy` (lm⋅s, cd⋅s⋅sr), unit conversion factor.

```Julia
julia> luminousenergy(IAU,Metric) # s⋅day⁻¹
$(luminousenergy(IAU,Metric))
```
"""
@pure luminousenergy(U::UnitSystem,S::UnitSystem) = unit(frequency(U,S)*(luminousefficacy(S)*planckreduced(S))/(luminousefficacy(U)*planckreduced(U)))

"""
$(convertext(:luminousexposure,"luminance(U,S)*time(U,S)"))

Integrated `luminance` along `time` (lx⋅s, lm⋅s⋅m⁻², cd⋅s⋅m⁻²⋅sr), unit conversion factor.

```Julia
julia> luminousexposure(CGS,Metric) # lx⋅ph⁻¹
$(luminousexposure(CGS,Metric))

julia> luminousexposure(IAU,Metric) # s⋅au²⋅day⁻¹⋅m⁻²
$(luminousexposure(IAU,Metric))

julia> luminousexposure(English,Metric) # ft²⋅m⁻²
$(luminousexposure(English,Metric))
```
"""
@pure luminousexposure(U::UnitSystem,S::UnitSystem) = unit(luminance(U,S)*time(U,S))

@doc """
$(convertext(:luminousefficacy,"luminousefficacy(S)/luminousefficacy(U)"))

Ratio of `luminousflux` to `power` or `luminousefficacy` (lm⋅W⁻¹), unit conversion factor.

```Julia
julia> luminousefficacy(CGS,Metric) # erg⋅s⁻¹⋅W⁻¹
$(luminousefficacy(CGS,Metric))

julia> luminousefficacy(English,Metric) # ft⋅lb⋅s⁻¹⋅W⁻¹
$(luminousefficacy(English,Metric))
```
""" luminousefficacy(::UnitSystem,::UnitSystem)
