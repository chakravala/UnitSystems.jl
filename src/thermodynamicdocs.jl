
#   This file is part of UnitSystems.jl
#   It is licensed under the MIT license
#   UnitSystems Copyright (C) 2020 Michael Reed
#       _           _                         _
#      | |         | |                       | |
#   ___| |__   __ _| | ___ __ __ ___   ____ _| | __ _
#  / __| '_ \ / _` | |/ / '__/ _` \ \ / / _` | |/ _` |
# | (__| | | | (_| |   <| | | (_| |\ V / (_| | | (_| |
#  \___|_| |_|\__,_|_|\_\_|  \__,_| \_/ \__,_|_|\__,_|
#
#   https://github.com/chakravala
#   https://crucialflow.com

# thermodynamics

@doc """
$(convertext(:temperature,"mass(U,S)*speed(U,S)^2/entropy(U,S)"))

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
""" temperature, Θ

@doc """
$(convertext(:entropy,"energy(U,S)/temperature(U,S)"))

Heat capacity or `energy` per `temperature` or `entropy` (J⋅K⁻¹), unit conversion factor.

```Julia
julia> entropy(Metric,SI2019) # K⋅K⁻¹
$(entropy(Metric,SI2019))

julia> entropy(CGS,Metric) # J⋅erg⁻¹
$(entropy(CGS,Metric))

julia> entropy(English,SI2019) # J⋅°R⋅K⁻¹⋅ft⁻¹⋅lb⁻¹
$(entropy(English,SI2019))

julia> entropy(Survey,English) # ftUS²⋅°R⋅°ft⁻²⋅°R⁻¹
$(entropy(Survey,English))
```
""" entropy

@doc """
$(convertext(:specificentropy,"specificenergy(U,S)/temperature(U,S)"))

Specific heat capacity or `specificentropy` (J⋅K⁻¹⋅kg⁻¹), unit conversion factor.

```Julia
julia> specificentropy(Metric,SI2019) # m²⋅K⋅K⁻¹⋅cm⁻²
$(specificentropy(Metric,SI2019))

julia> specificentropy(CGS,Metric) # m²⋅cm⁻²
$(specificentropy(CGS,Metric))

julia> specificentropy(English,Metric) # m²⋅°R⋅K⁻¹⋅ft⁻²
$(specificentropy(English,Metric))

julia> specificentropy(Survey,English) # ft²⋅°R⋅ftUS⁻²⋅°R⁻¹
$(specificentropy(Survey,English))
```
""" specificentropy

@doc """
$(convertext(:volumeheatcapacity,"entropy(U,S)/volume(U,S)"))

The `entropy` per `volume` or `volumeheatcapacity` (J⋅K⁻¹⋅m⁻³), unit conversion factor.

```Julia
julia> volumeheatcapacity(Metric,SI2019) # K⋅K⁻¹
$(volumeheatcapacity(Metric,SI2019))

julia> volumeheatcapacity(CGS,Metric) # J⋅cm³⋅erg⁻¹⋅m⁻³
$(volumeheatcapacity(CGS,Metric))

julia> volumeheatcapacity(English,SI2019) # J⋅ft²⋅°R⋅K⁻¹⋅lb⁻¹⋅m⁻³
$(volumeheatcapacity(English,SI2019))

julia> volumeheatcapacity(Survey,English) # ftUS⁵°R⋅°ft⁻⁵⋅°R⁻¹
$(volumeheatcapacity(Survey,English))
```
""" volumeheatcapacity

@doc """
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
""" thermalconductivity

@doc """
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
""" thermalconductance

@doc """
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
""" thermalresistance

@doc """
$(convertext(:thermalresistivity,"1/thermalconductivity(U,S)"))

Resistance to heat flow or `thermalresistance` (K⋅W⁻¹), unit conversion factor.

```Julia
julia> thermalresistance(Metric,SI2019) # K⋅K⁻¹
$(thermalresistance(Metric,SI2019))

julia> thermalresistance(CGS,Metric) # erg⋅s⁻¹⋅W⁻¹
$(thermalresistance(CGS,Metric))

julia> thermalresistance(English,Metric) # ft⋅lb⋅K⋅°R⁻¹⋅J⁻¹
$(thermalresistance(English,Metric))
```
""" thermalresistivity


@doc """
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
""" thermalexpansion

@doc """
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
""" lapserate

# molar

@doc """
$(convertext(:molarmass,"molarmass(S)/molarmass(U)"))

Molar mass or `mass` per `mole` (kg⋅mol⁻¹), unit conversion factor.

```Julia
julia> molarmass(CGS,Metric) # kg⋅g⁻¹
$(molarmass(CGS,Metric))

julia> molarmass(Metric,SI2019) # mol⋅mol⁻¹
$(molarmass(Metric,SI2019))
```
""" molarmass(::UnitSystem,::UnitSystem)

@doc """
$(convertext(:molality,"molarmass(U)/molarmass(S)"))

Molality or `mole` per `mass` (mol⋅kg⁻¹), unit conversion factor.

```Julia
julia> molality(CGS,Metric) # kg⋅g⁻¹
$(molality(CGS,Metric))

julia> molality(Metric,SI2019) # mol⋅mol⁻¹
$(molality(Metric,SI2019))
```
""" molality

@doc """
$(convertext(:molaramount,"mass(U,S)*molality(U,S)"))

Amount of molecular substance or `molaramount` (mol), unit conversion factor.

```Julia
julia> molaramount(SI2019,Metric) # mol⋅mol⁻¹
$(molaramount(SI2019,Metric))

julia> molaramount(British,SI2019) # mol⋅slug-mol⁻¹
$(molaramount(English,SI2019))

julia> molaramount(English,SI2019) # mol⋅lb-mol⁻¹
$(molaramount(English,SI2019))
```
""" molaramount, N
#julia> molaramount(English,SI2019)*ft/g₀ # mol⋅lb-mol⁻¹
#$(molaramount(English,SI2019)*ft/g₀)

@doc """
$(convertext(:molarity,"molaramount(U,S)/volume(U,S)"))

Molar concentration or `molaramount` per `volume` (mol⋅m⁻³), unit conversion factor.

```Julia
julia> molarity(CGS,Metric) # cm³⋅m⁻³
$(molarity(CGS,Metric))

julia> molarity(English,SI2019) # ft³⋅m⁻³
$(molarity(English,SI2019))
```
""" molarity

@doc """
$(convertext(:molarvolume,"volume(U,S)/molaramount(U,S)"))

Occupied `volume` per `molaramount` or `molarvolume` (m³⋅mol⁻¹), unit conversion factor.

```Julia
julia> molarvolume(CGS,Metric) # m³⋅cm⁻³
$(molarvolume(CGS,Metric))

julia> molarvolume(English,SI2019) # m³⋅ft⁻³
$(molarvolume(English,SI2019))
```
""" molarvolume

@doc """
$(convertext(:molarentropy,"entropy(U,S)/molaramount(U,S)"))

Molar heat capacity or `entropy` per `molaramount` (J⋅K⁻¹⋅mol⁻¹), unit conversion factor.

```Julia
julia> molarentropy(CGS,Metric) # J⋅erg⁻¹
$(molarentropy(CGS,Metric))

julia> molarentropy(English,SI2019) # J⋅°R⋅lb-mol⋅ft⁻¹⋅lb⁻¹⋅K⁻¹⋅mol⁻¹
$(molarentropy(English,SI2019))
```
""" molarentropy

@doc """
$(convertext(:molarenergy,"energy(U,S)/molaramount(U,S)"))

Gibbs free `energy` per `mole` or `molarenergy` (J⋅mol⁻¹), unit conversion factor.

```Julia
julia> molarenergy(CGS,Metric) # J⋅erg⁻¹
$(molarenergy(CGS,Metric))

julia> molarenergy(English,SI2019) # J⋅slug-mol⋅ft⁻¹⋅lb⁻¹⋅mol⁻¹
$(molarenergy(English,SI2019))
```
""" molarenergy

@doc """
$(convertext(:molarconductivity,"conductivity(U,S)*area(U,S)/molaramount(U,S)"))

Conductivity per `molarvolume` or `molarconductivity` (S⋅m²⋅mol⁻¹), unit conversion factor.

```Julia
julia> molarconductivity(EMU,Metric) # S⋅m²⋅abΩ⋅cm⁻²
$(molarconductivity(EMU,Metric))

julia> molarconductivity(ESU,Metric) # S⋅m²⋅statΩ⋅cm⁻²
$(molarconductivity(ESU,Metric))
```
""" molarconductivity

@doc """
$(convertext(:molarsusceptibility,"specificsusceptibility(U,S)*molarmass(U,S)"))

Magnetic/electric molar mass `susceptibility` (m³⋅mol⁻¹), unit conversion factor.

```Julia
julia> molarsusceptibility(CGS,Metric) # m³⋅cm⁻³
$(molarsusceptibility(CGS,Metric))

julia> molarsusceptibility(Metric,SI2019) # m³⋅mol⋅mol⁻¹⋅cm⁻³
$(molarsusceptibility(Metric,SI2019))
```
""" molarsusceptibility

@doc """
$(convertext(:catalysis,"molaramount(U,S)/time(U,S)"))

Catalytic activity or `molaramount` per `time` (kat, mol⋅s⁻¹), unit conversion factor.

```Julia
julia> catalysis(English,Metric) # kat⋅s⋅lb-mol⁻¹
$(catalysis(English,Metric))
```
""" catalysis

@doc """
$(convertext(:specificity,"volume(U,S)/molaramount(U,S)/time(U,S)"))

Catalytic efficiency or `volume` per `mole` per `time` (m³⋅mol⁻¹⋅s⁻¹), unit conversion factor.

```Julia
julia> specificity(CGS,Metric) # m³⋅cm⁻³
$(specificity(CGS,Metric))

julia> specificity(English,Metric) # m³⋅lb-mol⋅mol⁻¹⋅ft⁻³
$(specificity(English,Metric))
```
""" specificity

@doc """
$(convertext(:diffusionflux,"molaramount(U,S)*photonirradiance(U,S)"))

Molar diffusion flux or `molarmount` times `flux` (mol⋅s⁻¹⋅m⁻²), unit conversion factor.

```Julia
julia> diffusionflux(CGS,Metric) # cm²⋅m⁻²
$(diffusionflux(CGS,Metric))

julia> diffusionflux(English,Metric) # ft²⋅mol⋅lb-mol⁻¹⋅m⁻²
$(diffusionflux(English,Metric))
```
""" diffusionflux

# photometrics

@doc """
$(convertext(:luminousflux,"luminousenergy(U,S)*frequency(U,S)"))

Perceived power of light or `luminousflux` (lm, cd⋅rad⋅²), unit conversion factor.
""" luminousflux, J

@doc """
$(convertext(:luminousintensity,"luminousflux(U,S)/solidangle(U,S)"))

Perceived power of light or `luminousintensity` (cd, lm⋅rad⁻²), unit conversion factor.
""" luminousintensity

@doc """
$(convertext(:luminance,"luminousintensity(U,S)/area(U,S)"))

Luminous intensity per `area` or `luminance` (cd⋅m⁻², lm⋅m⁻²⋅rad⁻²), unit conversion factor.

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
""" luminance

@doc """
$(convertext(:illuminance,"luminousflux(U,S)/area(U,S)"))

Luminous flux per `area` or `luminance` (lx, lm⋅m⁻², cd⋅m⁻²⋅rad²), unit conversion factor.

```Julia
julia> illuminance(CGS,Metric) # lx⋅ph⁻¹
$(illuminance(CGS,Metric))

julia> illuminance(IAU,Metric) # lx⋅au²⋅lm⁻¹
$(illuminance(IAU,Metric))

julia> illuminance(English,Metric) # ft²⋅m⁻²
$(illuminance(English,Metric))

julia> 1/10.76 # lx⋅fc⁻¹
$(1/10.76)
```
""" illuminance

@doc """
$(convertext(:luminousenergy,"luminousflux(U,S)*time(U,S)"))

Perceived quantity of light or `luminousenergy` (lm⋅s, cd⋅s⋅sr), unit conversion factor.

```Julia
julia> luminousenergy(IAU,Metric) # s⋅day⁻¹
$(luminousenergy(IAU,Metric))
```
""" luminousenergy

@doc """
$(convertext(:luminousexposure,"illuminance(U,S)*time(U,S)"))

Integrated `luminance` along `time` (lx⋅s, lm⋅s⋅m⁻², cd⋅s⋅m⁻²⋅sr), unit conversion factor.

```Julia
julia> luminousexposure(CGS,Metric) # lx⋅ph⁻¹
$(luminousexposure(CGS,Metric))

julia> luminousexposure(IAU,Metric) # s⋅au²⋅day⁻¹⋅m⁻²
$(luminousexposure(IAU,Metric))

julia> luminousexposure(English,Metric) # ft²⋅m⁻²
$(luminousexposure(English,Metric))
```
""" luminousexposure

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
