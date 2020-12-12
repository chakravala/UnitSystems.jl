
#   This file is part of UnitSystems.jl. It is licensed under the MIT license
#   UnitSystems Copyright (C) 2020 Michael Reed

"""
    kilograms(m::Real) = $(slug)m

Converts mass `m` from slugs to kilogram (kg).
"""
@pure kilograms(m::Real,U::UnitSystem=English) = mass(m,Metric,U)

"""
    slugs(m::Real) = $(1/slug)m

Converts mass `m` from kilograms to slugs (slug).
"""
@pure slugs(m::Real,U::UnitSystem=Metric) = mass(m,English,U)

"""
    feet(d) = $(1/ft)d

Converts distance `d` from meters to feet (ft).
"""
@pure feet(d,U::UnitSystem=Metric) = length(d,English,U)

"""
    meters(d) = $(ft)d

Converts distance `d` from feet to meters (m).
"""
@pure meters(d,U::UnitSystem=English) = length(d,Metric,U)

# special exact cases

# IAU to SI
@pure length(U::UnitSystem{kB,ħ,𝘤} where {kB,ħ},S::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ}) = length(U,S,1/au)
@pure time(U::UnitSystem{kB,ħ,𝘤} where {kB,ħ},S::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ}) = time(U,S,1/day)

# SI to IAU
@pure length(U::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ},S::UnitSystem{kB,ħ,𝘤} where {kB,ħ}) = length(U,S,au)
@pure time(U::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ},S::UnitSystem{kB,ħ,𝘤} where {kB,ħ}) = time(U,S,day)

# IAU to CGS
@pure length(::UnitSystem{kB,ħ,100𝘤} where {kB,ħ},S::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ}) = length(U,S,1/au)
@pure time(::UnitSystem{kB,ħ,100𝘤} where {kB,ħ},S::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ}) = time(U,S,1/day)

# CGS to IAU
@pure length(U::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ},S::UnitSystem{kB,ħ,100𝘤} where {kB,ħ}) = length(U,S,au)
@pure time(U::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ},S::UnitSystem{kB,ħ,100𝘤} where {kB,ħ}) = time(U,S,day)

# IAU to English
@pure length(U::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ},S::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ}) = length(U,S,ft/au)
@pure time(U::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ},S::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ}) = time(U,S,1/day)

# English to IAU
@pure length(U::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ},S::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ}) = length(U,S,au/ft)
@pure time(U::UnitSystem{kB,ħ,day*𝘤/au} where {kB,ħ},S::UnitSystem{kB,ħ,𝘤/ft} where {kB,ħ}) = time(U,S,day)

# base unit conversions

convertext(unit,fun) = """
```Julia
$unit(U::UnitSystem,S::UnitSystem) = $fun
$unit(v::Real,U::UnitSystem,S::UnitSystem) = v/$unit(U,S)
```
"""

# derived unit conversions

#@pure heatcapacity(U::UnitSystem,S::UnitSystem) = (boltzmann(S)/boltzmann(U))/volume(U,S)
#@pure acoustic(U::UnitSystem,S::UnitSystem) = (planckreduced(S)^3*electronmass(U)^2*lightspeed(U)^4)/(planckreduced(U)^3*electronmass(S)^2*lightspeed(S)^4)

# spacetime

@pure length(U::UnitSystem,S::UnitSystem,l=1) = unit((planckreduced(S)*electronmass(U)*lightspeed(U))/(planckreduced(U)*electronmass(S)*lightspeed(S)),l)
@doc """
$(convertext(:length,"planck(U,S)/mass(U,S)/speed(U,S)"))

Extent of one-dimensional shape or `length` (m), unit conversion factor.

```Julia
julia> length(CGS,Metric) # m⋅cm⁻¹
$(length(CGS,Metric))

julia> length(IAU,Metric) # m⋅au⁻¹
$(length(IAU,Metric))

julia> length(English,Metric) # m⋅ft⁻¹
$(length(English,Metric))

julia> length(EnglishUS,English) # ft⋅ftUS⁻¹
$(length(EnglishUS,English))

julia> length(PlanckGauss,Metric) # m⋅ℓP⁻¹
$(length(PlanckGauss,Metric))
```
""" length, ft, ftUS, au

"""
$(convertext(:area,"length(U,S)^2"))

Extent of two-dimensional shape or `area` (m²), unit conversion factor.

```Julia
julia> area(CGS,Metric) # m²⋅cm⁻²
$(area(CGS,Metric))

julia> area(English,Metric) # m²⋅ft⁻²
$(area(English,Metric))

julia> area(EnglishUS,English) # ft²⋅ftUS⁻²
$(area(EnglishUS,English))
```
"""
@pure area(U::UnitSystem,S::UnitSystem) = unit(length(U,S)^2)

"""
$(convertext(:volume,"length(U,S)^3"))

Extent of three-dimensional shape or `volume` (m³), unit conversion factor.

```Julia
julia> volume(CGS,Metric) # m³⋅cm⁻³
$(volume(CGS,Metric))

julia> volume(English,Metric) # m³⋅ft⁻³
$(volume(English,Metric))

julia> volume(EnglishUS,English) # ft³⋅ftUS⁻³
$(volume(EnglishUS,English))
```
"""
@pure volume(U::UnitSystem,S::UnitSystem) = unit(length(U,S)^3)

"""
$(convertext(:wavenumber,"1/length(U,S)"))

Number of occurences per unit of space (m⁻¹), unit conversion factor.

```Julia
julia> wavenumber(CGS,Metric) # cm⋅m⁻¹
$(wavenumber(CGS,Metric))

julia> wavenumber(English,Metric) # ft⋅m⁻¹
$(wavenumber(English,Metric))
```
"""
@pure wavenumber(U::UnitSystem,S::UnitSystem) = unit(length(S,U))

"""
$(convertext(:fuelefficiency,"1/area(U,S)"))

Distance per volume or `fuelefficiency` (m⋅m⁻³, m⁻²), unit conversion factor.

```Julia
julia> fuelefficiency(CGS,Metric) # cm²⋅m⁻²
$(fuelefficiency(CGS,Metric))

julia> fuelefficiency(English,Metric) # ft²⋅m⁻²
$(fuelefficiency(English,Metric))
```
"""
@pure fuelefficiency(U::UnitSystem,S::UnitSystem) = area(S,U)

@pure time(U::UnitSystem,S::UnitSystem,t=1) = unit(length(U,S)/lightspeed(U,S),1)
@doc """
$(convertext(:time,"length(U,S)/speed(U,S)"))

Dimension along which events are ordered or `time` (s), unit conversion factor.

```Julia
julia> time(IAU,Metric) # s⋅day⁻¹
$(time(IAU,Metric))

julia> time(PlanckGauss,Metric) # s⋅tP⁻¹
$(time(PlanckGauss,Metric))
```
""" time, day

"""
$(convertext(:frequency,"1/time(U,S)"))

Number of occurences per unit of time (Hz or s⁻¹), unit conversion factor.

```Julia
julia> frequency(IAU,Metric) day⋅s⁻¹
$(frequency(IAU,Metric))
```
"""
@pure frequency(U::UnitSystem,S::UnitSystem) = time(S,U)

"""
$(convertext(:frequencydrift,"1/time(U,S)^2"))

Drift of `frequency` per `time` or `frequencydrift` (Hz⋅s⁻¹, s⁻²), unit conversion factor.
```Julia
julia> frequencydrift(IAU,Metric) day²⋅Hz⋅s⁻¹
$(frequencydrift(IAU,Metric))
```
"""
@pure frequencydrift(U::UnitSystem,S::UnitSystem) = unit(time(S,U)^2)

"""
$(convertext(:speed,"lightspeed(S)/lightspeed(U)"))

Velocity or `length` per `time` or `speed` (m⋅s⁻¹), unit conversion factor.

```Julia
julia> speed(CGS,Metric) # m⋅cm⁻¹
$(speed(CGS,Metric))

julia> speed(IAU,Metric) # m⋅day⋅s⁻¹⋅au⁻¹
$(speed(IAU,Metric))

julia> speed(English,Metric) # m⋅ft⁻¹
$(speed(English,Metric))

julia> speed(EnglishUS,English) # ft⋅ftUS⁻¹
$(speed(EnglishUS,English))
```
"""
@pure speed(U::UnitSystem,S::UnitSystem) = lightspeed(U,S)

@doc """
$(convertext(:acceleration,"speed(U,S)/time(U,S)"))

Specific force or `speed` per `time` or `acceleration` (m⋅s⁻²), unit conversion factor.

```Julia
julia> acceleration(CGS,Metric) # m⋅s⁻¹⋅gal⁻¹
$(acceleration(CGS,Metric))

julia> acceleration(IAU,Metric) # m⋅day²⋅s⁻²⋅au⁻¹
$(acceleration(IAU,Metric))

julia> acceleration(English,Metric) # m⋅ft⁻¹
$(acceleration(English,Metric))

julia> acceleration(EnglishUS,English) # ft⋅ftUS⁻¹
$(acceleration(EnglishUS,English))
```
"""
@pure acceleration(U::UnitSystem,S::UnitSystem) = unit(speed(U,S)/time(U,S))

"""
$(convertext(:jerk,"speed(U,S)/time(U,S)^2"))

Jolt or `acceleration` per `time` or `jerk` (m⋅s⁻³), unit conversion factor.

```Julia
julia> jerk(CGS,Metric) # m⋅cm⁻¹
$(jerk(CGS,Metric))

julia> jerk(IAU,Metric) # m⋅day³⋅s⁻³⋅au⁻¹
$(jerk(IAU,Metric))

julia> jerk(English,Metric) # m⋅ft⁻¹
$(jerk(English,Metric))

julia> jerk(EnglishUS,English) # ft⋅ftUS⁻¹
$(jerk(EnglishUS,English))
```
"""
@pure jerk(U::UnitSystem,S::UnitSystem) = unit(speed(U,S)/time(U,S)^2)

"""
$(convertext(:snap,"speed(U,S)/time(U,S)^3"))

Jounce or `jerk` per `time` or `snap` (m⋅s⁻⁴), unit conversion factor.

```Julia
julia> snap(CGS,Metric) # m⋅cm⁻¹
$(snap(CGS,Metric))

julia> snap(IAU,Metric) # m⋅day⁴⋅s⁻⁴⋅au⁻¹
$(snap(IAU,Metric))

julia> snap(English,Metric) # m⋅ft⁻¹
$(snap(English,Metric))

julia> snap(EnglishUS,English) # ft⋅ftUS⁻¹
$(snap(EnglishUS,English))
```
"""
@pure snap(U::UnitSystem,S::UnitSystem) = unit(speed(U,S)/time(U,S)^3)

"""
$(convertext(:volumeflow,"area(U,S)*speed(U,S)"))

Volumetric flow rate or `volumeflow` (m³⋅s⁻¹), unit conversion factor.

```Julia
julia> volumeflow(CGS,Metric) # m³⋅cm⁻³
$(volume(CGS,Metric))

julia> volumeflow(English,Metric) # m³⋅ft⁻³
$(volume(English,Metric))

julia> volumeflow(EnglishUS,English) # ft³⋅ftUS⁻³
$(volume(EnglishUS,English))
```
"""
@pure volumeflow(U::UnitSystem,S::UnitSystem) = unit(area(U,S)*speed(U,S))

"""
$(convertext(:specificenergy,"speed(U,S)^2"))

Massic energy or `energy` per `mass` or `specificenergy` (J⋅kg⁻¹), unit conversion factor.

```Julia
julia> specificenergy(CGS,Metric) # m²⋅cm⁻²
$(specificenergy(CGS,Metric))

julia> specificenergy(IAU,Metric) # m²⋅day²⋅s⁻²⋅au⁻²
$(specificenergy(IAU,Metric))

julia> specificenergy(English,Metric) # m²⋅ft⁻²
$(specificenergy(English,Metric))

julia> specificenergy(EnglishUS,English) # ft²⋅ftUS⁻²
$(specificenergy(EnglishUS,English))
```
"""
@pure specificenergy(U::UnitSystem,S::UnitSystem) = unit(speed(U,S)^2)

# kinematic

@doc """
$(convertext(:mass,"electronmass(S)/electronmass(U)"))

Inertal `mass` or resitance to aceleration or quantity of matter (kg), unit conversion factor.

```Julia
julia> mass(CGS,Metric) # kg⋅g⁻¹
$(mass(CGS,Metric))

julia> mass(CODATA,Metric) # kg⋅kg⁻¹
$(mass(CODATA,Metric))

julia> mass(Conventional,Metric) # kg⋅kg⁻¹
$(mass(Conventional,Metric))

julia> mass(English,Metric) # kg⋅slug⁻¹
$(mass(English,Metric))

julia> mass(IAU,Metric) # kg⋅m⊙⁻¹
$(mass(IAU,Metric))

julia> mass(PlanckGauss,Metric) # kg⋅mP⁻¹
$(mass(PlanckGauss,Metric))
```
""" mass, slug, mₛ

@pure energy(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)*specificenergy(U,S))
"""
$(convertext(:energy,"mass(U,S)*specificenergy(U,S)"))

Work or heat or `energy` (J, N⋅m, kg⋅m²⋅s⁻²), unit conversion factor.

```Julia
julia> energy(CGS,Metric) # J⋅erg⁻¹
$(energy(CGS,Metric))

julia> energy(CGS,English) # ft⋅lb⋅erg⁻¹
$(energy(CGS,English))

julia> energy(English,Metric) # J⋅ft⁻¹⋅lb⁻¹
$(energy(English,Metric))

julia> 0.001/3600 # J⋅kW⁻¹⋅h⁻¹
$(0.001/3600)

julia> 1/kcal # J⋅kcalₜₕ⁻¹
$(1/kcal)

julia> 1/BTUJ # J⋅BTU⁻¹
$(1/1055.036345118633)

julia> 1/charge(SI2019) # J⋅eV⁻¹
$(1/charge(SI2019))
```
""" energy, cal, kcal, BTUJ, BTUftlb

"""
$(convertext(:power,"energy(U,S)/time(U,S))"))

Radiant flux or `power` or `energy` per `time` (W, J⋅s⁻¹, kg⋅m²⋅s⁻³), unit conversion factor.

```Julia
julia> power(CGS,Metric) # W⋅s⋅erg⁻¹
$(power(CGS,Metric))

julia> power(English,Metric) # W⋅s⋅ft⁻¹⋅lb⁻¹
$(power(English,Metric))
```
"""
@pure power(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)/time(U,S))

@pure force(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)*acceleration(U,S))
"""
$(convertext(:force,"mass(U,S)*acceleration(U,S)"))

Weight or force or `mass` times `acceleration` (N, kg⋅m⋅s⁻²), unit conversion factor.

```Julia
julia> force(CGS,Metric) # N⋅dyn⁻¹
$(force(CGS,Metric))

julia> force(CGS,English) # lb⋅dyn⁻¹
$(force(CGS,English))

julia> force(English,Metric) # N⋅lb⁻¹
$(force(English,Metric))

julia> force(1/lbm,Metric,English) # pdl⋅N⁻¹
$(force(1/lbm,Metric,English))

julia> g₀ # kp⋅N⁻¹
$g₀
```
""" force, g₀, g0, lbm

@pure pressure(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)/length(U,S)/time(U,S)^2)
@doc """
$(convertext(:pressure,"energy(U,S)/volume(U,S)"))

Pressure or stress or `force` per `area` (Pa, N⋅m⁻², kg⋅m⁻¹⋅s⁻²), unit conversion factor.

```Julia
julia> pressure(CGS,Metric) # Pa⋅Ba⁻¹
$(pressure(CGS,Metric))

julia> 1/atm # Pa⋅atm⁻¹
$(1/atm)

julia> pressure(English,Metric) # Pa⋅ft²⋅lb⁻¹
$(pressure(English,Metric))

julia> pressure(Metric,English)/12^2 # psi⋅Pa⁻¹
$(pressure(Metric,English)/12^2)
```
""" pressure, atm

# mechanical

"""
$(convertext(:momentum,"mass(U,S)*speed(U,S)"))

Linear `momentum` or `mass` times `speed` (N⋅s, kg⋅m⋅s⁻¹), unit conversion factor.

```Julia
julia> momentum(CGS,Metric) # N⋅dyn⁻¹
$(momentum(CGS,Metric))

julia> momentum(CGS,English) # lb⋅dyn⁻¹
$(momentum(CGS,English))

julia> momentum(English,Metric) # N⋅lb⁻¹
$(momentum(English,Metric))
```
"""
@pure momentum(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)*speed(U,S))

"""
$(convertext(:angularmomentum,"momentum(U,S)*length(U,S)"))

Rotational momentum or `angularmomentum` (N⋅m⋅s, kg⋅m²⋅s⁻¹), unit conversion factor.

```Julia
julia> momentum(CGS,Metric) # N⋅m⋅dyn⁻¹⋅cm⁻¹
$(momentum(CGS,Metric))

julia> momentum(CGS,English) # lb⋅ft⋅dyn⁻¹⋅cm⁻¹
$(momentum(CGS,English))

julia> momentum(English,Metric) # N⋅m⋅lb⁻¹⋅ft⁻¹
$(momentum(English,Metric))
```
"""
@pure angularmomentum(U::UnitSystem,S::UnitSystem) = unit(momentum(U,S)*length(U,S))

"""
$(convertext(:yank,"mass(U,S)*jerk(U,S)"))

Rate of change of `force` or `yank` (N⋅s⁻¹, kg⋅m⋅s⁻³), unit conversion factor.

```Julia
julia> yank(CGS,Metric) # N⋅dyn⁻¹
$(yank(CGS,Metric))

julia> yank(CGS,English) # lb⋅dyn⁻¹
$(yank(CGS,English))

julia> yank(English,Metric) # N⋅lb⁻¹⋅
$(yank(English,Metric))
```
"""
@pure yank(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)*jerk(U,S))

"""
$(convertext(:areadensity,"mass(U,S)/area(U,S)"))

Surface or `areadensity` or `mass` per `area` (kg⋅m⁻²), unit conversion factor.

```Julia
julia> areadensity(CGS,Metric) # kg⋅cm²⋅g⁻¹⋅m⁻²
$(areadensity(CGS,Metric))

julia> areadensity(CGS,English) # slug⋅cm²⋅g⁻¹⋅ft⁻²
$(areadensity(CGS,English))

julia> areadensity(English,Metric) # kg⋅ft²⋅slug⁻¹⋅m⁻²
$(areadensity(English,Metric))
```
"""
@pure areadensity(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)/area(U,S))

"""
$(convertext(:density,"mass(U,S)/volume(U,S)"))

Specific mass or `mass` per `volume` or `density` (kg⋅m⁻³), unit conversion factor.

```Julia
julia> density(CGS,Metric) # kg⋅cm³⋅g⁻¹⋅m⁻³
$(density(CGS,Metric))

julia> density(CGS,English) # slug⋅cm³⋅g⁻¹⋅ft⁻³
$(density(CGS,English))

julia> density(English,Metric) # kg⋅ft³⋅slug⁻¹⋅m⁻³
$(density(English,Metric))
```
"""
@pure density(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)/volume(U,S))

"""
$(convertext(:specificvolume,"volume(U,S)/mass(U,S)"))

Reciprocal `density` or `volume` per `mass` or `specificvolume` (m³⋅kg), unit conversion factor.

```Julia
julia> specificvolume(CGS,Metric) # g⋅m³⋅kg⁻¹⋅cm⁻³
$(specificvolume(CGS,Metric))

julia> specificvolume(CGS,English) # kg⋅ft³⋅slug⁻¹⋅cm⁻³
$(specificvolume(CGS,English))

julia> specificvolume(English,Metric) # slug⋅m³⋅kg⁻¹⋅ft⁻³
$(specificvolume(English,Metric))
```
"""
@pure specificvolume(U::UnitSystem,S::UnitSystem) = unit(volume(U,S)/mass(U,S))

"""
$(convertext(:action,"energy(U,S)*time(U,S)"))

Integrated `momentum` over `length` or `action` (J⋅s, N⋅m⋅s), unit conversion factor.

```Julia
julia> action(CGS,Metric) # J⋅erg⁻¹
$(action(CGS,Metric))

julia> action(CGS,English) # ft⋅lb⋅erg⁻¹
$(action(CGS,English))

julia> action(English,Metric) # J⋅ft⁻¹⋅lb⁻¹
$(action(English,Metric))
```
"""
@pure action(U::UnitSystem,S::UnitSystem) = unit(momentum(U,S)*length(U,S))

"""
$(convertext(:stiffness,"energy(U,S)/area(U,S)"))

Amount of `force` per `length` or `stiffness` (N⋅m⁻¹, J⋅m⁻², kg⋅s⁻²), unit conversion factor.

```Julia
julia> stiffness(CGS,Metric) # kg⋅g⁻¹
$(stiffness(CGS,Metric))

julia> stiffness(CGS,English) # slug⋅g⁻¹
$(stiffness(CGS,English))

julia> stiffness(English,Metric) # kg⋅slug⁻¹
$(stiffness(English,Metric))
```
"""
@pure stiffness(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)/area(U,S))

"""
$(convertext(:irradiance,"power(U,S)/area(U,S)"))

Heat flux density or `power` per `area` or `irradiance` (W⋅m⁻², kg⋅s⁻³), unit conversion factor.

```Julia
julia> irradiance(CGS,Metric) # kg⋅g⁻¹
$(irradiance(CGS,Metric))

julia> irradiance(CGS,English) # slug⋅g⁻¹
$(irradiance(CGS,English))

julia> irradiance(English,Metric) # kg⋅slug⁻¹
$(irradiance(English,Metric))
```
"""
@pure irradiance(U::UnitSystem,S::UnitSystem) = unit(power(U,S)/area(U,S))

"""
$(convertext(:diffusivity,"planck(U,S)/mass(U,S)"))

Thermal `diffusivity` or kinematic viscostiy (m²⋅s⁻¹), unit conversion factor.

```Julia
julia> diffusivity(CGS,Metric) # m²⋅cm⁻²
$(diffusivity(CGS,Metric))

julia> diffusivity(English,Metric) # m²⋅ft⁻²
$(diffusivity(English,Metric))

julia> diffusivity(EnglishUS,English) # ft²⋅ftUS⁻²
$(diffusivity(EnglishUS,English))
```
"""
@pure diffusivity(U::UnitSystem,S::UnitSystem) = unit((planckreduced(S)*electronmass(U))/(planckreduced(U)*electronmass(S)))

"""
$(convertext(:viscosity,"mass(U,S)/length(U,S)/time(U,S)^2"))

Resistance to deformation or `viscosity` (Pa⋅s, kg⋅m⁻¹⋅s⁻¹), unit conversion factor.

```Julia
julia> viscosity(CGS,Metric) # Pa⋅Ba⁻¹
$(viscosity(CGS,Metric))

julia> viscosity(English,Metric) # Pa⋅ft²⋅lb⁻¹
$(viscosity(English,Metric))
```
"""
@pure viscosity(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)/length(U,S)/time(U,S))

"""
$(convertext(:lineardensity,"mass(U,S)/length(U,S)"))

Amount of `lineardensity` or `mass` per `length` (kg⋅m⁻¹), unit conversion factor.

```Julia
julia> lineardensity(CGS,Metric) # kg⋅cm¹⋅g⁻¹⋅m⁻¹
$(lineardensity(CGS,Metric))

julia> lineardensity(CGS,English) # slug⋅cm¹⋅g⁻¹⋅ft⁻¹
$(lineardensity(CGS,English))

julia> lineardensity(English,Metric) # kg⋅ft¹⋅slug⁻¹⋅m⁻¹
$(lineardensity(English,Metric))
```
"""
@pure lineardensity(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)/length(U,S))

"""
$(convertext(:massflow,"mass(U,S)/time(U,S)"))

Rate of `massflow` or `mass` per `time` (kg⋅s⁼¹), unit conversion factor.

```Julia
julia> massflow(CGS,Metric) # kg⋅g⁻¹
$(massflow(CGS,Metric))

julia> massflow(CODATA,Metric) # kg⋅kg⁻¹
$(massflow(CODATA,Metric))

julia> massflow(Conventional,Metric) # kg⋅kg⁻¹
$(massflow(Conventional,Metric))

julia> massflow(English,Metric) # kg⋅slug⁻¹
$(massflow(English,Metric))
```
"""
@pure massflow(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)/time(U,S))

"""
$(convertext(:radiantflux,"power(U,S)/length(U,S)"))

Spectral power or `radiantflux` (W⋅m⁻¹), unit conversion factor.

```Julia
julia> radiantflux(CGS,Metric) # kg⋅m⋅g⁻¹⋅cm⁻¹
$(radiantflux(CGS,Metric))

julia> radiantflux(CGS,English) # slug⋅ft⋅g⁻¹⋅cm⁻¹
$(radiantflux(CGS,English))

julia> radiantflux(English,Metric) # kg⋅m⋅slug⁻¹⋅ft⁻¹
$(radiantflux(English,Metric))
```
"""
@pure radiantflux(U::UnitSystem,S::UnitSystem) = unit(power(U,S)/length(U,S))

"""
$(convertext(:powerdensity,"power(U,S)/volume(U,S)"))

Spectral irradiance (volume) or `powerdensity` (W⋅m⁻³), unit conversion factor.

```Julia
julia> powerdensity(CGS,Metric) # kg⋅cm⋅g⁻¹⋅m⁻¹
$(powerdensity(CGS,Metric))

julia> powerdensity(CGS,English) # slug⋅cm⋅g⁻¹⋅ft⁻¹
$(powerdensity(CGS,English))

julia> powerdensity(English,Metric) # kg⋅ft⋅slug⁻¹⋅m⁻¹
$(powerdensity(English,Metric))
```
"""
@pure powerdensity(U::UnitSystem,S::UnitSystem) = unit(power(U,S)/volume(U,S))

"""
$(convertext(:compressibility,"1/pressure(U,S)"))

Relative volume change or `compressibility` (Pa⁻¹), unit conversion factor.

```Julia
julia> compressibility(CGS,Metric) # Ba⋅Pa⁻¹
$(compressibility(CGS,Metric))

julia> compressibility(English,Metric) # lb⋅ft⁻²⋅Pa⁻¹
$(compressibility(English,Metric))

julia> compressibility(Metric,English)/12^2 # Pa⋅psi⁻¹
$(compressibility(Metric,English)/12^2)
```
"""
@pure compressibility(U::UnitSystem,S::UnitSystem) = pressure(S,U)

"""
$(convertext(:fluence,"energy(U,S)/area(U,S"))

Radiant exposure or `energy` per `area` or `fluence` (J⋅m⁻²), unit conversion factor.

```Julia
julia> fluence(CGS,Metric) # kg⋅g⁻¹
$(mass(CGS,Metric))

julia> fluence(CODATA,Metric) # kg⋅kg⁻¹
$(mass(CODATA,Metric))

julia> fluence(Conventional,Metric) # kg⋅kg⁻¹
$(fluence(Conventional,Metric))

julia> fluence(English,Metric) # kg⋅slug⁻¹
$(fluence(English,Metric))
```
"""
@pure fluence(U::UnitSystem,S::UnitSystem) = unit(energy(U,S)/area(U,S))

"""
$(convertext(:rotationalinertia,"mass(U,S)*area(U,S)"))

Moment of inertia or `rotationalinertia` (kg⋅m²), unit conversion factor.

```Julia
julia> rotationalinertia(CGS,Metric) # kg⋅m²⋅g⁻¹⋅cm⁻²
$(rotationalinertia(CGS,Metric))

julia> rotationalinertia(CGS,English) # slug⋅ft²⋅g⁻¹⋅cm⁻²
$(rotationalinertia(CGS,English))

julia> rotationalinertia(English,Metric) # kg⋅m²⋅slug⁻¹⋅ft⁻²
$(rotationalinertia(English,Metric))
```
"""
@pure rotationalinertia(U::UnitSystem,S::UnitSystem) = unit(mass(U,S)*area(U,S))
