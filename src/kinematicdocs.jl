
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

@doc """
$(convertext(:angle,"angle(U,S)"))

Extent of one-dimensional angle or `angle` (rad), unit conversion factor.

```Julia
julia> angle(CGS,Metric) # rad⋅rad⁻¹
$(angle(CGS,Metric))
```
""" angle, A

@doc """
$(convertext(:solidangle,"angle(U,S)^2"))

Extent of two-dimensional angle or `solidangle` (rad²), unit conversion factor.

```Julia
julia> solidangle(CGS,Metric) # rad²⋅rad⁻²
$(solidangle(CGS,Metric))
```
""" solidangle

# spacetime

@doc """
$(convertext(:area,"length(U,S)^2"))

Extent of two-dimensional shape or `area` (m²), unit conversion factor.

```Julia
julia> area(CGS,Metric) # m²⋅cm⁻²
$(area(CGS,Metric))

julia> area(English,Metric) # m²⋅ft⁻²
$(area(English,Metric))

julia> area(Survey,English) # ft²⋅ftUS⁻²
$(area(Survey,English))
```
""" area

@doc """
$(convertext(:volume,"length(U,S)^3"))

Extent of three-dimensional shape or `volume` (m³), unit conversion factor.

```Julia
julia> volume(CGS,Metric) # m³⋅cm⁻³
$(volume(CGS,Metric))

julia> volume(English,Metric) # m³⋅ft⁻³
$(volume(English,Metric))

julia> volume(Survey,English) # ft³⋅ftUS⁻³
$(volume(Survey,English))
```
""" volume

@doc """
$(convertext(:wavenumber,"1/length(U,S)"))

Number of occurences per unit of space (m⁻¹), unit conversion factor.

```Julia
julia> wavenumber(CGS,Metric) # cm⋅m⁻¹
$(wavenumber(CGS,Metric))

julia> wavenumber(English,Metric) # ft⋅m⁻¹
$(wavenumber(English,Metric))
```
""" wavenumber

@doc """
$(convertext(:angularwavenumber,"angle(U,S)/length(U,S)"))

Number of occurences per unit of space (m⁻¹), unit conversion factor.

```Julia
julia> angularwavenumber(CGS,Metric) # cm⋅m⁻¹
$(angularwavenumber(CGS,Metric))

julia> angularwavenumber(English,Metric) # ft⋅m⁻¹
$(angularwavenumber(English,Metric))
```
""" angularwavenumber

@doc """
$(convertext(:fuelefficiency,"1/area(U,S)"))

Distance per volume or `fuelefficiency` (m⋅m⁻³, m⁻²), unit conversion factor.

```Julia
julia> fuelefficiency(CGS,Metric) # cm²⋅m⁻²
$(fuelefficiency(CGS,Metric))

julia> fuelefficiency(English,Metric) # ft²⋅m⁻²
$(fuelefficiency(English,Metric))
```
""" fuelefficiency

@doc """
$(convertext(:numberdensity,"1/volume(U,S)"))

Number per `volume` or `numberdensity` (m⁻³ or ft⁻³), unit conversion factor.

```Julia
julia> numberdensity(CGS,Metric) # cm³⋅m⁻³
$(numberdensity(CGS,Metric))

julia> numberdensity(English,Metric) # ft³⋅m⁻³
$(numberdensity(English,Metric))
```
""" numberdensity

@doc """
$(convertext(:frequency,"1/time(U,S)"))

Number of occurences per unit of time (Hz or s⁻¹), unit conversion factor.

```Julia
julia> frequency(IAU,Metric) day⋅s⁻¹
$(frequency(IAU,Metric))
```
""" frequency

@doc """
$(convertext(:angularfrequency,"angle(U,S)/time(U,S)"))

Circular radian frequency (rad⋅Hz or rad⋅s⁻¹), unit conversion factor.

```Julia
julia> angularfrequency(IAU,Metric) day⋅s⁻¹
$(frequency(IAU,Metric))
```
""" angularfrequency

@doc """
$(convertext(:frequencydrift,"1/time(U,S)^2"))

Drift of `frequency` per `time` or `frequencydrift` (Hz⋅s⁻¹, s⁻²), unit conversion factor.
```Julia
julia> frequencydrift(IAU,Metric) day²⋅Hz⋅s⁻¹
$(frequencydrift(IAU,Metric))
```
""" frequencydrift

@doc """
$(convertext(:speed,"lightspeed(S)/lightspeed(U)"))

Velocity or `length` per `time` or `speed` (m⋅s⁻¹), unit conversion factor.

```Julia
julia> speed(CGS,Metric) # m⋅cm⁻¹
$(speed(CGS,Metric))

julia> speed(IAU,Metric) # m⋅day⋅s⁻¹⋅au⁻¹
$(speed(IAU,Metric))

julia> speed(English,Metric) # m⋅ft⁻¹
$(speed(English,Metric))

julia> speed(Survey,English) # ft⋅ftUS⁻¹
$(speed(Survey,English))
```
""" speed

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

julia> acceleration(Survey,English) # ft⋅ftUS⁻¹
$(acceleration(Survey,English))
```
""" acceleration

@doc """
$(convertext(:jerk,"speed(U,S)/time(U,S)^2"))

Jolt or `acceleration` per `time` or `jerk` (m⋅s⁻³), unit conversion factor.

```Julia
julia> jerk(CGS,Metric) # m⋅cm⁻¹
$(jerk(CGS,Metric))

julia> jerk(IAU,Metric) # m⋅day³⋅s⁻³⋅au⁻¹
$(jerk(IAU,Metric))

julia> jerk(English,Metric) # m⋅ft⁻¹
$(jerk(English,Metric))

julia> jerk(Survey,English) # ft⋅ftUS⁻¹
$(jerk(Survey,English))
```
""" jerk

@doc """
$(convertext(:snap,"speed(U,S)/time(U,S)^3"))

Jounce or `jerk` per `time` or `snap` (m⋅s⁻⁴), unit conversion factor.

```Julia
julia> snap(CGS,Metric) # m⋅cm⁻¹
$(snap(CGS,Metric))

julia> snap(IAU,Metric) # m⋅day⁴⋅s⁻⁴⋅au⁻¹
$(snap(IAU,Metric))

julia> snap(English,Metric) # m⋅ft⁻¹
$(snap(English,Metric))

julia> snap(Survey,English) # ft⋅ftUS⁻¹
$(snap(Survey,English))
```
""" snap

@doc """
$(convertext(:crackle,"speed(U,S)/time(U,S)^4"))

A `snap` per `time` or `crackle` (m⋅s⁻⁵), unit conversion factor.

```Julia
julia> crackle(CGS,Metric) # m⋅cm⁻¹
$(crackle(CGS,Metric))

julia> crackle(IAU,Metric) # m⋅day⁵⋅s⁻⁵⋅au⁻¹
$(crackle(IAU,Metric))

julia> crackle(English,Metric) # m⋅ft⁻¹
$(crackle(English,Metric))

julia> crackle(Survey,English) # ft⋅ftUS⁻¹
$(crackle(Survey,English))
```
""" crackle

@doc """
$(convertext(:pop,"speed(U,S)/time(U,S)^5"))

A `crackle` per `time` or `pop` (m⋅s⁻⁶), unit conversion factor.

```Julia
julia> pop(CGS,Metric) # m⋅cm⁻¹
$(pop(CGS,Metric))

julia> pop(IAU,Metric) # m⋅day⁶⋅s⁻⁶⋅au⁻¹
$(pop(IAU,Metric))

julia> pop(English,Metric) # m⋅ft⁻¹
$(pop(English,Metric))

julia> pop(Survey,English) # ft⋅ftUS⁻¹
$(pop(Survey,English))
```
""" pop

@doc """
$(convertext(:volumeflow,"area(U,S)*speed(U,S)"))

Volumetric flow rate or `volumeflow` (m³⋅s⁻¹), unit conversion factor.

```Julia
julia> volumeflow(CGS,Metric) # m³⋅cm⁻³
$(volume(CGS,Metric))

julia> volumeflow(English,Metric) # m³⋅ft⁻³
$(volume(English,Metric))

julia> volumeflow(Survey,English) # ft³⋅ftUS⁻³
$(volume(Survey,English))
```
""" volumeflow

@doc """
$(convertext(:specificenergy,"speed(U,S)^2"))

Massic energy or `energy` per `mass` or `specificenergy` (J⋅kg⁻¹), unit conversion factor.

```Julia
julia> specificenergy(CGS,Metric) # m²⋅cm⁻²
$(specificenergy(CGS,Metric))

julia> specificenergy(IAU,Metric) # m²⋅day²⋅s⁻²⋅au⁻²
$(specificenergy(IAU,Metric))

julia> specificenergy(English,Metric) # m²⋅ft⁻²
$(specificenergy(English,Metric))

julia> specificenergy(Survey,English) # ft²⋅ftUS⁻²
$(specificenergy(Survey,English))
```
""" specificenergy

# kinematic

@doc """
$(convertext(:mass,"electronmass(S)/electronmass(U)"))

Inertal `mass` or matter quantity or resistance to aceleration (kg), unit conversion factor.

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
""" mass, M

@doc """
$(convertext(:inertia,"mass(U,S)/gravity(U,S)"))

Inertal `mass` or matter quantity or resistance to aceleration (kg), unit conversion factor.

```Julia
julia> inertia(CGS,Metric) # kg⋅g⁻¹
$(inertia(CGS,Metric))

julia> inertia(CODATA,Metric) # kg⋅kg⁻¹
$(inertia(CODATA,Metric))

julia> inertia(Conventional,Metric) # kg⋅kg⁻¹
$(inertia(Conventional,Metric))

julia> inertia(English,Metric) # kg⋅slug⁻¹
$(inertia(English,Metric))

julia> inertia(IAU,Metric) # kg⋅m⊙⁻¹
$(inertia(IAU,Metric))

julia> inertia(PlanckGauss,Metric) # kg⋅mP⁻¹
$(inertia(PlanckGauss,Metric))
```
""" inertia

@doc """
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

julia> 1/elementarycharge(SI2019) # J⋅eV⁻¹
$(inv(elementarycharge(SI2019)))
```
""" energy
#julia> inv(kilocalorie(SI2019) # J⋅kcalₜₕ⁻¹
#$(inv(kilocalorie(SI2019)))
#julia> 1/BTUJ # J⋅BTU⁻¹
#$(1/1055.036345118633)

@doc """
$(convertext(:power,"energy(U,S)/time(U,S))"))

Radiant flux or `power` or `energy` per `time` (W, J⋅s⁻¹, kg⋅m²⋅s⁻³), unit conversion factor.

```Julia
julia> power(CGS,Metric) # W⋅s⋅erg⁻¹
$(power(CGS,Metric))

julia> power(English,Metric) # W⋅s⋅ft⁻¹⋅lb⁻¹
$(power(English,Metric))
```
""" power

@doc """
$(convertext(:force,"inertia(U,S)*acceleration(U,S)"))

Weight or force or `inertia` times `acceleration` (N, kg⋅m⋅s⁻²), unit conversion factor.

```Julia
julia> force(CGS,Metric) # N⋅dyn⁻¹
$(force(CGS,Metric))

julia> force(CGS,English) # lb⋅dyn⁻¹
$(force(CGS,English))

julia> force(English,Metric) # N⋅lb⁻¹
$(force(English,Metric))

julia> force(FPS,Metric) # pdl⋅N⁻¹
$(force(FPS,Metric))

julia> force(MetricEngineering,Metric) # kp⋅N⁻¹
$(force(MetricEngineering,Metric))
```
""" force, F

@doc """
$(convertext(:specificforce,"acceleration(U,S)/gravity(U,S)"))

Weight or `force` per `mass` or `gforce` (N/kg, m⋅s⁻²), unit conversion factor.
```Julia
julia> specificforce(CGS,Metric)
$(specificforce(CGS,Metric))

julia> specificforce(MetricEngineering,Metric)
$(specificforce(MetricEngineering,Metric))

julia> specificforce(English,Metric)
$(specificforce(English,Metric))
```
""" specificforce

@doc """
$(convertext(:gravityforce,"acceleration(U,S)/specificforce(U,S)"))

Reference `acceleration` per `specificforce` (𝟏, F⁻¹MLT⁻²), unit conversion factor.
```Julia
julia> gravityforce(Metric,CGS)
$(gravityforce(Metric,CGS))

julia> gravityforce(Metric,MetricEngineering)
$(gravityforce(Metric,MetricEngineering))

julia> gravityforce(Metric,English)
$(gravityforce(Metric,English))
```
""" gravityforce

@doc """
$(convertext(:pressure,"force(U,S)/area(U,S)"))

Pressure or stress or `force` per `area` (Pa, N⋅m⁻², kg⋅m⁻¹⋅s⁻²), unit conversion factor.

```Julia
julia> pressure(CGS,Metric) # Pa⋅Ba⁻¹
$(pressure(CGS,Metric))

julia> 1/atm # Pa⋅atm⁻¹
$(inv(atm))

julia> pressure(English,Metric) # Pa⋅ft²⋅lb⁻¹
$(pressure(English,Metric))

julia> pressure(Metric,IPS) # psi⋅Pa⁻¹
$(pressure(Metric,IPS))
```
""" pressure

# mechanical

@doc """
$(convertext(:impulse,"force(U,S)*time(U,S)"))

Linear `impulse` or `force` times `time` (N⋅s, kg⋅m⋅s⁻¹), unit conversion factor.

```Julia
julia> impulse(CGS,Metric) # N⋅dyn⁻¹
$(impulse(CGS,Metric))

julia> impulse(CGS,English) # lb⋅dyn⁻¹
$(impulse(CGS,English))

julia> impulse(English,Metric) # N⋅lb⁻¹
$(impulse(English,Metric))
```
""" impulse

@doc """
$(convertext(:momentum,"mass(U,S)*speed(U,S)"))

Linear `momentum` or `mass` times `speed` (N⋅s, kg⋅m⋅s⁻¹), unit conversion factor.

```Julia
julia> momentum(CGS,Metric) # N⋅dyn⁻¹
$(momentum(CGS,Metric))

julia> momentum(CGS,English) # lb⋅dyn⁻¹
$(momentum(CGS,English))

julia> momentum(British,Metric) # N⋅lb⁻¹
$(momentum(British,Metric))
```
""" momentum

@doc """
$(convertext(:angularmomentum,"momentum(U,S)*lengt(U,S)*angle(U,S)"))

Rotational momentum or `angularmomentum` (N⋅m⋅s, kg⋅m²⋅s⁻¹), unit conversion factor.

```Julia
julia> momentum(CGS,Metric) # N⋅m⋅dyn⁻¹⋅cm⁻¹
$(momentum(CGS,Metric))

julia> momentum(CGS,English) # lb⋅ft⋅dyn⁻¹⋅cm⁻¹
$(momentum(CGS,English))

julia> momentum(British,Metric) # N⋅m⋅lb⁻¹⋅ft⁻¹
$(momentum(British,Metric))
```
""" angularmomentum

@doc """
$(convertext(:yank,"mass(U,S)*jerk(U,S)"))

Rate of change of `force` or `yank` (N⋅s⁻¹, kg⋅m⋅s⁻³), unit conversion factor.

```Julia
julia> yank(CGS,Metric) # N⋅dyn⁻¹
$(yank(CGS,Metric))

julia> yank(CGS,English) # lb⋅dyn⁻¹
$(yank(CGS,English))

julia> yank(British,Metric) # N⋅lb⁻¹⋅
$(yank(British,Metric))
```
""" yank

@doc """
$(convertext(:areadensity,"mass(U,S)/area(U,S)"))

Surface or `areadensity` or `mass` per `area` (kg⋅m⁻²), unit conversion factor.

```Julia
julia> areadensity(CGS,Metric) # kg⋅cm²⋅g⁻¹⋅m⁻²
$(areadensity(CGS,Metric))

julia> areadensity(CGS,English) # lb⋅cm²⋅g⁻¹⋅ft⁻²
$(areadensity(CGS,English))

julia> areadensity(English,Metric) # kg⋅ft²⋅lb⁻¹⋅m⁻²
$(areadensity(English,Metric))

julia> areadensity(British,Metric) # kg⋅ft²⋅slug⁻¹⋅m⁻²
$(areadensity(British,Metric))
```
""" areadensity

@doc """
$(convertext(:density,"mass(U,S)/volume(U,S)"))

Specific mass or `mass` per `volume` or `density` (kg⋅m⁻³), unit conversion factor.

```Julia
julia> density(CGS,Metric) # kg⋅cm³⋅g⁻¹⋅m⁻³
$(density(CGS,Metric))

julia> density(CGS,Brtish) # slug⋅cm³⋅g⁻¹⋅ft⁻³
$(density(CGS,British))

julia> density(English,Metric) # kg⋅ft³⋅lb⁻¹⋅m⁻³
$(density(English,Metric))
```
""" density

@doc """
$(convertext(:specificweight,"force(U,S)/volume(U,S)"))

Specific weight or `force` per `volume` (N⋅m⁻³ or lb⋅ft⁻³), unit conversion factor.
```Julia
julia> specificweight(CGS,Metric) # N⋅cm³⋅dyn⁻¹⋅m⁻³
$(specificweight(CGS,Metric))

julia> specificweight(CGS,Brtish) # lb⋅cm³⋅dyn⁻¹⋅ft⁻³
$(specificweight(CGS,British))

julia> specificweight(English,Metric) # N⋅ft³⋅lb⁻¹⋅m⁻³
$(specificweight(English,Metric))
```
""" specificweight

@doc """
$(convertext(:specificvolume,"volume(U,S)/mass(U,S)"))

Reciprocal `density` or `volume` per `mass` or `specificvolume` (m³⋅kg), unit conversion factor.

```Julia
julia> specificvolume(CGS,Metric) # g⋅m³⋅kg⁻¹⋅cm⁻³
$(specificvolume(CGS,Metric))

julia> specificvolume(CGS,British) # kg⋅ft³⋅slug⁻¹⋅cm⁻³
$(specificvolume(CGS,British))

julia> specificvolume(English,Metric) # lb⋅m³⋅kg⁻¹⋅ft⁻³
$(specificvolume(English,Metric))
```
""" specificvolume

@doc """
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
""" action

#=@doc """
$(convertext(:stiffness,"force(U,S)/length(U,S)"))

Amount of `force` per `length` or `stiffness` (N⋅m⁻¹, J⋅m⁻², kg⋅s⁻²), unit conversion factor.

```Julia
julia> stiffness(CGS,Metric) # kg⋅g⁻¹
$(stiffness(CGS,Metric))

julia> stiffness(CGS,English) # lb⋅g⁻¹
$(stiffness(CGS,English))

julia> stiffness(English,Metric) # kg⋅lb⁻¹
$(stiffness(English,Metric))
```
""" stiffness=#

@doc """
$(convertext(:irradiance,"power(U,S)/area(U,S)"))

Heat flux density or irradiance or `power` per `area` (W⋅m⁻², kg⋅s⁻³), unit conversion factor.

```Julia
julia> irradiance(CGS,Metric) # kg⋅g⁻¹
$(irradiance(CGS,Metric))

julia> irradiance(CGS,English) # lb⋅g⁻¹
$(irradiance(CGS,English))

julia> irradiance(English,Metric) # kg⋅lb⁻¹
$(irradiance(English,Metric))
```
""" irradiance, intensity

@doc """
$(convertext(:radiance,"irradiance(U,S)/solidangle(U,S)"))

Radiance or `irradiance` per `solidangle` (W⋅m⁻²⋅sr⁻¹, kg⋅s⁻³⋅sr⁻¹), unit conversion factor.

```Julia
julia> radiance(CGS,Metric) # kg⋅g⁻¹
$(radiance(CGS,Metric))

julia> radiance(CGS,English) # lb⋅g⁻¹
$(radiance(CGS,English))

julia> radiance(English,Metric) # kg⋅lb⁻¹
$(radiance(English,Metric))
```
""" radiance

@doc """
$(convertext(:radiantintensity,"power(U,S)/solidangle(U,S)"))

Radiant intensity or `power` per `solidangle` (W⋅sr⁻¹, W⋅rad⁻²), unit conversion factor.

```Julia
julia> radiantintensity(CGS,Metric) # W⋅s⋅erg⁻¹
$(radiantintensity(CGS,Metric))

julia> radiantintensity(English,Metric) # W⋅s⋅ft⁻¹⋅lb⁻¹
$(radiantintensity(English,Metric))
```
""" radiantintensity

@doc """
$(convertext(:spectralexposure,"force(U,S)/speed(U,S)"))

Spectral exposure or `fluence` per `frequency` (N⋅s⋅m⁻¹, J⋅s⋅m⁻²), unit conversion factor.

```Julia
julia> spectralexposure(CGS,Metric) # kg⋅g⁻¹
$(spectralexposure(CGS,Metric))

julia> spectralexposure(CGS,English) # lb⋅g⁻¹
$(spectralexposure(CGS,English))

julia> spectralexposure(CODATA,Metric) # kg⋅kg⁻¹
$(spectralexposure(CODATA,Metric))

julia> spectralexposure(Conventional,Metric) # kg⋅kg⁻¹
$(spectralexposure(Conventional,Metric))

julia> spectralexposure(English,Metric) # kg⋅lb⁻¹
$(spectralexposure(English,Metric))
```
""" spectralexposure

@doc """
$(convertext(:diffusivity,"speed(U,S)*length(U,S)"))

Thermal `diffusivity` or kinematic viscostiy (m²⋅s⁻¹), unit conversion factor.

```Julia
julia> diffusivity(CGS,Metric) # m²⋅cm⁻²
$(diffusivity(CGS,Metric))

julia> diffusivity(English,Metric) # m²⋅ft⁻²
$(diffusivity(English,Metric))

julia> diffusivity(Survey,English) # ft²⋅ftUS⁻²
$(diffusivity(Survey,English))
```
""" diffusivity

@doc """
$(convertext(:viscosity,"inertia(U,S)/length(U,S)/time(U,S)"))

Resistance to deformation or `viscosity` (Pa⋅s, kg⋅m⁻¹⋅s⁻¹), unit conversion factor.

```Julia
julia> viscosity(CGS,Metric) # Pa⋅Ba⁻¹
$(viscosity(CGS,Metric))

julia> viscosity(English,Metric) # Pa⋅ft²⋅lb⁻¹
$(viscosity(English,Metric))

julia> viscosity(British,Metric) # Pa⋅ft²⋅lb⁻¹
$(viscosity(British,Metric))
```
""" viscosity

@doc """
$(convertext(:lineardensity,"mass(U,S)/length(U,S)"))

Amount of `lineardensity` or `mass` per `length` (kg⋅m⁻¹), unit conversion factor.

```Julia
julia> lineardensity(CGS,Metric) # kg⋅cm¹⋅g⁻¹⋅m⁻¹
$(lineardensity(CGS,Metric))

julia> lineardensity(CGS,British) # slug⋅cm¹⋅g⁻¹⋅ft⁻¹
$(lineardensity(CGS,British))

julia> lineardensity(English,Metric) # kg⋅ft¹⋅lb⁻¹⋅m⁻¹
$(lineardensity(English,Metric))
```
""" lineardensity

@doc """
$(convertext(:massflow,"mass(U,S)/time(U,S)"))

Rate of `massflow` or `mass` per `time` (kg⋅s⁻¹), unit conversion factor.

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
""" massflow

@doc """
$(convertext(:spectralflux,"power(U,S)/length(U,S)"))

Spectral power or `power` per wave `length` (W⋅m⁻¹), unit conversion factor.

```Julia
julia> spectralflux(CGS,Metric) # kg⋅m⋅g⁻¹⋅cm⁻¹
$(spectralflux(CGS,Metric))

julia> spectralflux(CGS,English) # lb⋅ft⋅g⁻¹⋅cm⁻¹
$(spectralflux(CGS,English))

julia> spectralflux(English,Metric) # kg⋅m⋅lb⁻¹⋅ft⁻¹
$(spectralflux(English,Metric))
```
""" spectralflux

@doc """
$(convertext(:powerdensity,"power(U,S)/volume(U,S)"))

Spectral irradiance (volume) or `powerdensity` (W⋅m⁻³), unit conversion factor.

```Julia
julia> powerdensity(CGS,Metric) # kg⋅cm⋅g⁻¹⋅m⁻¹
$(powerdensity(CGS,Metric))

julia> powerdensity(CGS,English) # lb⋅cm⋅g⁻¹⋅ft⁻¹
$(powerdensity(CGS,English))

julia> powerdensity(English,Metric) # kg⋅ft⋅lb⁻¹⋅m⁻¹
$(powerdensity(English,Metric))
```
""" powerdensity

@doc """
$(convertext(:compressibility,"1/pressure(U,S)"))

Relative volume change or `compressibility` (Pa⁻¹), unit conversion factor.

```Julia
julia> compressibility(CGS,Metric) # Ba⋅Pa⁻¹
$(compressibility(CGS,Metric))

julia> compressibility(English,Metric) # lb⋅ft⁻²⋅Pa⁻¹
$(compressibility(English,Metric))

julia> compressibility(Metric,IPS) # Pa⋅psi⁻¹
$(compressibility(Metric,IPS))
```
""" compressibility

@doc """
$(convertext(:fluence,"energy(U,S)/area(U,S"))

Radiant exposure or `force` per `length` or stiffness (N⋅m⁻¹, J⋅m⁻²), unit conversion factor.

```Julia
julia> fluence(CGS,Metric) # kg⋅g⁻¹
$(mass(CGS,Metric))

julia> fluence(CGS,English) # lb⋅g⁻¹
$(fluence(CGS,English))

julia> fluence(CODATA,Metric) # kg⋅kg⁻¹
$(mass(CODATA,Metric))

julia> fluence(Conventional,Metric) # kg⋅kg⁻¹
$(fluence(Conventional,Metric))

julia> fluence(English,Metric) # kg⋅lb⁻¹
$(fluence(English,Metric))
```
""" fluence

@doc """
$(convertext(:rotationalinertia,"mass(U,S)*area(U,S)"))

Moment of inertia or `rotationalinertia` (kg⋅m²), unit conversion factor.

```Julia
julia> rotationalinertia(CGS,Metric) # kg⋅m²⋅g⁻¹⋅cm⁻²
$(rotationalinertia(CGS,Metric))

julia> rotationalinertia(CGS,British) # slug⋅ft²⋅g⁻¹⋅cm⁻²
$(rotationalinertia(CGS,British))

julia> rotationalinertia(English,Metric) # kg⋅m²⋅lb⁻¹⋅ft⁻²
$(rotationalinertia(English,Metric))
```
""" rotationalinertia

# acoustic

@doc """
$(convertext(:soundexposure,"pressure(U,S)^2*time(U,S)"))

Square of `pressure` by `time` or `soundexposure` (Pa²⋅s, N²⋅m⁻⁴), unit conversion factor.

```Julia
julia> soundexposure(CGS,Metric) # Pa²⋅Ba⁻²
$(soundexposure(CGS,Metric))

julia> soundexposure(English,Metric) # Pa²⋅ft⁴⋅lb⁻²
$(soundexposure(English,Metric))
```
""" soundexposure

@doc """
$(convertext(:specificimpedance,"pressure(U,S)/speed(U,S)"))

Characteristic specific acoustic impedance (Rayl, Pa⋅s⋅m⁻¹), unit conversion factor.

```Julia
julia> specificimpedance(CGS,Metric) # Pa⋅cm⋅m⁻¹⋅Ba⁻¹
$(specificimpedance(CGS,Metric))

julia> specificimpedance(English,Metric) # Pa⋅ft³⋅m⁻¹⋅lb⁻¹
$(specificimpedance(English,Metric))
```
""" specificimpedance

@doc """
$(convertext(:impedance,"specificimpedance(U,S)/area(U,S)"))

Acoustic `impedance` (Rayl⋅m⁻², Pa⋅s⋅m⁻³, kg⋅s⁻¹⋅m⁻⁴), unit conversion factor.

```Julia
julia> impedance(CGS,Metric) # Pa⋅cm³⋅m⁻³⋅Ba⁻¹
$(impedance(CGS,Metric))

julia> impedance(English,Metric) # Pa⋅ft⁵⋅m⁻³⋅lb⁻¹
$(impedance(English,Metric))
```
""" impedance(U::UnitSystem,S::UnitSystem)

@doc """
$(convertext(:admittance,"area(U,S)/specificimpedance(U,S)"))

Acoustic `admittance` (m²⋅Rayl⁻¹, m³⋅s⁻¹⋅Pa⁻¹, m⁴⋅s⋅kg⁻¹), unit conversion factor.

```Julia
julia> admittance(CGS,Metric) # Ba⋅m³⋅cm⁻³⋅Pa⁻¹
$(admittance(CGS,Metric))

julia> admittance(English,Metric) # lb⋅m³⋅ft⁻⁵⋅Pa⁻¹
$(admittance(English,Metric))
```
""" admittance

@doc """
$(convertext(:compliance,"time(U,S)^2/mass(U,S)"))

Acoustic `compliance` is reciprocal of `fluence` (m⋅N⁻¹, m³⋅Pa⁻¹), unit conversion factor.

```Julia
julia> compliance(CGS,Metric) # kg⋅g⁻¹
$(compliance(CGS,Metric))

julia> compliance(CGS,English) # slug⋅g⁻¹
$(compliance(CGS,English))

julia> compliance(English,Metric) # kg⋅lb⁻¹
$(compliance(English,Metric))
```
""" compliance

@doc """
$(convertext(:inertance,"mass(U,S)/length(U,S)^4"))

Acoustic mass or `inertance` (kg⋅m⁴, Pa⋅s²⋅m⁻³), unit conversion factor.

```Julia
julia> inertance(CGS,Metric) # kg⋅cm⁴⋅g⁻¹⋅m⁻⁴
$(inertance(CGS,Metric))

julia> inertance(CGS,English) # slug⋅cm⁴⋅g⁻¹⋅ft⁻⁴
$(inertance(CGS,English))

julia> inertance(English,Metric) # kg⋅ft⁴⋅lb⁻¹⋅m⁻⁴
$(inertance(English,Metric))
```
""" inertance
