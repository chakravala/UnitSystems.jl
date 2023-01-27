
#   This file is part of UnitSystems.jl
#   It is licensed under the MIT license
#   UnitSystems Copyright (C) 2021 Michael Reed
#       _           _                         _
#      | |         | |                       | |
#   ___| |__   __ _| | ___ __ __ ___   ____ _| | __ _
#  / __| '_ \ / _` | |/ / '__/ _` \ \ / / _` | |/ _` |
# | (__| | | | (_| |   <| | | (_| |\ V / (_| | | (_| |
#  \___|_| |_|\__,_|_|\_\_|  \__,_| \_/ \__,_|_|\__,_|
#
#   https://github.com/chakravala
#   https://crucialflow.com

export deka,hecto,kilo,mega,giga,tera,peta,exa,zetta,yotta,sixty
export deci,centi,milli,micro,nano,pico,femto,atto,zepto,yocto
export byte,kibi,mebi,gibi,tebi,pebi,exbi,zebi,yobi

export slug, ft, KJ1990, KJ2014, RK1990, RK2014, mₑ1990, mₑ2014, temp, units, °R, T₀, eV
export slugs, kilograms, lbm, meters, feet, rankine, kelvin, moles, molecules, universal
export Universe, UnitSystem, US, universe, HOUR, DAY, th, lc, mc, tcq, lcq, mcq, ς, φ, 𝟔𝟎
export similitude, 𝟙, F, M, L, T, Q, Θ, N, J, A, Λ, C, sackurtetrode,atomicmass,intensity
export °R, τ, 𝟏𝟎, 𝟐, 𝟑, 𝟓, nm, 𝟏, mₑ, μ₀, Mᵤ, Rᵤ, αG, GG, slug, ħ, μₚₑ, αL, 𝟕, 𝟏𝟏, 𝟏𝟗, 𝟒𝟑
export two, three, five, seven, eleven, nineteen, fourtythree, tau, eulergamma, golden

const EMU2019,ESU2019,stiffness,intensity,atomicmass = EMU,ESU,fluence,irradiance,dalton

# == Metric is different
const eV = electronvolt(SI2019)
const κ = einstein(SI2019)
const σ = stefan(SI2019) #
const μB = magneton(SI2019) #
const ε₀ = vacuumpermittivity(SI2019) #
const kₑ = electrostatic(SI2019) #
const mₚ = protonmass(SI2019)
const Da = dalton(SI2019)
const 𝔉 = faraday(SI2019) #
const Φ₀ = magneticfluxquantum(SI2019) #
const Z₀ = vacuumimpedance(SI2019) #
const G₀ = conductancequantum(SI2019) #
const Eₕ = hartree(SI2019)
const a₀ = bohr(SI2019)
const rₑ = electronradius(SI2019)
const RH,Ry = R∞*mₚ/(electronmass(SI2019)+mₚ),𝘩*𝘤*R∞

const ℓP = length(PlanckGauss,SI2019)
const tP = time(PlanckGauss,SI2019)
const TP = temperature(PlanckGauss,SI2019)

const lS = length(Stoney,SI2019)
const tS = time(Stoney,SI2019)
const mS = mass(Stoney,SI2019)
const qS = charge(Stoney,SI2019)

const lA = length(Hartree,SI2019)
const tA = time(Hartree,SI2019)
const mA = mass(Hartree,SI2019)
const qA = charge(Hartree,SI2019)

const lQCD = length(QCD,SI2019)
const tQCD = time(QCD,SI2019)
const mQCD = mass(QCD,SI2019)

# non standard units

const BTU = thermalunit(British)
const BTUftlb = thermalunit(British) # BTU⋅ft⁻¹⋅lb⁻¹
const BTUJ = thermalunit(SI2019) # BTU⋅J⁻¹
const HP = horsepower(Metric)
const gal = gallon(Metric)
const kcal = kilocalorie(SI2019)
const cal = calorie(SI2019)
const universal = molargas

# constant aliases

const US,mpe, mep, meu, mpu, ainv, aG = UnitSystem,μₚₑ, μₑₚ, μₑᵤ, μₚᵤ, αinv, αG
const Mu,Ru,SB,hh,cc,m0,e0,ke,me,mp,mu,mᵤ,ee,FF,Z0,G0,Eh,a0,re,g0,lP,aL,ϵ₀ = Mᵤ,Rᵤ,σ,𝘩,𝘤,μ₀,ε₀,kₑ,mₑ,mₚ,Da,Da,𝘦,𝔉,Z₀,G₀,Eₕ,a₀,rₑ,g₀,ℓP,αL,ε₀
export κ, G, GG, NA, kB, Rᵤ, σ, 𝘩, ħ, 𝘤, μ₀, ε₀, kₑ, mₑ, mₚ, Da, 𝘦, 𝔉, Φ₀, Z₀, G₀, Eₕ, R∞, a₀, rₑ, KJ, RK, Ru, SB, hh, cc, m0, e0, ke, me, mp, mu, mᵤ, ee, FF, Z0, G0, Eh, a0, re, μB
export αG, αinv, μₚₑ, μₑᵤ, μₚᵤ, mpe, meu, mpu, mP, δμ₀, Mᵤ, Mu, RH, Ry, ΔνCs, Kcd, ainv
export cal, kcal, calₜₕ, kcalₜₕ, calᵢₜ, kcalᵢₜ, ℓP, g₀, g0, atm, lbm, BTUJ, BTUftlb, aG
export lP, tP, TP, lS, tS, mS, qS, lA, tA, mA, qA, lQCD, tQCD, mQCD, ϵ₀, αL, aL, μₑₚ, 𝘦ₙ
export GM☉, GME, GMJ, LD, JD, lb, fur, ftUS, Rᵤ2014, Ωᵢₜ, Vᵢₜ, em, mi, 𝘦ᵣ

# engineering unit systems docs

@doc """
$(systext(:Metric,"MetricSystem(milli,𝟐*τ/𝟏𝟎^7)"))

Standard `Metric` system based on exact `molarmass` and `vacuumpermeability`.

```Julia
julia> boltzmann(Metric) # J⋅K⁻¹
$(boltzmann(Metric))

julia> planckreduced(Metric) # J⋅s⋅rad⁻¹
$(planckreduced(Metric))

julia> lightspeed(Metric) # m⋅s⁻¹
$(lightspeed(Metric))

julia> vacuumpermeability(Metric) # H⋅m⁻¹
$(vacuumpermeability(Metric))

julia> electronmass(Metric) # kg
$(electronmass(Metric))

julia> molarmass(Metric) # kg⋅mol⁻¹
$(molarmass(Metric))

julia> luminousefficacy(Metric) # lm⋅W⁻¹
$(luminousefficacy(Metric))
```
""" Metric, MKS

@doc """
$(systext(:SI2019,"MetricSystem(Mᵤ,μ₀)"))

Systeme International d'Unites based on approximate `molarmass` and `vacuumpermeability`.

```Julia
julia> boltzmann(SI2019) # J⋅K⁻¹
$(boltzmann(SI2019))

julia> planckreduced(SI2019) # J⋅s⋅rad⁻¹
$(planckreduced(SI2019))

julia> lightspeed(SI2019) # m⋅s⁻¹
$(lightspeed(SI2019))

julia> vacuumpermeability(SI2019) # H⋅m⁻¹
$(vacuumpermeability(SI2019))

julia> electronmass(SI2019) # kg
$(electronmass(SI2019))

julia> molarmass(SI2019) # kg⋅mol⁻¹
$(molarmass(SI2019))

julia> luminousefficacy(SI2019) # lm⋅W⁻¹
$(luminousefficacy(SI2019))
```
""" SI2019, SI

@doc """
$(systext(:MetricTurn,"MetricSystem(milli,𝟐*τ/𝟏𝟎^7,Rᵤ,𝟏,𝟏/τ)"))

Standard `MetricTurn` system based on exact `molarmass` and `vacuumpermeability`.

```Julia
julia> boltzmann(MetricTurn) # J⋅K⁻¹
$(boltzmann(MetricTurn))

julia> planckreduced(MetricTurn) # J⋅s⋅τ⁻¹
$(planckreduced(MetricTurn))

julia> lightspeed(MetricTurn) # m⋅s⁻¹
$(lightspeed(MetricTurn))

julia> vacuumpermeability(MetricTurn) # H⋅m⁻¹
$(vacuumpermeability(MetricTurn))

julia> electronmass(MetricTurn) # kg
$(electronmass(MetricTurn))

julia> molarmass(MetricTurn) # kg⋅mol⁻¹
$(molarmass(MetricTurn))

julia> luminousefficacy(MetricTurn) # lm⋅W⁻¹
$(luminousefficacy(MetricTurn))
```
""" MetricTurn

@doc """
$(systext(:MetricDegree,"MetricSystem(milli,𝟐*τ/𝟏𝟎^7,Rᵤ,𝟏,𝟐^3*𝟑^2*𝟓/τ)"))

Standard `MetricDegree` system based on exact `molarmass` and `vacuumpermeability`.

```Julia
julia> boltzmann(MetricDegree) # J⋅K⁻¹
$(boltzmann(MetricDegree))

julia> planckreduced(MetricDegree) # J⋅s⋅deg⁻¹
$(planckreduced(MetricDegree))

julia> lightspeed(MetricDegree) # m⋅s⁻¹
$(lightspeed(MetricDegree))

julia> vacuumpermeability(MetricDegree) # H⋅m⁻¹
$(vacuumpermeability(MetricDegree))

julia> electronmass(MetricDegree) # kg
$(electronmass(MetricDegree))

julia> molarmass(MetricDegree) # kg⋅mol⁻¹
$(molarmass(MetricDegree))

julia> luminousefficacy(MetricDegree) # lm⋅W⁻¹
$(luminousefficacy(MetricDegree))
```
""" MetricDegree

@doc """
$(systext(:MetricArcminute,"MetricSystem(milli,𝟐*τ/𝟏𝟎^7,Rᵤ,𝟏,𝟐^5*𝟑^3*𝟓^2/τ)"))

Standard `MetricArcminute` system based on exact `molarmass` and `vacuumpermeability`.

```Julia
julia> boltzmann(MetricArcminute) # J⋅K⁻¹
$(boltzmann(MetricArcminute))

julia> planckreduced(MetricArcminute) # J⋅s⋅amin⁻¹
$(planckreduced(MetricArcminute))

julia> lightspeed(MetricArcminute) # m⋅s⁻¹
$(lightspeed(MetricArcminute))

julia> vacuumpermeability(MetricArcminute) # H⋅m⁻¹
$(vacuumpermeability(MetricArcminute))

julia> electronmass(MetricArcminute) # kg
$(electronmass(MetricArcminute))

julia> molarmass(MetricArcminute) # kg⋅mol⁻¹
$(molarmass(MetricArcminute))

julia> luminousefficacy(MetricArcminute) # lm⋅W⁻¹
$(luminousefficacy(MetricArcminute))
```
""" MetricArcminute

@doc """
$(systext(:MetricArcsecond,"MetricSystem(milli,𝟐*τ/𝟏𝟎^7,Rᵤ,𝟏,𝟐^7*𝟑^4*3/τ)"))

Standard `MetricArcsecond` system based on exact `molarmass` and `vacuumpermeability`.

```Julia
julia> boltzmann(MetricArcsecond) # J⋅K⁻¹
$(boltzmann(MetricArcsecond))

julia> planckreduced(MetricArcsecond) # J⋅s⋅asec⁻¹
$(planckreduced(MetricArcsecond))

julia> lightspeed(MetricArcsecond) # m⋅s⁻¹
$(lightspeed(MetricArcsecond))

julia> vacuumpermeability(MetricArcsecond) # H⋅m⁻¹
$(vacuumpermeability(MetricArcsecond))

julia> electronmass(MetricArcsecond) # kg
$(electronmass(MetricArcsecond))

julia> molarmass(MetricArcsecond) # kg⋅mol⁻¹
$(molarmass(MetricArcsecond))

julia> luminousefficacy(MetricArcsecond) # lm⋅W⁻¹
$(luminousefficacy(MetricArcsecond))
```
""" MetricArcsecond

@doc """
$(systext(:MetricGradian,"MetricSystem(milli,𝟐*τ/𝟏𝟎^7,Rᵤ,𝟏,𝟐^4*𝟓^2/τ)"))

Standard `MetricGradian` system based on exact `molarmass` and `vacuumpermeability`.

```Julia
julia> boltzmann(MetricGradian) # J⋅K⁻¹
$(boltzmann(MetricGradian))

julia> planckreduced(MetricGradian) # J⋅s⋅gon⁻¹
$(planckreduced(MetricGradian))

julia> lightspeed(MetricGradian) # m⋅s⁻¹
$(lightspeed(MetricGradian))

julia> vacuumpermeability(MetricGradian) # H⋅m⁻¹
$(vacuumpermeability(MetricGradian))

julia> electronmass(MetricGradian) # kg
$(electronmass(MetricGradian))

julia> molarmass(MetricGradian) # kg⋅mol⁻¹
$(molarmass(MetricGradian))

julia> luminousefficacy(MetricGradian) # lm⋅W⁻¹
$(luminousefficacy(MetricGradian))
```
""" MetricGradian

@doc """
$(systext(:Engineering,"MetricSystem(milli,𝟐*τ/𝟏𝟎^7,Rᵤ,g₀)"))

Standard Metric `Engineering` system based on kilogram and kilopond (kilogram-force) units.

```Julia
julia> boltzmann(Engineering) # kgf⋅m⋅K⁻¹
$(boltzmann(Engineering))

julia> planckreduced(Engineering) # kgf⋅m⋅s⋅rad⁻¹
$(planckreduced(Engineering))

julia> lightspeed(Engineering) # m⋅s⁻¹
$(lightspeed(Engineering))

julia> vacuumpermeability(Engineering) # kgf⋅s²⋅C⁻²
$(vacuumpermeability(Engineering))

julia> electronmass(Engineering) # kg
$(electronmass(Engineering))

julia> molarmass(Engineering) # kg⋅mol⁻¹
$(molarmass(Engineering))

julia> luminousefficacy(Engineering) # lm⋅s⋅m⁻¹⋅kgf⁻¹
$(luminousefficacy(Engineering))

julia> gravity(Engineering) # kg⋅m⋅kgf⁻¹⋅s⁻²
$(gravity(Engineering))
```
""" Engineering, ME

#=@doc """
$(systext(:SI2019Engineering,"MetricSystem(Mᵤ,μ₀,Rᵤ,g₀)"))

Systeme International d'Unites based on kilogram and kilopond (kilogram-force) units.

```Julia
julia> boltzmann(SI2019Engineering) # kgf⋅m⋅K⁻¹
$(boltzmann(SI2019Engineering))

julia> planckreduced(SI2019Engineering) # kgf⋅m⋅s⋅rad⁻¹
$(planckreduced(SI2019Engineering))

julia> lightspeed(SI2019Engineering) # m⋅s⁻¹
$(lightspeed(SI2019Engineering))

julia> vacuumpermeability(SI2019Engineering) # kgf⋅s²⋅C⁻²
$(vacuumpermeability(SI2019Engineering))

julia> electronmass(SI2019Engineering) # kg
$(electronmass(SI2019Engineering))

julia> molarmass(SI2019Engineering) # kg⋅mol⁻¹
$(molarmass(SI2019Engineering))

julia> luminousefficacy(SI2019Engineering) # lm⋅s⋅m⁻¹⋅kgf⁻¹
$(luminousefficacy(SI2019Engineering))

julia> gravity(SI2019Engineering) # kg⋅m⋅kgf⁻¹⋅s⁻²
$(gravity(SI2019Engineering))
```
""" SI2019Engineering, SIE=#

@doc """
$(systext(:SI1976,"MetricSystem(milli,𝟐*τ/𝟏𝟎^7,8.31432)"))

Reference `UnitSystem` with universal gas constant of `8.31432` from 1976 standard atmosphere.

```Julia
julia> boltzmann(SI1976) # J⋅K⁻¹
$(boltzmann(SI1976))

julia> planckreduced(SI1976) # J⋅s⋅rad⁻¹
$(planckreduced(SI1976))

julia> lightspeed(SI1976) # m⋅s⁻¹
$(lightspeed(SI1976))

julia> vacuumpermeability(SI1976) # H⋅m⁻¹
$(vacuumpermeability(SI1976))

julia> electronmass(SI1976) # kg
$(electronmass(SI1976))

julia> molarmass(SI1976) # kg⋅mol⁻¹
$(molarmass(SI1976))

julia> luminousefficacy(SI1976) # lm⋅W⁻¹
$(luminousefficacy(SI1976))
```
""" SI1976

@doc """
$(systext(:CODATA,"ConventionalSystem(RK2014,KJ2014,Rᵤ2014)"))

Reference `UnitSystem` based on Committee on Data of the International Science Council.

```Julia
julia> josephson(CODATA) # Hz⋅V⁻¹
$(josephson(CODATA))

julia> klitzing(CODATA) # Ω
$(klitzing(CODATA))

julia> boltzmann(CODATA) # J⋅K⁻¹
$(boltzmann(CODATA))

julia> planckreduced(CODATA) # J⋅s⋅rad⁻¹
$(planckreduced(CODATA))

julia> lightspeed(CODATA) # m⋅s⁻¹
$(lightspeed(CODATA))

julia> vacuumpermeability(CODATA) # H⋅m⁻¹
$(vacuumpermeability(CODATA))

julia> electronmass(CODATA) # kg
$(electronmass(CODATA))

julia> molarmass(CODATA) # kg⋅mol⁻¹
$(molarmass(CODATA))

julia> luminousefficacy(CODATA) # lm⋅W⁻¹
$(luminousefficacy(CODATA))
```
""" CODATA

@doc """
$(systext(:Conventional,"ConventionalSystem(RK1990,KJ2014)"))

Conventional electronic `UnitSystem` with 1990 tuned `josephson` and `klitzing` constants.

```Julia
julia> josephson(Conventional) # Hz⋅V⁻¹
$(josephson(Conventional))

julia> klitzing(Conventional) # Ω
$(klitzing(Conventional))

julia> boltzmann(Conventional) # J⋅K⁻¹
$(boltzmann(Conventional))

julia> planckreduced(Conventional) # J⋅s⋅rad⁻¹
$(planckreduced(Conventional))

julia> lightspeed(Conventional) # m⋅s⁻¹
$(lightspeed(Conventional))

julia> vacuumpermeability(Conventional) # H⋅m⁻¹
$(vacuumpermeability(Conventional))

julia> electronmass(Conventional) # kg
$(electronmass(Conventional))

julia> molarmass(Conventional) # kg⋅mol⁻¹
$(molarmass(Conventional))

julia> luminousefficacy(Conventional) # lm⋅W⁻¹
$(luminousefficacy(Conventional))
```
""" Conventional

@doc """
$(systext(:International,"ElectricSystem(Metric,$Ωᵢₜ,$Vᵢₜ)"))

International `UnitSystem` with United States measurements of `Ωᵢₜ` and `Vᵢₜ`.

```Julia
julia> resistance(International,Metric) # Ω⋅Ω⁻¹
$(resistance(International,Metric))

julia> electricpotential(International,Metric) # V⋅V⁻¹
$(electricpotential(International,Metric))

julia> boltzmann(International) # J⋅K⁻¹
$(boltzmann(International))

julia> planckreduced(International) # J⋅s⋅rad⁻¹
$(planckreduced(International))

julia> lightspeed(International) # m⋅s⁻¹
$(lightspeed(International))

julia> vacuumpermeability(International) # H⋅m⁻¹
$(vacuumpermeability(International))

julia> electronmass(International) # kg
$(electronmass(International))

julia> molarmass(International) # kg⋅mol⁻¹
$(molarmass(International))

julia> luminousefficacy(International) # lm⋅W⁻¹
$(luminousefficacy(International))
```
""" International

@doc """
$(systext(:InternationalMean,"ElectricSystem(Metric,1.00049,1.00034)"))

International `UnitSystem` with mean measurements of `Ωᵢₜ` and `Vᵢₜ`.

```Julia
julia> resistance(InternationalMean,Metric) # Ω⋅Ω⁻¹
$(resistance(InternationalMean,Metric))

julia> electricpotential(InternationalMean,Metric) # V⋅V⁻¹
$(electricpotential(InternationalMean,Metric))

julia> boltzmann(InternationalMean) # J⋅K⁻¹
$(boltzmann(InternationalMean))

julia> planckreduced(InternationalMean) # J⋅s⋅rad⁻¹
$(planckreduced(InternationalMean))

julia> lightspeed(InternationalMean) # m⋅s⁻¹
$(lightspeed(InternationalMean))

julia> vacuumpermeability(InternationalMean) # H⋅m⁻¹
$(vacuumpermeability(InternationalMean))

julia> electronmass(InternationalMean) # kg
$(electronmass(InternationalMean))

julia> molarmass(InternationalMean) # kg⋅mol⁻¹
$(molarmass(InternationalMean))

julia> luminousefficacy(International) # lm⋅W⁻¹
$(luminousefficacy(International))
```
""" InternationalMean

@doc """
$(systext(:Meridian,"EntropySystem(Metric,𝟏,em,em^3,𝟏,τ/𝟐^6/𝟓^7,milli)"))

Modern ideal `Meridian` system defined by France's original `earthmeter` definition.

```Julia
julia> greatcircle(Meridian) # em
$(greatcircle(Meridian))

julia> boltzmann(Meridian) # eJ⋅K⁻¹
$(boltzmann(Meridian))

julia> planckreduced(Meridian) # eJ⋅s⋅rad⁻¹
$(planckreduced(Meridian))

julia> lightspeed(Meridian) # em⋅s⁻¹
$(lightspeed(Meridian))

julia> vacuumpermeability(Meridian) # kegf⋅s²⋅eC⁻²
$(vacuumpermeability(Meridian))

julia> electronmass(Meridian) # keg
$(electronmass(Meridian))

julia> molarmass(Meridian) # keg⋅eg-mol⁻¹
$(molarmass(Meridian))

julia> luminousefficacy(Meridian) # lm⋅W⁻¹
$(luminousefficacy(Meridian))
```
""" Meridian

#=@doc """
$(systext(:MeridianEngineering,"EntropySystem(Engineering,𝟏,em,em^3,𝟏,τ/𝟐^6/𝟓^7/g₀^2,milli)"))

Modern ideal engineering `UnitSystem` variant of the original French `Meridian` system.

```Julia
julia> greatcircle(MeridianEngineering) # em
$(greatcircle(MeridianEngineering))

julia> boltzmann(MeridianEngineering) # kegf⋅em⋅K⁻¹
$(boltzmann(MeridianEngineering))

julia> planckreduced(MeridianEngineering) # kegf⋅em⋅s⋅rad⁻¹
$(planckreduced(MeridianEngineering))

julia> lightspeed(MeridianEngineering) # em⋅s⁻¹
$(lightspeed(MeridianEngineering))

julia> vacuumpermeability(MeridianEngineering) # kegf⋅s²⋅eC⁻²
$(vacuumpermeability(MeridianEngineering))

julia> electronmass(MeridianEngineering) # keg
$(electronmass(MeridianEngineering))

julia> molarmass(MeridianEngineering) # keg⋅eg-mol⁻¹
$(molarmass(MeridianEngineering))

julia> luminousefficacy(MeridianEngineering) # lm⋅s⋅m⁻¹⋅kgf⁻¹
$(luminousefficacy(MeridianEngineering))

julia> gravity(MeridianEngineering) # keg⋅em⋅kegf⁻¹⋅s⁻²
$(gravity(MeridianEngineering))
```
""" MeridianEngineering=#


cgstext(US,AMP,cgs=eval(US)) = """
```Julia
julia> boltzmann($US) # erg⋅K⁻¹
$(boltzmann(cgs))

julia> planckreduced($US) # erg⋅s⋅rad⁻¹
$(planckreduced(cgs))

julia> lightspeed($US) # cm⋅s⁻¹
$(lightspeed(cgs))

julia> vacuumpermeability($US) # $(US ∈ (:CGSe,:ESU) ? "stat" : "ab")H⋅cm⁻¹
$(vacuumpermeability(cgs))

julia> electronmass($US) # g
$(electronmass(cgs))

julia> molarmass($US) # g⋅mol⁻¹
$(molarmass(cgs))

julia> luminousefficacy($US) # lm⋅s⋅erg⁻¹
$(luminousefficacy(cgs))

julia> rationalization($US)
$(rationalization(cgs))
```
"""

for U ∈ (:CGSm,:CGSe,:EMU,:ESU)
    (EU,AMP) = QuoteNode.(U ∉ (:CGSe,:ESU) ? (:EMU,:Bi) : (:ESU,:statA))
@eval @doc """
$(systext($(QuoteNode(U)),"GaussSystem(Metric,$($EU≠:EMU ? "(𝟏𝟎*𝘤)^-2" : "𝟏"),𝟐*τ)"))

Centimetre-gram-second `UnitSystem` variant based on `$($EU)` (non-rationalized).

$(cgstext($(QuoteNode(U)),$AMP))
""" $U

#=U ∉ (:CGSm,:CGSe) && @eval @doc """
    $(Symbol($(QuoteNode(U)),:2019)) = EntropySystem(SI2019,𝟏,0.01,0.001,𝟏,$($EU≠:EMU ? "1e3*μ₀/𝘤^2" : "1e7*μ₀"))

Centimetre-gram-second `UnitSystem` variant of tuned `SI2019` based on `$($EU)` (rationalized).

$(cgstext(Symbol($(QuoteNode(U)),:2019),$AMP))
""" $(Symbol(U,:2019))=#
end

#=@doc """
$(systext(:Thomson,"GaussSystem(Metric,𝟏,𝟐*τ,𝟏/𝟐)"))

Centimetre-gram-second `UnitSystem` variant `Thomson` (EMU-Lorentz, non-rationalized).

```Julia
julia> boltzmann(Thomson) # erg⋅K⁻¹
$(boltzmann(Thomson))

julia> planckreduced(Thomson) # erg⋅s⋅rad⁻¹
$(planckreduced(Thomson))

julia> lightspeed(Thomson) # cm⋅s⁻¹
$(lightspeed(Thomson))

julia> vacuumpermeability(Thomson) # abH⋅cm⁻¹
$(vacuumpermeability(Thomson))

julia> electronmass(Thomson) # g
$(electronmass(Thomson))

julia> molarmass(Thomson) # g⋅mol⁻¹
$(molarmass(Thomson))

julia> luminousefficacy(Thomson) # lm⋅s⋅erg⁻¹
$(luminousefficacy(Thomson))

julia> rationalization(Thomson)
$(rationalization(Thomson))

julia> lorentz(Thomson)
$(lorentz(Thomson))
```
""" Thomson=#

@doc """
$(systext(:Gauss,"GaussSystem(Metric,𝟏,𝟐*τ,𝟏𝟎^-2/𝘤)"))

Centimetre-gram-second `UnitSystem` variant `CGS` (Gauss-Lorentz, non-rationalized).

```Julia
julia> boltzmann(Gauss) # erg⋅K⁻¹
$(boltzmann(Gauss))

julia> planckreduced(Gauss) # erg⋅s⋅rad⁻¹
$(planckreduced(Gauss))

julia> lightspeed(Gauss) # cm⋅s⁻¹
$(lightspeed(Gauss))

julia> vacuumpermeability(Gauss) # statH⋅cm⁻¹
$(vacuumpermeability(Gauss))

julia> electronmass(Gauss) # g
$(electronmass(Gauss))

julia> molarmass(Gauss) # g⋅mol⁻¹
$(molarmass(Gauss))

julia> luminousefficacy(Gauss) # lm⋅s⋅erg⁻¹
$(luminousefficacy(Gauss))

julia> rationalization(Gauss)
$(rationalization(Gauss))

julia> lorentz(Gauss)
$(lorentz(Gauss))
```
""" Gauss, CGS

@doc """
$(systext(:LorentzHeaviside,"GaussSystem(Metric,𝟏,𝟏,centi/𝘤)"))

Centimetre-gram-second `UnitSystem` variant `HLU` (Heaviside-Lorentz, rationalized).

```Julia
julia> boltzmann(LorentzHeaviside) # erg⋅K⁻¹
$(boltzmann(LorentzHeaviside))

julia> planckreduced(LorentzHeaviside) # erg⋅s⋅rad⁻¹
$(planckreduced(LorentzHeaviside))

julia> lightspeed(LorentzHeaviside) # cm⋅s⁻¹
$(lightspeed(LorentzHeaviside))

julia> vacuumpermeability(HLU) # hlH⋅cm⁻¹
$(vacuumpermeability(LorentzHeaviside))

julia> electronmass(LorentzHeaviside) # g
$(electronmass(LorentzHeaviside))

julia> molarmass(LorentzHeaviside) # g⋅mol⁻¹
$(molarmass(LorentzHeaviside))

julia> luminousefficacy(LorentzHeaviside) # lm⋅s⋅erg⁻¹
$(luminousefficacy(LorentzHeaviside))

julia> rationalization(LorentzHeaviside)
$(rationalization(LorentzHeaviside))

julia> lorentz(LorentzHeaviside)
$(lorentz(LorentzHeaviside))
```
""" LorentzHeaviside, HLU

#=@doc """
$(systext(:Kennelly,"GaussSystem(Metric,𝟏𝟎^-7,𝟐*τ,𝟏,𝟏,𝟏)"))

Kennelly ? variant `UnitSystem` of the standard `Metric` units ???

```Julia
julia> boltzmann(Kennelly) # J⋅K⁻¹
$(boltzmann(Kennelly))

julia> planckreduced(Kennelly) # J⋅s⋅rad⁻¹
$(planckreduced(Kennelly))

julia> lightspeed(Kennelly) # m⋅s⁻¹
$(lightspeed(Kennelly))

julia> vacuumpermeability(Kennelly) # H⋅m⁻¹
$(vacuumpermeability(Kennelly))

julia> electronmass(Kennelly) # kg
$(electronmass(Kennelly))

julia> molarmass(Kennelly) # kg⋅mol⁻¹
$(molarmass(Kennelly))

julia> luminousefficacy(Kennelly) # lm⋅W⁻¹
$(luminousefficacy(Kennelly))

julia> rationalization(Kennelly)
$(rationalization(Kennelly))
```
""" Kennelly=#

@doc """
$(systext(:Gravitational,"EntropySystem(Metric,𝟏,𝟏,g₀)"))

Standard `Gravitational` system based on `hyl` and `kilopond` units.

```Julia
julia> boltzmann(Gravitational) # kgf⋅m⋅K⁻¹
$(boltzmann(Gravitational))

julia> planckreduced(Gravitational) # kgf⋅m⋅s⋅rad⁻¹
$(planckreduced(Gravitational))

julia> lightspeed(Gravitational) # m⋅s⁻¹
$(lightspeed(Gravitational))

julia> vacuumpermeability(Gravitational) # H⋅m⁻¹
$(vacuumpermeability(Gravitational))

julia> electronmass(Gravitational) # hyl
$(electronmass(Gravitational))

julia> molarmass(Gravitational) # hyl⋅mol⁻¹
$(molarmass(Gravitational))

julia> luminousefficacy(Gravitational) # lm⋅s⋅m⁻¹⋅kgf⁻¹
$(luminousefficacy(Gravitational))
```
""" Gravitational, GM

#=@doc """
$(systext(:GraviationalSI2019,"EntropySystem(SI2019,𝟏,𝟏,g₀)"))

Gravitational Systeme International d'Unites based on `hyl` and `kilopond` units.

```Julia
julia> boltzmann(GravitationalSI2019) # kgf⋅m⋅K⁻¹
$(boltzmann(GravitationalSI2019))

julia> planckreduced(GravitationalSI2019) # kgf⋅m⋅s⋅rad⁻¹
$(planckreduced(GravitationalSI2019))

julia> lightspeed(GravitationalSI2019) # m⋅s⁻¹
$(lightspeed(GravitationalSI2019))

julia> vacuumpermeability(GravitationalSI2019) # kgf⋅s²⋅C⁻²
$(vacuumpermeability(GravitationalSI2019))

julia> electronmass(GravitationalSI2019) # hyl
$(electronmass(GravitationalSI2019))

julia> molarmass(GravitationalSI2019) # hyl⋅mol⁻¹
$(molarmass(GravitationalSI2019))

julia> luminousefficacy(Gravitational) # lm⋅s⋅m⁻¹⋅kgf⁻¹
$(luminousefficacy(Gravitational))
```
""" GravitationalSI2019, GSI, GSI2019=#

#=@doc """
$(systext(:GravitationalMeridian,"EntropySystem(Metric,𝟏,em,g₀*em^3,𝟏,τ/𝟐^6/𝟓^7/g₀,milli)"))

Gravitational `UnitSystem` variant of the original French `Meridian` unit defintion.

```Julia
julia> greatcircle(GravitationalMeridian) # em
$(greatcircle(GravitationalMeridian))

julia> boltzmann(GravitationalMeridian) # kegf⋅em⋅K⁻¹
$(boltzmann(GravitationalMeridian))

julia> planckreduced(GravitationalMeridian) # kegf⋅em⋅s⋅rad⁻¹
$(planckreduced(GravitationalMeridian))

julia> lightspeed(GravitationalMeridian) # em⋅s⁻¹
$(lightspeed(GravitationalMeridian))

julia> vacuumpermeability(GravitationalMeridian) # kegf⋅s²⋅eC⁻²
$(vacuumpermeability(GravitationalMeridian))

julia> electronmass(GravitationalMeridian) # ehyl
$(electronmass(GravitationalMeridian))

julia> molarmass(GravitationalMeridian) # ehyl⋅eg-mol⁻¹
$(molarmass(GravitationalMeridian))

julia> luminousefficacy(GravitationalMeridian) # lm⋅s⋅em⁻¹⋅kegf⁻¹
$(luminousefficacy(GravitationalMeridian))
```
""" GravitationalMeridian=#

@doc """
$(systext(:MTS,"EntropySystem(SI2019,𝟏,𝟏,kilo)"))

Metre-tonne-second `UnitSystem` variant of `Metric` system.

```Julia
julia> boltzmann(MTS) # kJ⋅K⁻¹
$(boltzmann(MTS))

julia> planckreduced(MTS) # kJ⋅s⋅rad⁻¹
$(planckreduced(MTS))

julia> lightspeed(MTS) # m⋅s⁻¹
$(lightspeed(MTS))

julia> vacuumpermeability(MTS) # kH⋅m⁻¹
$(vacuumpermeability(MTS))

julia> electronmass(MTS) # t
$(electronmass(MTS))

julia> molarmass(MTS) # t⋅mol⁻¹
$(molarmass(MTS))

julia> luminousefficacy(MTS) # lm⋅kW⁻¹
$(luminousefficacy(MTS))
```
""" MTS

@doc """
$(systext(:KKH,"EntropySystem(Metric,HOUR,kilo,𝟏)"))

Kilometer-kilogram-hour `UnitSystem` variant of `Metric` system.

```Julia
julia> boltzmann(KKH) # kg⋅km²⋅h⁻²⋅K⁻¹
$(boltzmann(KKH))

julia> planckreduced(KKH) # kg⋅km²⋅h⁻¹
$(planckreduced(KKH))

julia> lightspeed(KKH) # km⋅hr⁻¹
$(lightspeed(KKH))

julia> vacuumpermeability(KKH) # kg⋅km⋅C⁻²
$(vacuumpermeability(KKH))

julia> electronmass(KKH) # kg
$(electronmass(KKH))

julia> molarmass(KKH) # kg⋅mol⁻¹
$(molarmass(KKH))

julia> luminousefficacy(KKH) # lm⋅h³⋅kg⁻¹⋅km⁻²
$(luminousefficacy(KKH))
```
""" KKH

@doc """
$(systext(:IAU☉,"EntropySystem(Metric,DAY,au,GM☉/G)"))

Solar `UnitSystem` defined by International Astronomical Union and `solarmass`.

```Julia
julia> boltzmann(IAU) # M⊙⋅au²⋅D⁻²⋅K⁻¹
$(boltzmann(IAU))

julia> planckreduced(IAU) # M⊙⋅au²⋅D⁻¹⋅rad⁻¹
$(planckreduced(IAU))

julia> lightspeed(IAU) # au⋅D⁻¹
$(lightspeed(IAU))

julia> vacuumpermeability(IAU) # M⊙⋅au²⋅C⁻²
$(vacuumpermeability(IAU))

julia> electronmass(IAU) # M⊙
$(electronmass(IAU))

julia> molarmass(IAU) # M☉⋅mol⁻¹
$(molarmass(IAU))

julia> luminousefficacy(IAU) # lm⋅D³⋅M☉⁻¹⋅au⁻²
$(luminousefficacy(IAU))

julia> gaussgravitation(IAU) # D⁻¹
$(gaussgravitation(IAU))
```
""" IAU☉, IAU

@doc """
$(systext(:IAUE,"EntropySystem(Metric,DAY,LD,GME/G)"))

Astronomical (Earth) `UnitSystem` defined by `lunardistance` around the `earthmass`.

```Julia
julia> boltzmann(IAUE) # ME⋅LD²⋅D⁻²⋅K⁻¹
$(boltzmann(IAUE))

julia> planckreduced(IAUE) # ME⋅LD²⋅D⁻¹⋅rad⁻¹
$(planckreduced(IAUE))

julia> lightspeed(IAUE) # LD⋅D⁻¹
$(lightspeed(IAUE))

julia> vacuumpermeability(IAUE) # ME⋅LD²⋅C⁻²
$(vacuumpermeability(IAUE))

julia> electronmass(IAUE) # ME
$(electronmass(IAUE))

julia> molarmass(IAUE) # ME⋅mol⁻¹
$(molarmass(IAUE))

julia> luminousefficacy(IAUE) # lm⋅D³⋅ME⁻¹⋅LD⁻²
$(luminousefficacy(IAUE))

julia> turn(IAU)/gaussianmonth(IAU) # D⁻¹
$(turn(IAU)/gaussianmonth(IAU))
```
""" IAUE

@doc """
$(systext(:IAUJ,"EntropySystem(Metric,DAY,JD,GMJ/G)"))

Astronomical (Jupiter) `UnitSystem` defined by `jupiterdistance` around the `solarmass`.

```Julia
julia> boltzmann(IAUJ) # MJ⋅JD²⋅D⁻²⋅K⁻¹
$(boltzmann(IAUJ))

julia> planckreduced(IAUJ) # MJ⋅JD²⋅D⁻¹⋅rad⁻¹
$(planckreduced(IAUJ))

julia> lightspeed(IAUJ) # JD⋅D⁻¹
$(lightspeed(IAUJ))

julia> vacuumpermeability(IAUJ) # MJ⋅JD²⋅C⁻²
$(vacuumpermeability(IAUJ))

julia> electronmass(IAUJ) # MJ
$(electronmass(IAUJ))

julia> molarmass(IAUJ) # MJ⋅mol⁻¹
$(molarmass(IAUJ))

julia> luminousefficacy(IAUJ) # lm⋅D³⋅MJ⁻¹⋅JD⁻²
$(luminousefficacy(IAUJ))

julia> sqrt(gravitation(IAUJ)*solarmass(IAUJ)/jupiterdistance(IAUJ)^3) # D⁻¹
$(sqrt(gravitation(IAUJ)*solarmass(IAUJ)/jupiterdistance(IAUJ)^3))
```
""" IAUJ

#=@doc """
$(systext(:Astronomical,"EntropySystem(Metric,𝟏,𝟏,𝟏/gravitation(Metric))"))

Astronomical `UnitSystem` defined by making the `newton` gravitational constant 1.

```Julia
julia> boltzmann(Astronomical)
$(boltzmann(Astronomical))

julia> planckreduced(Astronomical)
$(planckreduced(Astronomical))

julia> lightspeed(Astronomical)
$(lightspeed(Astronomical))

julia> vacuumpermeability(Astronomical)
$(vacuumpermeability(Astronomical))

julia> electronmass(Astronomical)
$(electronmass(Astronomical))

julia> molarmass(Astronomical)
$(molarmass(Astronomical))

julia> luminousefficacy(Astronomical)
$(luminousefficacy(Astronomical))

julia> gravitation(Astronomical)
$(gravitation(Astronomical))
```
""" Astronomical=#

@doc """
$(systext(:Hubble,"AstronomicalSystem(Metric,th,𝘤*th,mₑ)"))

Hubble `UnitSystem` defined by `hubble` parameter.

```Julia
julia> boltzmann(Hubble)
$(boltzmann(Hubble))

julia> planckreduced(Hubble)
$(planckreduced(Hubble))

julia> lightspeed(Hubble)
$(lightspeed(Hubble))

julia> vacuumpermeability(Hubble)
$(vacuumpermeability(Hubble))

julia> electronmass(Hubble)
$(electronmass(Hubble))

julia> molarmass(Hubble)
$(molarmass(Hubble))

julia> luminousefficacy(Hubble)
$(luminousefficacy(Hubble))

julia> hubble(Hubble)
$(hubble(Hubble))

julia> cosmological(Hubble)
$(cosmological(Hubble))
```
""" Hubble

@doc """
$(systext(:Cosmological,"AstronomicalSystem(Metric,lc/𝘤,lc,mc)"))

Cosmological scale `UnitSystem` defined by `darkenergydensity`.

```Julia
julia> boltzmann(Cosmological)
$(boltzmann(Cosmological))

julia> planckreduced(Cosmological)
$(planckreduced(Cosmological))

julia> lightspeed(Cosmological)
$(lightspeed(Cosmological))

julia> vacuumpermeability(Cosmological)
$(vacuumpermeability(Cosmological))

julia> electronmass(Cosmological)
$(electronmass(Cosmological))

julia> molarmass(Cosmological)
$(molarmass(Cosmological))

julia> luminousefficacy(Cosmological)
$(luminousefficacy(Cosmological))

julia> hubble(Cosmological)
$(hubble(Cosmological))

julia> cosmological(Cosmological)
$(cosmological(Cosmological))
```
""" Cosmological

@doc """
$(systext(:CosmologicalQuantum,"AstronomicalSystem(Metric,tcq,lcq,mcq)"))

Cosmological quantum scale `UnitSystem` defined by `darkenergydensity`.

```Julia
julia> boltzmann(CosmologicalQuantum)
$(boltzmann(CosmologicalQuantum))

julia> planckreduced(CosmologicalQuantum)
$(planckreduced(CosmologicalQuantum))

julia> lightspeed(CosmologicalQuantum)
$(lightspeed(CosmologicalQuantum))

julia> vacuumpermeability(CosmologicalQuantum)
$(vacuumpermeability(CosmologicalQuantum))

julia> electronmass(CosmologicalQuantum)
$(electronmass(CosmologicalQuantum))

julia> molarmass(CosmologicalQuantum)
$(molarmass(CosmologicalQuantum))

julia> luminousefficacy(CosmologicalQuantum)
$(luminousefficacy(CosmologicalQuantum))
```
""" CosmologicalQuantum

@doc """
$(systext(:British,"RankineSystem(Metric,ft,slug)"))

British Gravitational `UnitSystem` historically used by Britain and United States.

```Julia
julia> boltzmann(British) # ft⋅lb⋅°R⁻¹
$(boltzmann(British))

julia> planckreduced(British) # ft⋅lb⋅s⋅rad⁻¹
$(planckreduced(British))

julia> lightspeed(British) # ft⋅s⁻¹
$(lightspeed(British))

julia> vacuumpermeability(British) # slug⋅ft⋅C⁻²
$(vacuumpermeability(British))

julia> electronmass(British) # slugs
$(electronmass(British))

julia> molarmass(British) # slug⋅slug-mol⁻¹
$(molarmass(British))

julia> luminousefficacy(British) # lm⋅s⋅ft⁻¹⋅lb⁻¹
$(luminousefficacy(British))
```
""" British, BritishGravitational, BG

@doc """
$(systext(:English,"RankineSystem(Metric,ft,lb,g₀/ft)"))

English Engineering `UnitSystem` historically used in the United States of America.

```Julia
julia> boltzmann(English) # ft⋅lbf⋅°R⁻¹
$(boltzmann(English))

julia> planckreduced(English) # ft⋅lbf⋅s⋅rad⁻¹
$(planckreduced(English))

julia> lightspeed(English) # ft⋅s⁻¹
$(lightspeed(English))

julia> vacuumpermeability(English) # lbm⋅ft⋅C⁻²
$(vacuumpermeability(English))

julia> electronmass(English) # lbm
$(electronmass(English))

julia> molarmass(English) # lbm⋅lb-mol⁻¹
$(molarmass(English))

julia> luminousefficacy(English) # lm⋅s⋅ft⁻¹⋅lbf⁻¹
$(luminousefficacy(English))

julia> gravity(English) # lbm⋅ft⋅lbf⁻¹⋅s⁻²
$(gravity(English))
```
""" English, EnglishEngineering, EE

@doc """
$(systext(:Survey,"RankineSystem(Metric,ftUS,lb,g₀/ftUS)"))

English Engineering `UnitSystem` based on the geophysical US survey foot (1200/3937).

```Julia
julia> boltzmann(Survey) # ftUS⋅lbf⋅°R⁻¹
$(boltzmann(Survey))

julia> planckreduced(Survey) # ftUS⋅lbf⋅s⋅rad⁻¹
$(planckreduced(Survey))

julia> lightspeed(Survey) # ftUS⋅s⁻¹
$(lightspeed(Survey))

julia> vacuumpermeability(Survey) # lbm⋅ftUS⋅C⁻²
$(vacuumpermeability(Survey))

julia> electronmass(Survey) # lbm
$(electronmass(Survey))

julia> molarmass(Survey) # lbm⋅lb-mol⁻¹
$(molarmass(Survey))

julia> luminousefficacy(Survey) # lm⋅s⋅ft⁻¹⋅lbf⁻¹
$(luminousefficacy(Survey))

julia> gravity(Survey) # lbm⋅ftUS⋅lbf⁻¹⋅s⁻²
$(gravity(Survey))
```
""" Survey, EnglishUS

@doc """
$(systext(:FPS,"RankineSystem(Metric,ft,lb)"))

Absolute English `UnitSystem` based on the foot, pound, second, and poundal.

```Julia
julia> boltzmann(FPS) # ft⋅pdl⋅°R⁻¹
$(boltzmann(FPS))

julia> planckreduced(FPS) # ft⋅pdl⋅s⋅rad⁻¹
$(planckreduced(FPS))

julia> lightspeed(FPS) # ft⋅s⁻¹
$(lightspeed(FPS))

julia> vacuumpermeability(FPS) # lb⋅ft⋅C⁻²
$(vacuumpermeability(FPS))

julia> electronmass(FPS) # lb
$(electronmass(FPS))

julia> molarmass(FPS) # lb⋅lb-mol⁻¹
$(molarmass(FPS))

julia> luminousefficacy(FPS) # lm⋅s³⋅lb⁻¹⋅ft⁻²
$(luminousefficacy(FPS))
```
""" FPS, AbsoluteEnglish, AE

@doc """
$(systext(:IPS,"RankineSystem(Metric,ft/𝟐^2/𝟑,lb*g₀*𝟐^2*𝟑/ft)"))

British Gravitational `UnitSystem` historically used in the United States of America.

```Julia
julia> boltzmann(IPS) # in⋅lb⋅°R⁻¹
$(boltzmann(IPS))

julia> planckreduced(IPS) # in⋅lb⋅s⋅rad⁻¹
$(planckreduced(IPS))

julia> lightspeed(IPS) # in⋅s⁻¹
$(lightspeed(IPS))

julia> vacuumpermeability(IPS) # slinch⋅in⋅C⁻²
$(vacuumpermeability(IPS))

julia> electronmass(IPS) # slinch
$(electronmass(IPS))

julia> molarmass(IPS) # slinch⋅slinch-mol⁻¹
$(molarmass(IPS))

julia> luminousefficacy(IPS) # lm⋅s⋅in⁻¹⋅lb⁻¹
$(luminousefficacy(IPS))
```
""" IPS

@doc """
$(systext(:FFF,"EntropySystem(Metric,𝟐*𝟕*DAY,fur,𝟐*𝟑^2*𝟓*lb,°R,0,𝟏)"))

Furlong–firkin–fortnight `FFF` is a humorous `UnitSystem` based on unusal impractical units.

```Julia
julia> boltzmann(FFF) # fir⋅fur²⋅ftn⁻²⋅F⁻¹
$(boltzmann(FFF))

julia> planckreduced(FFF) # fir⋅fur²⋅ftn⁻¹⋅rad⁻¹
$(planckreduced(FFF))

julia> lightspeed(FFF) # fur⋅ftn⁻¹
$(lightspeed(FFF))

julia> vacuumpermeability(FFF) # fir⋅fur⋅Inf⁻²
$(vacuumpermeability(FFF))

julia> electronmass(FFF) # fir
$(electronmass(FFF))

julia> molarmass(FFF) # fir⋅fir-mol⁻¹
$(molarmass(FFF))

julia> luminousefficacy(FFF) # lm⋅ftn³⋅fir⁻¹⋅fur⁻²
$(luminousefficacy(FFF))
```
""" FFF

@doc """
$(systext(:MPH,"EntropySystem(FPS,HOUR,mi,𝟏)"))

Mile-pound-hour specification based on `FPS` absolute `UnitSystem`.

```Julia
julia> boltzmann(MPH) # lbf⋅mi²⋅hr⁻²⋅F⁻¹
$(boltzmann(MPH))

julia> planckreduced(MPH) # lbf⋅mi²⋅hr⁻¹⋅rad⁻¹
$(planckreduced(MPH))

julia> lightspeed(MPH) # mi⋅hr⁻¹
$(lightspeed(MPH))

julia> vacuumpermeability(MPH) # lbm⋅mi⋅C⁻²
$(vacuumpermeability(MPH))

julia> electronmass(MPH) # lbm
$(electronmass(MPH))

julia> molarmass(MPH) # lbm⋅lb-mol⁻¹
$(molarmass(MPH))

julia> luminousefficacy(MPH) # lm⋅h³⋅lb⁻¹⋅mi⁻²
$(luminousefficacy(MPH))
```
""" MPH

@doc """
$(systext(:Nautical,"EntropySystem(Metric,HOUR,nm,em^3,𝟏,τ*𝟑^3/𝟐^10/𝟓^12,milli)"))

Nautical miles, kilo-earthgram, hour specification based on `Meridian` definition.

```Julia
julia> greatcircle(Nautical) # nm
$(greatcircle(Nautical))

julia> boltzmann(Nautical) # keg⋅nm²⋅hr⁻²⋅K⁻¹
$(boltzmann(Nautical))

julia> planckreduced(Nautical) # keg⋅nm²⋅hr⁻¹⋅rad⁻¹
$(planckreduced(Nautical))

julia> lightspeed(Nautical) # nm⋅hr⁻¹
$(lightspeed(Nautical))

julia> vacuumpermeability(Nautical) # keg⋅nm⋅eC⁻²
$(vacuumpermeability(Nautical))

julia> electronmass(Nautical) # keg
$(electronmass(Nautical))

julia> molarmass(Nautical) # keg⋅eg-mol⁻¹
$(molarmass(Nautical))

julia> luminousefficacy(Nautical) # lm⋅h³⋅keg⁻¹⋅nm⁻²
$(luminousefficacy(Nautical))
```
""" Nautical

# natural unit system docs

textunits(U,S) = """
```Julia
julia> boltzmann($S)
$(boltzmann(U))

julia> planckreduced($S)
$(planckreduced(U))

julia> lightspeed($S)
$(lightspeed(U))

julia> vacuumpermeability($S)
$(vacuumpermeability(U))

julia> electronmass($S)
$(electronmass(U))
```
"""

@doc """
$(systext(:Planck,"UnitSystem(𝟏,𝟏,𝟏,𝟏,√(𝟐*τ*αG))"))

Planck `UnitSystem` with the `electronmass` value `√(4π*αG)` using gravitational coupling.

$(textunits(Planck,:Planck))
""" Planck

@doc """
$(systext(:PlanckGauss,"UnitSystem(𝟏,𝟏,𝟏,𝟐*τ,√αG)"))

Planck (Gauss) `UnitSystem` with `permeability` of `4π` and `electronmass` coupling `√αG`.

$(textunits(PlanckGauss,:PlanckGauss))

The well known `PlanckGauss` values for `length`, `time`, `mass`, and `temperature` are:
```Julia
julia> length(PlanckGauss,SI2019) # ℓP
$(length(PlanckGauss,SI2019))

julia> time(PlanckGauss,SI2019) # tP
$(time(PlanckGauss,SI2019))

julia> mass(PlanckGauss,SI2019) # mP
$(mass(PlanckGauss,SI2019))

julia> temperature(PlanckGauss,SI2019) # TP
$(temperature(PlanckGauss,SI2019))
```
""" PlanckGauss

@doc """
$(systext(:Stoney,"UnitSystem(𝟏,𝟏/α,𝟏,𝟐*τ,√(αG/α))"))

Stoney `UnitSystem` with `permeability` of `4π` and `electronmass` coupling `√(αG/α)`.

$(textunits(Stoney,:Stoney))

The well known `Stoney` values for `length`, `time`, `mass`, and `charge` are:
```Julia
julia> length(Stoney,SI2019) # lS
$(length(Stoney,SI2019))

julia> time(Stoney,SI2019) # tS
$(time(Stoney,SI2019))

julia> mass(Stoney,SI2019) # mS
$(mass(Stoney,SI2019))

julia> charge(Stoney,SI2019) # qS
$(charge(Stoney,SI2019))
```
""" Stoney

@doc """
$(systext(:Hartree,"UnitSystem(𝟏,𝟏,𝟏/α,𝟐*τ*α^2,𝟏)"))

Hartree atomic `UnitSystem` based on `bohr` radius and `elementarycharge` scale.

$(textunits(Hartree,:Hartree))

The well known `Hartree` atomic unit values for `length`, `time`, `mass`, and `charge` are:
```Julia
julia> length(Hartree,SI2019) # lA
$(length(Hartree,SI2019))

julia> time(Hartree,SI2019) # tA
$(time(Hartree,SI2019))

julia> mass(Hartree,SI2019) # mA
$(mass(Hartree,SI2019))

julia> charge(Hartree,SI2019) # qA
$(charge(Hartree,SI2019))
```
""" Hartree

@doc """
$(systext(:Rydberg,"UnitSystem(𝟏,𝟏,𝟐/α,τ/𝟐*α^2,𝟏/𝟐)"))

Rydberg `UnitSystem` with `lightspeed` of `𝟐/α` and `permeability` of `π*α^2`.

$(textunits(Rydberg,:Rydberg))

The well known `Rydberg` atomic unit values for `length`, `time`, `mass`, and `charge` are:
```Julia
julia> length(Rydberg,SI2019) # lR
$(length(Rydberg,SI2019))

julia> time(Rydberg,SI2019) # tR
$(time(Rydberg,SI2019))

julia> mass(Rydberg,SI2019) # mR
$(mass(Rydberg,SI2019))

julia> charge(Rydberg,SI2019) # qR
$(charge(Rydberg,SI2019))
```
""" Rydberg

@doc """
$(systext(:Schrodinger,"UnitSystem(𝟏,𝟏,𝟏/α,𝟐*τ*α^2,√(αG/α))"))

Schrodinger `UnitSystem` with `permeability` of `4π/αinv^2` and `electronmass` of `√(αG*αinv)`.

$(textunits(Schrodinger,:Schrodinger))
""" Schrodinger

@doc """
$(systext(:Electronic,"UnitSystem(𝟏,𝟏/α,𝟏,𝟐*τ,𝟏)"))

Electronic `UnitSystem` with `planckreduced` of `1/α` and `permeability` of `4π`.

$(textunits(Electronic,:Electronic))
""" Electronic

@doc """
$(systext(:Natural,"UnitSystem(𝟏,𝟏,𝟏,𝟏,𝟏)"))

Natural `UnitSystem` with all primary constants having unit value.

$(textunits(Natural,:Natural))

The well known `Natural` values for `length`, `time`, `mass`, and `charge` are:
```Julia
julia> length(Natural,SI2019)
$(length(Natural,SI2019))

julia> time(Natural,SI2019)
$(time(Natural,SI2019))

julia> mass(Natural,SI2019)
$(mass(Natural,SI2019))

julia> charge(Natural,SI2019)
$(charge(Natural,SI2019))
```
""" Natural

@doc """
$(systext(:NaturalGauss,"UnitSystem(𝟏,𝟏,𝟏,𝟐*τ,𝟏)"))

Natural (Gauss) `UnitSystem` with the Gaussian `permeability` value of `4π`.

$(textunits(NaturalGauss,:NaturalGauss))
""" NaturalGauss

@doc """
$(systext(:QCD,"UnitSystem(𝟏,𝟏,𝟏,𝟏,𝟏/μₚₑ)"))

Qunatum chromodynamics `UnitSystem` based on the `protonmass` scale.

$(textunits(QCD,:QCD))

The well known `QCD` values for `length`, `time`, `mass`, and `charge` are:
```Julia
julia> length(QCD,SI2019) # lQCD
$(length(QCD,SI2019))

julia> time(QCD,SI2019) # tQCD
$(time(QCD,SI2019))

julia> mass(QCD,SI2019) # mQCD
$(mass(QCD,SI2019))

julia> charge(QCD,SI2019) # qQCD
$(charge(QCD,SI2019))
```
""" QCD

@doc """
$(systext(:QCDGauss,"UnitSystem(𝟏,𝟏,𝟏,𝟐*τ,𝟏/μₚₑ)"))

Qunatum chromodynamics (Gauss) `UnitSystem` based on the `protonmass` scale.

$(textunits(QCDGauss,:QCDGauss))

The well known `QCDGauss` values for `length`, `time`, `mass`, and `charge` are:
```Julia
julia> length(QCDGauss,SI2019) # lQCD
$(length(QCDGauss,SI2019))

julia> time(QCDGauss,SI2019) # tQCD
$(time(QCDGauss,SI2019))

julia> mass(QCDGauss,SI2019) # mQCD
$(mass(QCDGauss,SI2019))

julia> charge(QCDGauss,SI2019) # qQCD
$(charge(QCDGauss,SI2019))
```
""" QCDGauss

@doc """
$(systext(:QCDoriginal,"UnitSystem(𝟏,𝟏,𝟏,𝟐*τ*α,𝟏/μₚₑ)"))

Qunatum chromodynamics (original) `UnitSystem` scaled by `protonmass` and `elementarycharge`.

$(textunits(QCDoriginal,:QCDoriginal))

The well known `QCDoriginal` values for `length`, `time`, `mass`, and `charge` are:
```Julia
julia> length(QCDoriginal,SI2019) # lQCD
$(length(QCDoriginal,SI2019))

julia> time(QCDoriginal,SI2019) # tQCD
$(time(QCDoriginal,SI2019))

julia> mass(QCDoriginal,SI2019) # mQCD
$(mass(QCDoriginal,SI2019))

julia> charge(QCDoriginal,SI2019) # qQCD
$(charge(QCDoriginal,SI2019))
```
""" QCDoriginal
