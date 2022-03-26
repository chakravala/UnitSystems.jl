
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

export slug, ft, KJ1990, KJ2014, RK1990, RK2014, mₑ1990, mₑ2014, temp, units, °R
export slugs, kilograms, lbm, meters, feet, rankine, kelvin, moles, molecules, universal
export Universe, UnitSystem, US, universe, HOUR, DAY, th, lc, mc, tcq, lcq, mcq
export similitude, 𝟙, F, M, L, T, Q, Θ, N, J, A, Λ, C

# == Metric is different
const κ = einstein(SI2019)
const σ = stefan(SI2019) #
const μB = magneton(SI2019) #
const ε₀ = vacuumpermittivity(SI2019) #
const kₑ = coulomb(SI2019) #
const mₚ = protonmass(SI2019)
const mᵤ = atomicmass(SI2019)
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
const universal = universalgas

# constant aliases

const US,mpe, meu, mpu, ainv, aG = UnitSystem,μₚₑ, μₑᵤ, μₚᵤ, αinv, αG
const Mu,Ru,SB,hh,cc,m0,e0,ke,me,mp,mu,ee,FF,Z0,G0,Eh,a0,re,g0,lP,aL,ϵ₀ = Mᵤ,Rᵤ,σ,𝘩,𝘤,μ₀,ε₀,kₑ,mₑ,mₚ,mᵤ,𝘦,𝔉,Z₀,G₀,Eₕ,a₀,rₑ,g₀,ℓP,αL,ε₀
export κ, G, GG, NA, kB, Rᵤ, σ, 𝘩, ħ, 𝘤, μ₀, ε₀, kₑ, mₑ, mₚ, mᵤ, 𝘦, 𝔉, Φ₀, Z₀, G₀, Eₕ, R∞, a₀, rₑ, KJ, RK, Ru, SB, hh, cc, m0, e0, ke, me, mp, mu, ee, FF, Z0, G0, Eh, a0, re, μB
export αG, αinv, μₚₑ, μₑᵤ, μₚᵤ, mpe, meu, mpu, mP, δμ₀, Mᵤ, Mu, RH, Ry, ΔνCs, Kcd, ainv
export cal, kcal, calₜₕ, kcalₜₕ, calᵢₜ, kcalᵢₜ, ℓP, g₀, g0, atm, lbm, BTUJ, BTUftlb, aG
export lP, tP, TP, lS, tS, mS, qS, lA, tA, mA, qA, lQCD, tQCD, mQCD, ϵ₀, αL, aL

# engineering unit systems docs

@doc """
    Metric = MetricSystem(𝟏𝟎^-3,𝟐*τ/𝟏𝟎^7)

Systeme International d'Unites (the SI units) adopted as the preferred `UnitSystem`.

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
```
""" Metric, MKS

@doc """
    SI2019 = MetricSystem(Mᵤ,μ₀)

Systeme International d'Unites (the SI units) with `μ₀` for a tuned `charge` exactly.

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
```
""" SI2019, SI

@doc """
    MetricEngineering = MetricSystem(𝟏𝟎^-3,𝟐*τ/𝟏𝟎^7,Rᵤ,g₀)

Systeme International d'Unites (the SI units) based on kilogram and kilopond units.

```Julia
julia> boltzmann(MetricEngineering) # J⋅K⁻¹
$(boltzmann(MetricEngineering))

julia> planckreduced(MetricEngineering) # J⋅s⋅rad⁻¹
$(planckreduced(MetricEngineering))

julia> lightspeed(MetricEngineering) # m⋅s⁻¹
$(lightspeed(MetricEngineering))

julia> vacuumpermeability(MetricEngineering) # H⋅m⁻¹
$(vacuumpermeability(MetricEngineering))

julia> electronmass(MetricEngineering) # kg
$(electronmass(MetricEngineering))

julia> molarmass(MetricEngineering) # kg⋅mol⁻¹
$(molarmass(MetricEngineering))

julia> gravity(MetricEngineering) # kg⋅m⋅N⁻¹⋅s⁻²
$(gravity(MetricEngineering))
```
""" MetricEngineering, ME

@doc """
    SI2019Engineering = MetricSystem(Mᵤ,μ₀,Rᵤ,g₀)

Systeme International d'Unites (the SI units) based on kilogram and kilopond units.

```Julia
julia> boltzmann(SI2019Engineering) # J⋅K⁻¹
$(boltzmann(SI2019Engineering))

julia> planckreduced(SI2019Engineering) # J⋅s⋅rad⁻¹
$(planckreduced(SI2019Engineering))

julia> lightspeed(SI2019Engineering) # m⋅s⁻¹
$(lightspeed(SI2019Engineering))

julia> vacuumpermeability(SI2019Engineering) # H⋅m⁻¹
$(vacuumpermeability(SI2019Engineering))

julia> electronmass(SI2019Engineering) # kg
$(electronmass(SI2019Engineering))

julia> molarmass(SI2019Engineering) # kg⋅mol⁻¹
$(molarmass(SI2019Engineering))

julia> gravity(SI2019Engineering) # kg⋅m⋅N⁻¹⋅s⁻²
$(gravity(SI2019Engineering))
```
""" SI2019Engineering, SIE

@doc """
    SI1976 = MetricSystem(𝟏𝟎^-3,𝟐*τ/𝟏𝟎^7,8.31432)

Systeme International d'Unites (the SI units) with universal gas constant of `8.31432`.

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
```
""" SI1976

@doc """
    CODATA = ConventionalSystem(RK2014,KJ2014,Rᵤ2014)

Metric `UnitSystem` based on Committee on Data of the International Science Council.

```Julia
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
```
""" CODATA

@doc """
    Conventional = ConventionalSystem(RK1990,KJ2014)

Conventional electronic `UnitSystem` with 1990 tuned `josephson` and `klitzing` constants.

```Julia
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
```
""" Conventional

@doc """
    International = ElectricSystem(Metric,Ωᵢₜ,Vᵢₜ)

International `UnitSystem` with United States measurements of `Ωᵢₜ` and `Vᵢₜ`.

```Julia
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
```
""" International

@doc """
    InternationalMean = ElectricSystem(Metric,1.00049,1.00034)

International `UnitSystem` with mean measurements of `Ωᵢₜ` and `Vᵢₜ`.

```Julia
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
```
""" InternationalMean

cgstext(US,AMP,cgs=eval(US)) = """
```Julia
julia> boltzmann($US) # erg⋅K⁻¹
$(boltzmann(cgs))

julia> planckreduced($US) # erg⋅s⋅rad⁻¹
$(planckreduced(cgs))

julia> lightspeed($US) # cm⋅s⁻¹
$(lightspeed(cgs))

julia> vacuumpermeability($US) # statH⋅cm⁻¹
$(vacuumpermeability(cgs))

julia> electronmass($US) # g
$(electronmass(cgs))

julia> molarmass($US) # g⋅mol⁻¹
$(molarmass(cgs))

julia> rationalization($US)
$(rationalization(cgs))
```
"""

for U ∈ (:CGSm,:CGSe,:EMU,:ESU)
    (EU,AMP) = QuoteNode.(U ∉ (:CGSe,:ESU) ? (:EMU,:Bi) : (:ESU,:statA))
@eval @doc """
    $($(QuoteNode(U))) = GaussSystem(Metric,$($EU≠:EMU ? "(𝟏𝟎*𝘤)^-2" : 𝟏),𝟐*τ)

Centimetre-gram-second `UnitSystem` variant based on `$($EU)` (non-rationalized).

$(cgstext($(QuoteNode(U)),$AMP))
""" $U

U ∉ (:CGSm,:CGSe) && @eval @doc """
    $(Symbol($(QuoteNode(U)),:2019)) = EntropySystem(SI2019,𝟏,0.01,0.001,𝟏,$($EU≠:EMU ? "1e3*μ₀/𝘤^2" : "1e7*μ₀"))

Centimetre-gram-second `UnitSystem` variant of tuned `SI2019` based on `$($EU)` (rationalized).

$(cgstext(Symbol($(QuoteNode(U)),:2019),$AMP))
""" $(Symbol(U,:2019))
end

@doc """
    Thomson = GaussSystem(Metric,𝟏,𝟐*τ,𝟏/𝟐)

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

julia> rationalization(Thomson)
$(rationalization(Thomson))

julia> lorentz(Thomson)
$(lorentz(Thomson))
```
""" Thomson

@doc """
    Gauss = GaussSystem(Metric,𝟏,𝟐*τ,𝟏𝟎^-2/𝘤)

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

julia> rationalization(Gauss)
$(rationalization(Gauss))

julia> lorentz(Gauss)
$(lorentz(Gauss))
```
""" Gauss, CGS

@doc """
    LorentzHeaviside = GaussSystem(Metric,𝟏,𝟏,𝟏𝟎^-2/𝘤)

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

julia> rationalization(LorentzHeaviside)
$(rationalization(LorentzHeaviside))

julia> lorentz(LorentzHeaviside)
$(lorentz(LorentzHeaviside))
```
""" LorentzHeaviside, HLU

@doc """
    Kennelly = GaussSystem(Metric,𝟏𝟎^-7,𝟐*τ,𝟏,𝟏,𝟏)

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

julia> rationalization(Kennelly)
$(rationalization(Kennelly))
```
""" Kennelly

@doc """
    GravitationalMetric = EntropySystem(Metric,𝟏,𝟏,g₀)

Systeme International d'Unites (the SI units) based on hyl and kilopond units.

```Julia
julia> boltzmann(GravitationalMetric) # J⋅K⁻¹
$(boltzmann(GravitationalMetric))

julia> planckreduced(GravitationalMetric) # J⋅s⋅rad⁻¹
$(planckreduced(GravitationalMetric))

julia> lightspeed(GravitationalMetric) # m⋅s⁻¹
$(lightspeed(GravitationalMetric))

julia> vacuumpermeability(GravitationalMetric) # H⋅m⁻¹
$(vacuumpermeability(GravitationalMetric))

julia> electronmass(GravitationalMetric) # hyl
$(electronmass(GravitationalMetric))

julia> molarmass(GravitationalMetric) # hyl⋅mol⁻¹
$(molarmass(GravitationalMetric))
```
""" GravitationalMetric, GM

@doc """
    GraviationalSI2019 = EntropySystem(SI2019,𝟏,𝟏,g₀)

Systeme International d'Unites (the SI units) based on hyl and kilopond units.

```Julia
julia> boltzmann(GravitationalSI2019) # J⋅K⁻¹
$(boltzmann(GravitationalSI2019))

julia> planckreduced(GravitationalSI2019) # J⋅s⋅rad⁻¹
$(planckreduced(GravitationalSI2019))

julia> lightspeed(GravitationalSI2019) # m⋅s⁻¹
$(lightspeed(GravitationalSI2019))

julia> vacuumpermeability(SI2019Engineering) # H⋅m⁻¹
$(vacuumpermeability(GravitationalSI2019))

julia> electronmass(SI2019Engineering) # hyl
$(electronmass(GravitationalSI2019))

julia> molarmass(SI2019) # hyl⋅mol⁻¹
$(molarmass(SI2019Engineering))
```
""" GravitationalSI2019, GSI, GSI2019

@doc """
    MTS = EntropySystem(SI2019,𝟏,𝟏,𝟏𝟎^3)

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
```
""" MTS

@doc """
    KKH = EntropySystem(Metric,HOUR,𝟏𝟎^3,𝟏)

Kilometer-kilogram-hour `UnitSystem` variant of `Metric` system.

```Julia
julia> boltzmann(KKH)
$(boltzmann(KKH))

julia> planckreduced(KKH)
$(planckreduced(KKH))

julia> lightspeed(KKH) # km⋅hr⁻¹
$(lightspeed(KKH))

julia> vacuumpermeability(KKH)
$(vacuumpermeability(KKH))

julia> electronmass(KKH) # kg
$(electronmass(KKH))

julia> molarmass(KKH) # kg⋅mol⁻¹
$(molarmass(KKH))
```
""" KKH

@doc """
    IAU☉ = EntropySystem(Metric,DAY,au,GM☉/G)

Astronomical (solar) `UnitSystem` defined by International Astronomical Union.

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

julia> newton(IAU)
$(newton(IAU))
```
""" IAU☉, IAU

@doc """
    IAUE = EntropySystem(Metric,DAY,au,GME/G)

Astronomical (Earth) `UnitSystem` defined by International Astronomical Union.

```Julia
julia> boltzmann(IAUE) # ME⋅au²⋅D⁻²⋅K⁻¹
$(boltzmann(IAUE))

julia> planckreduced(IAUE) # ME⋅au²⋅D⁻¹⋅rad⁻¹
$(planckreduced(IAUE))

julia> lightspeed(IAUE) # au⋅D⁻¹
$(lightspeed(IAUE))

julia> vacuumpermeability(IAUE) # ME⋅au²⋅C⁻²
$(vacuumpermeability(IAUE))

julia> electronmass(IAUE) # ME
$(electronmass(IAUE))

julia> molarmass(IAUE) # ME⋅mol⁻¹
$(molarmass(IAUE))
```
""" IAUE

@doc """
    IAUJ = EntropySystem(Metric,DAY,au,GMJ/G)

Astronomical (Jupiter) `UnitSystem` defined by International Astronomical Union.

```Julia
julia> boltzmann(IAUJ) # MJ⋅au²⋅D⁻²⋅K⁻¹
$(boltzmann(IAUJ))

julia> planckreduced(IAUJ) # MJ⋅au²⋅D⁻¹⋅rad⁻¹
$(planckreduced(IAUJ))

julia> lightspeed(IAUJ) # au⋅D⁻¹
$(lightspeed(IAUJ))

julia> vacuumpermeability(IAUJ) # MJ⋅au²⋅C⁻²
$(vacuumpermeability(IAUJ))

julia> electronmass(IAU) # MJ
$(electronmass(IAUJ))

julia> molarmass(IAUJ) # MJ⋅mol⁻¹
$(molarmass(IAUJ))
```
""" IAUJ

@doc """
    Astronomical = AstronomicalSystem(Metric)

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

julia> newton(Astronomical)
$(newton(Astronomical))
```
""" Astronomical

@doc """
    Hubble = EntropySystem(Metric,th,𝘤*th,𝟏)

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

julia> hubble(Hubble)
$(hubble(Hubble))
```
""" Hubble

@doc """
    Cosmological = EntropySystem(Metric,lc/𝘤,lc,mc)

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
```
""" Cosmological

@doc """
    CosmologicalQuantum = EntropySystem(Metric,tcq,lcq,mcq)

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

julia> molarmass(Cosmological)
$(molarmass(Cosmological))
```
""" CosmologicalQuantum

@doc """
    British = RankineSystem(Metric,ft,slug)

British Gravitational `UnitSystem` historically used by Britain and United States.

```Julia
julia> boltzmann(British) # ft⋅lb⋅°R⁻¹
$(boltzmann(British))

julia> planckreduced(British) # ft⋅lb⋅s⋅rad⁻¹
$(planckreduced(British))

julia> lightspeed(British) # ft⋅s⁻¹
$(lightspeed(British))

julia> vacuumpermeability(British) # slug⋅ft²⋅?⁻²
$(vacuumpermeability(British))

julia> electronmass(British) # slugs
$(electronmass(British))

julia> molarmass(British) # slug⋅slug-mol⁻¹
$(molarmass(British))
```
""" British, BritishGravitational, BG

@doc """
    British2019 = RankineSystem(SI2019,ft,slug)

British Gravitational `UnitSystem` historically used by Britain and United States.

```Julia
julia> boltzmann(British2019) # ft⋅lb⋅°R⁻¹
$(boltzmann(British2019))

julia> planckreduced(British2019) # ft⋅lb⋅s⋅rad⁻¹
$(planckreduced(British2019))

julia> lightspeed(British2019) # ft⋅s⁻¹
$(lightspeed(British2019))

julia> vacuumpermeability(British2019) # slug⋅ft²⋅?⁻²
$(vacuumpermeability(British2019))

julia> electronmass(British2019) # slugs
$(electronmass(British2019))

julia> molarmass(British2019) # slug⋅slug-mol⁻¹
$(molarmass(British2019))
```
""" British2019, BritishGravitational2019, BG2019

@doc """
    English = RankineSystem(Metric,ft,lb,g₀/ft)

English Engineering `UnitSystem` historically used in the United States of America.

```Julia
julia> boltzmann(English) # ft⋅lbf⋅°R⁻¹
$(boltzmann(English))

julia> planckreduced(English) # ft⋅lbf⋅s⋅rad⁻¹
$(planckreduced(English))

julia> lightspeed(English) # ft⋅s⁻¹
$(lightspeed(English))

julia> vacuumpermeability(English) # lbm⋅ft²⋅?⁻²
$(vacuumpermeability(English))

julia> electronmass(English) # lbm
$(electronmass(English))

julia> molarmass(English) # lbm⋅lb-mol⁻¹
$(molarmass(English))

julia> gravity(English) # lbm⋅ft⋅lbf⁻¹⋅s⁻²
$(gravity(English))
```
""" English, EnglishEngineering, EE

@doc """
    English2019 = RankineSystem(SI2019,ft,lb,g₀/ft)

English Engineering `UnitSystem` historically used in the United States of America.

```Julia
julia> boltzmann(English2019) # ftUS⋅lbf⋅°R⁻¹
$(boltzmann(English2019))

julia> planckreduced(English2019) # ft⋅lbf⋅s⋅rad⁻¹
$(planckreduced(English2019))

julia> lightspeed(English2019) # ft⋅s⁻¹
$(lightspeed(English2019))

julia> vacuumpermeability(English2019) # lbm⋅ft²⋅?⁻²
$(vacuumpermeability(English2019))

julia> electronmass(English2019) # lbm
$(electronmass(English2019))

julia> molarmass(English2019) # lbm⋅lb-mol⁻¹
$(molarmass(English2019))

julia> gravity(English2019) # lbm⋅ft⋅lbf⁻¹⋅s⁻²
$(gravity(English2019))
```
""" English2019, EnglishEngineering2019, EE2019

@doc """
    Survey = RankineSystem(Metric,ftUS,lb,g₀/ftUS)

English Engineering `UnitSystem` based on the geophysical US survey foot (1200/3937).

```Julia
julia> boltzmann(Survey) # ftUS⋅lbf⋅°R⁻¹
$(boltzmann(Survey))

julia> planckreduced(Survey) # ftUS⋅lbf⋅s⋅rad⁻¹
$(planckreduced(Survey))

julia> lightspeed(Survey) # ftUS⋅s⁻¹
$(lightspeed(Survey))

julia> vacuumpermeability(Survey) # lbm⋅ftUS²⋅?⁻²
$(vacuumpermeability(Survey))

julia> electronmass(Survey) # lbm
$(electronmass(Survey))

julia> molarmass(Survey) # lbm⋅lb-mol⁻¹
$(molarmass(Survey))

julia> gravity(Survey) # lbm⋅ftUS⋅lbf⁻¹⋅s⁻²
$(gravity(Survey))
```
""" Survey, EnglishUS

@doc """
    Survey2019 = RankineSystem(SI2019,ftUS,lb,g₀/ftUS)

English Engineering `UnitSystem` based on the geophysical US survey foot (1200/3937).

```Julia
julia> boltzmann(Survey2019) # ftUS⋅lbf⋅°R⁻¹
$(boltzmann(Survey2019))

julia> planckreduced(Survey2019) # ftUS⋅lbf⋅s⋅rad⁻¹
$(planckreduced(Survey2019))

julia> lightspeed(Survey2019) # ftUS⋅s⁻¹
$(lightspeed(Survey2019))

julia> vacuumpermeability(Survey2019) # lbm⋅ftUS²⋅?⁻²
$(vacuumpermeability(Survey2019))

julia> electronmass(Survey2019) # lbm
$(electronmass(Survey2019))

julia> molarmass(Survey2019) # lbm⋅lb-mol⁻¹
$(molarmass(Survey2019))

julia> gravity(Survey2019) # lbm⋅ftUS⋅lbf⁻¹⋅s⁻²
$(gravity(Survey2019))
```
""" Survey2019

@doc """
    FPS = RankineSystem(Metric,ft,lb)

Absolute English `UnitSystem` based on the foot, pound, second, and poundal.

```Julia
julia> boltzmann(FPS) # ft⋅pdl⋅°R⁻¹
$(boltzmann(FPS))

julia> planckreduced(FPS) # ft⋅pdl⋅s⋅rad⁻¹
$(planckreduced(FPS))

julia> lightspeed(FPS) # ft⋅s⁻¹
$(lightspeed(FPS))

julia> vacuumpermeability(FPS) # lb⋅ft²⋅?⁻²
$(vacuumpermeability(FPS))

julia> electronmass(FPS) # lb
$(electronmass(FPS))

julia> molarmass(FPS) # lb⋅lb-mol⁻¹
$(molarmass(FPS))
```
""" FPS, AbsoluteEnglish, AE

@doc """
    FPS2019 = RankineSystem(SI2019,ft,lb)

Absolute English `UnitSystem` based on the foot, pound, second, and poundal.

```Julia
julia> boltzmann(FPS2019) # ft⋅pdl⋅°R⁻¹
$(boltzmann(FPS2019))

julia> planckreduced(FPS2019) # ft⋅pdl⋅s⋅rad⁻¹
$(planckreduced(FPS2019))

julia> lightspeed(FPS2019) # ft⋅s⁻¹
$(lightspeed(FPS2019))

julia> vacuumpermeability(FPS2019) # lb⋅ft²⋅?⁻²
$(vacuumpermeability(FPS2019))

julia> electronmass(FPS2019) # lb
$(electronmass(FPS2019))

julia> molarmass(FPS2019) # lb⋅lb-mol⁻¹
$(molarmass(FPS2019))
```
""" FPS2019, AE2019, AbsoluteEnglish2019

@doc """
    FFF = EntropySystem(Metric,𝟐*𝟕*DAY,fur,𝟐*𝟑^2*𝟓*lb,°R,0,𝟏)

Furlong–firkin–fortnight `FFF` is a humorous `UnitSystem` based on unusal impractical units.

```Julia
julia> boltzmann(FFF) # fir⋅fur²⋅ftn⁻²⋅F⁻¹
$(boltzmann(FFF))

julia> planckreduced(FFF) # fir⋅fur²⋅ftn⁻¹⋅rad⁻¹
$(planckreduced(FFF))

julia> lightspeed(FFF) # fur⋅ftn⁻¹
$(lightspeed(FFF))

julia> vacuumpermeability(FFF) # fir⋅fur²⋅Inf⁻²
$(vacuumpermeability(FFF))

julia> electronmass(FFF) # fir
$(electronmass(FFF))

julia> molarmass(FFF) # fir⋅fir-mol⁻¹
$(molarmass(FFF))
```
""" FFF

@doc """
    MPH = EntropySystem(English,HOUR,𝟐^5*𝟑*𝟓*𝟏𝟏,𝟏)

Miles, pound, hour specification based on `English` Engineering `UnitSystem`.

```Julia
julia> boltzmann(MPH) # lbf⋅mi²⋅hr⁻²⋅F⁻¹
$(boltzmann(MPH))

julia> planckreduced(MPH) # lbf⋅mi²⋅hr⁻¹⋅rad⁻¹
$(planckreduced(MPH))

julia> lightspeed(MPH) # mi⋅hr⁻¹
$(lightspeed(MPH))

julia> vacuumpermeability(MPH) # lbm⋅mi²⋅Inf⁻²
$(vacuumpermeability(MPH))

julia> electronmass(MPH) # lbm
$(electronmass(MPH))

julia> molarmass(MPH) # lbm⋅lb-mol⁻¹
$(molarmass(MPH))
```
""" MPH

@doc """
    Nautical = EntropySystem(English,HOUR,𝟐^6*𝟓*𝟏𝟗,𝟏)

Nautical miles, pound, hour specification based on `English` Engineering `UnitSystem`.

```Julia
julia> boltzmann(Nautical) # lbf⋅nm²⋅hr⁻²⋅F⁻¹
$(boltzmann(Nautical))

julia> planckreduced(Nautical) # lbf⋅nm²⋅hr⁻¹⋅rad⁻¹
$(planckreduced(Nautical))

julia> lightspeed(Nautical) # nm⋅hr⁻¹
$(lightspeed(Nautical))

julia> vacuumpermeability(Nautical) # lbm⋅nm²⋅Inf⁻²
$(vacuumpermeability(Nautical))

julia> electronmass(Nuatical) # lbm
$(electronmass(Nautical))

julia> molarmass(Nautical) # lbm⋅lb-mol⁻¹
$(molarmass(Nautical))
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
    Planck = UnitSystem(𝟏,𝟏,𝟏,𝟏,√(𝟐*τ*αG))

Planck `UnitSystem` with the `electronmass` value `√(4π*αG)` using gravitational coupling.

$(textunits(Planck,:Planck))
""" Planck

@doc """
    PlanckGauss = UnitSystem(𝟏,𝟏,𝟏,𝟐*τ,√αG)

Planck (Gauss) `UnitSystem` with `permeability` of `4π` and `electronmass` coupling `√αG`.

$(textunits(PlanckGauss,:PlanckGauss))

The well known `PlanckGauss` values for `length`, `time`, `mass`, and `temperature` are:
```Julia
julia> length(PlanckGauss,Metric) # ℓP
$(length(PlanckGauss,Metric))

julia> time(PlanckGauss,Metric) # tP
$(time(PlanckGauss,Metric))

julia> mass(PlanckGauss,Metric) # mP
$(mass(PlanckGauss,Metric))

julia> temperature(PlanckGauss,Metric) # TP
$(temperature(PlanckGauss,Metric))
```
""" PlanckGauss

@doc """
    Stoney = UnitSystem(𝟏,𝟏/α,𝟏,𝟐*τ,√(αG/α)}

Stoney `UnitSystem` with `permeability` of `4π` and `electronmass` coupling `√(αG/α)`.

$(textunits(Stoney,:Stoney))

The well known `Stoney` values for `length`, `time`, `mass`, and `charge` are:
```Julia
julia> length(Stoney,Metric) # lS
$(length(Stoney,Metric))

julia> time(Stoney,Metric) # tS
$(time(Stoney,Metric))

julia> mass(Stoney,Metric) # mS
$(mass(Stoney,Metric))

julia> charge(Stoney,Metric) # qS
$(charge(Stoney,Metric))
```
""" Stoney

@doc """
    Hartree = UnitSystem(𝟏,𝟏,𝟏/α,𝟐*τ*α^2,𝟏)

Hartree atomic `UnitSystem` with `lightspeed` of `αinv` and `permeability` of `𝟐*τ*α^2`.

$(textunits(Hartree,:Hartree))

The well known `Hartree` atomic unit values for `length`, `time`, `mass`, and `charge` are:
```Julia
julia> length(Hartree,Metric) # lA
$(length(Hartree,Metric))

julia> time(Hartree,Metric) # tA
$(time(Hartree,Metric))

julia> mass(Hartree,Metric) # mA
$(mass(Hartree,Metric))

julia> charge(Hartree,Metric) # qA
$(charge(Hartree,Metric))
```
""" Hartree

@doc """
    Rydberg = UnitSystem(𝟏,𝟏,𝟐/α,τ/𝟐*α^2,𝟏/𝟐)

Rydberg `UnitSystem` with `lightspeed` of `𝟐/α` and `permeability` of `π*α^2`.

$(textunits(Rydberg,:Rydberg))
""" Rydberg

@doc """
    Schrodinger = UnitSystem(𝟏,𝟏,𝟏/α,𝟐*τ*α^2,√(αG/α))

Schrodinger `UnitSystem` with `permeability` of `4π/αinv^2` and `electronmass` of `√(αG*αinv)`.

$(textunits(Schrodinger,:Schrodinger))
""" Schrodinger

@doc """
    Electronic = UnitSystem(𝟏,𝟏/α,𝟏,𝟐*τ,𝟏}

Electronic `UnitSystem` with `planckreduced` of `1/α` and `permeability` of `4π`.

$(textunits(Electronic,:Electronic))
""" Electronic

@doc """
    Natural = UnitSystem(𝟏,𝟏,𝟏,𝟏,𝟏)

Natural `UnitSystem` with all primary constants having unit value.

$(textunits(Natural,:Natural))

The well known `Natural` values for `length`, `time`, `mass`, and `charge` are:
```Julia
julia> length(Natural,Metric)
$(length(Natural,Metric))

julia> time(Natural,Metric)
$(time(Natural,Metric))

julia> mass(Natural,Metric)
$(mass(Natural,Metric))

julia> charge(Natural,Metric)
$(charge(Natural,Metric))
```
""" Natural

@doc """
    NaturalGauss = UnitSystem(𝟏,𝟏,𝟏,𝟐*τ,𝟏}

Natural (Gauss) `UnitSystem` with the Gaussian `permeability` value of `4π`.

$(textunits(NaturalGauss,:NaturalGauss))
""" NaturalGauss

@doc """
    QCD = UnitSystem(𝟏,𝟏,𝟏,𝟏,𝟏/μₚₑ)

Qunatum chromodynamics `UnitSystem` with `electronmass` of `𝟏/μₚₑ` or `𝟏/$μₚₑ`.

$(textunits(QCD,:QCD))

The well known `QCD` values for `length`, `time`, `mass`, and `charge` are:
```Julia
julia> length(QCD,Metric) # lQCD
$(length(QCD,Metric))

julia> time(QCD,Metric) # tQCD
$(time(QCD,Metric))

julia> mass(QCD,Metric) # mQCD
$(mass(QCD,Metric))

julia> charge(QCD,Metric)
$(charge(QCD,Metric))
```
""" QCD

@doc """
    QCDGauss = UnitSystem(𝟏,𝟏,𝟏,𝟐*τ,𝟏/μₚₑ)

Qunatum chromodynamics (Gauss) `UnitSystem` with `electronmass` of `𝟏/μₚₑ`.

$(textunits(QCDGauss,:QCDGauss))

The well known `QCDGauss` values for `length`, `time`, `mass`, and `charge` are:
```Julia
julia> length(QCDGauss,Metric) # lQCD
$(length(QCDGauss,Metric))

julia> time(QCDGauss,Metric) # tQCD
$(time(QCDGauss,Metric))

julia> mass(QCDGauss,Metric) # mQCD
$(mass(QCDGauss,Metric))

julia> charge(QCDGauss,Metric)
$(charge(QCDGauss,Metric))
```
""" QCDGauss

@doc """
    QCDoriginal = UnitSystem(𝟏,𝟏,𝟏,𝟐*τ*α,𝟏/μₚₑ)

Qunatum chromodynamics (original) `UnitSystem` with `permeability` of `4π*α`.

$(textunits(QCDoriginal,:QCDoriginal))

The well known `QCDoriginal` values for `length`, `time`, `mass`, and `charge` are:
```Julia
julia> length(QCDoriginal,Metric) # lQCD
$(length(QCDoriginal,Metric))

julia> time(QCDoriginal,Metric) # tQCD
$(time(QCDoriginal,Metric))

julia> mass(QCDoriginal,Metric) # mQCD
$(mass(QCDoriginal,Metric))

julia> charge(QCDoriginal,Metric)
$(charge(QCDoriginal,Metric))
```
""" QCDoriginal
