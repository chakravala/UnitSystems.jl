
#   This file is part of UnitSystems.jl. It is licensed under the MIT license
#   UnitSystems Copyright (C) 2021 Michael Reed

export Universe, coupling, finestructure, electronunit, protonunit, protonelectron
#const Mu,Ru,SB,hh,cc,m0,e0,ke,me,mp,mu,ee,FF,Z0,G0,Eh,a0,re,g0,lP,ϵ₀,mB = Mᵤ,Rᵤ,σ,𝘩,𝘤,μ₀,ε₀,kₑ,mₑ,mₚ,mᵤ,𝘦,𝔉,Z₀,G₀,Eₕ,a₀,rₑ,g₀,ℓP,ε₀,μB
export slug, ft, KJ1990, KJ2014, RK1990, RK2014, mₑ1990, mₑ2014, temp, units
export slugs, kilograms, lbm, meters, feet, rankine, kelvin, moles, molecules
export UnitSystem, US, SI, MKS, CGS, CGS2019, CGSm, CGSe, HLU, FFF

# == Metric is different
const κ = einstein(SI2019)
const σ = stefan(SI2019) #
const μB = magneton(SI2019) #
const ε₀ = permittivity(SI2019) #
const kₑ = coulomb(SI2019) #
const mₚ = protonmass(SI2019)
const mᵤ = atomicmass(SI2019)
const Mᵤ = molarmass(SI2019)
const 𝔉 = faraday(SI2019) #
const Φ₀ = magneticflux(SI2019) #
const Z₀ = impedance(SI2019) #
const G₀ = conductance(SI2019) #
const Eₕ = hartree(SI2019)
const a₀ = bohr(SI2019)
const rₑ = electronradius(SI2019)
const RK = klitzing(SI2019) #
const KJ = josephson(SI2019) #
const RH,Ry = R∞*mₚ/(mₑ+mₚ),𝘩*𝘤*R∞

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

const BTUftlb = 3600/0.5778thermalconductivity(English) # BTU⋅ft⁻¹⋅lb⁻¹
const BTUJ = energy(English)*BTUftlb # BTU⋅J⁻¹

# constant aliases

const mpe, meu, mpu, ainv, aG = μₚₑ, μₑᵤ, μₚᵤ, αinv, αG
const Mu,Ru,SB,hh,cc,m0,e0,ke,me,mp,mu,ee,FF,Z0,G0,Eh,a0,re,g0,lP,aL,ϵ₀ = Mᵤ,Rᵤ,σ,𝘩,𝘤,μ₀,ε₀,kₑ,mₑ,mₚ,mᵤ,𝘦,𝔉,Z₀,G₀,Eₕ,a₀,rₑ,g₀,ℓP,αL,ε₀
export κ, GG, NA, kB, Rᵤ, σ, 𝘩, ħ, 𝘤, μ₀, ε₀, kₑ, mₑ, mₚ, mᵤ, 𝘦, 𝔉, Φ₀, Z₀, G₀, Eₕ, R∞, a₀, rₑ, KJ, RK, Ru, SB, hh, cc, m0, e0, ke, me, mp, mu, ee, FF, Z0, G0, Eh, a0, re, μB
export αG, αinv, μₚₑ, μₑᵤ, μₚᵤ, mpe, meu, mpu, mP, δμ₀, Mᵤ, Mu, RH, Ry, ΔνCs, Kcd, ainv
export cal, kcal, calₜₕ, kcalₜₕ, calᵢₜ, kcalᵢₜ, ℓP, g₀, g0, atm, lbm, BTUJ, BTUftlb, aG
export lP, tP, TP, lS, tS, mS, qS, lA, tA, mA, qA, lQCD, tQCD, mQCD, ϵ₀, αL, aL

# engineering unit systems docs

cgstext(US,AMP,cgs=eval(US)) = """
```Julia
julia> boltzmann($US) # erg⋅K⁻¹
$(boltzmann(cgs))

julia> planckreduced($US) # erg⋅s⋅rad⁻¹
$(planckreduced(cgs))

julia> lightspeed($US) # cm⋅s⁻¹
$(lightspeed(cgs))

julia> permeability($US) # statH⋅cm⁻¹
$(permeability(cgs))

julia> electronmass($US) # g
$(electronmass(cgs))

julia> rationalization($US)
$(rationalization(cgs))
```
"""

for U ∈ (:CGSm,:CGSe,:EMU,:ESU)
    (EU,AMP) = QuoteNode.(U ∉ (:CGSe,:ESU) ? (:EMU,:Bi) : (:ESU,:statA))
@eval @doc """
    $($(QuoteNode(U)))::UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,$($EU≠:EMU ? "(100𝘤)^-2" : 1),1000mₑ,4π}

Centimetre-gram-second `UnitSystem` variant based on `$($EU)` (non-rationalized).

$(cgstext($(QuoteNode(U)),$AMP))
""" $U

U ∉ (:CGSm,:CGSe) && @eval @doc """
    $(Symbol($(QuoteNode(U)),:2019))::UnitSystem{1e7*kB,1e7*ħ,100𝘤,$($EU≠:EMU ? "1e3*μ₀/𝘤^2" : "1e7*μ₀"),1000mₑ}

Centimetre-gram-second `UnitSystem` variant of tuned `SI2019` based on `$($EU)` (rationalized).

$(cgstext(Symbol($(QuoteNode(U)),:2019),$AMP))
""" $(Symbol(U,:2019))
end

@doc """
    Thomson::UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,1,1000mₑ,4π,1/2}

Centimetre-gram-second `UnitSystem` variant `Thomson` (EMU-Lorentz, non-rationalized).

```Julia
julia> boltzmann(Thomson) # erg⋅K⁻¹
$(boltzmann(Thomson))

julia> planckreduced(Thomson) # erg⋅s⋅rad⁻¹
$(planckreduced(Thomson))

julia> lightspeed(Thomson) # cm⋅s⁻¹
$(lightspeed(Thomson))

julia> permeability(Thomson) # abH⋅cm⁻¹
$(permeability(Thomson))

julia> electronmass(Thomson) # g
$(electronmass(Thomson))

julia> rationalization(Thomson)
$(rationalization(Thomson))

julia> lorentz(Thomson)
$(lorentz(Thomson))
```
""" Thomson

@doc """
    Gauss::UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,1,1000mₑ,4π,0.01/𝘤}

Centimetre-gram-second `UnitSystem` variant `CGS` (Gauss-Lorentz, non-rationalized).

```Julia
julia> boltzmann(Gauss) # erg⋅K⁻¹
$(boltzmann(Gauss))

julia> planckreduced(Gauss) # erg⋅s⋅rad⁻¹
$(planckreduced(Gauss))

julia> lightspeed(Gauss) # cm⋅s⁻¹
$(lightspeed(Gauss))

julia> permeability(Gauss) # statH⋅cm⁻¹
$(permeability(Gauss))

julia> electronmass(Gauss) # g
$(electronmass(Gauss))

julia> rationalization(Gauss)
$(rationalization(Gauss))

julia> lorentz(Gauss)
$(lorentz(Gauss))
```
""" Gauss, CGS

@doc """
    LorentzHeaviside::UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,1,1000mₑ,1,0.01/𝘤}

Centimetre-gram-second `UnitSystem` variant `HLU` (Heaviside-Lorentz, rationalized).

```Julia
julia> boltzmann(LorentzHeaviside) # erg⋅K⁻¹
$(boltzmann(LorentzHeaviside))

julia> planckreduced(LorentzHeaviside) # erg⋅s⋅rad⁻¹
$(planckreduced(LorentzHeaviside))

julia> lightspeed(LorentzHeaviside) # cm⋅s⁻¹
$(lightspeed(LorentzHeaviside))

julia> permeability(HLU) # hlH⋅cm⁻¹
$(permeability(LorentzHeaviside))

julia> electronmass(LorentzHeaviside) # g
$(electronmass(LorentzHeaviside))

julia> rationalization(LorentzHeaviside)
$(rationalization(LorentzHeaviside))

julia> lorentz(LorentzHeaviside)
$(lorentz(LorentzHeaviside))
```
""" LorentzHeaviside, HLU

@doc """
    MTS::UnitSystem{1e6*Rᵤ*mₑ/μₑᵤ,1000ħ,𝘤,4π/1e4,mₑ/1000}

Metre-tonne-second `UnitSystem` variant of `Metric` system.

```Julia
julia> boltzmann(MTS) # kJ⋅K⁻¹
$(boltzmann(MTS))

julia> planckreduced(MTS) # kJ⋅s⋅rad⁻¹
$(planckreduced(MTS))

julia> lightspeed(MTS) # m⋅s⁻¹
$(lightspeed(MTS))

julia> permeability(MTS) # kH⋅m⁻¹
$(permeability(MTS))

julia> electronmass(MTS) # t
$(electronmass(MTS))
```
""" MTS

@doc """
    Metric::UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,ħ,𝘤,4π*1e-7,mₑ}

Systeme International d'Unites (the SI units) adopted as the preferred `UnitSystem`.

```Julia
julia> boltzmann(Metric) # J⋅K⁻¹
$(boltzmann(Metric))

julia> planckreduced(Metric) # J⋅s⋅rad⁻¹
$(planckreduced(Metric))

julia> lightspeed(Metric) # m⋅s⁻¹
$(lightspeed(Metric))

julia> permeability(Metric) # H⋅m⁻¹
$(permeability(Metric))

julia> electronmass(Metric) # kg
$(electronmass(Metric))
```
""" Metric

@doc """
    SI2019::UnitSystem{kB,ħ,𝘤,μ₀,mₑ}

Systeme International d'Unites (the SI units) with `μ₀` for a tuned `charge` exactly.

```Julia
julia> boltzmann(SI2019) # J⋅K⁻¹
$(boltzmann(SI2019))

julia> planckreduced(SI2019) # J⋅s⋅rad⁻¹
$(planckreduced(SI2019))

julia> lightspeed(SI2019) # m⋅s⁻¹
$(lightspeed(SI2019))

julia> permeability(SI2019) # H⋅m⁻¹
$(permeability(SI2019))

julia> electronmass(SI2019) # kg
$(electronmass(SI2019))
```
""" SI2019, SI

@doc """
    CODATA::UnitSystem{Rᵤ*mₑ2014/μₑᵤ/0.001,2/RK2014/KJ2014^2/π,𝘤,2RK2014/𝘤/αinv,mₑ2014}

Metric `UnitSystem` based on Committee on Data of the International Science Council.

```Julia
julia> boltzmann(CODATA) # J⋅K⁻¹
$(boltzmann(CODATA))

julia> planckreduced(CODATA) # J⋅s⋅rad⁻¹
$(planckreduced(CODATA))

julia> lightspeed(CODATA) # m⋅s⁻¹
$(lightspeed(CODATA))

julia> permeability(CODATA) # H⋅m⁻¹
$(permeability(CODATA))

julia> electronmass(CODATA) # kg
$(electronmass(CODATA))
```
""" CODATA

@doc """
    Conventional::UnitSystem{Rᵤ*mₑ1990/μₑᵤ/0.001,2/RK1990/KJ1990^2/π,𝘤,2RK1990/𝘤/αinv,mₑ1990}

Conventional electronic `UnitSystem` with 1990 tuned `josephson` and `klitzing` constants.

```Julia
julia> boltzmann(Conventional) # J⋅K⁻¹
$(boltzmann(Conventional))

julia> planckreduced(Conventional) # J⋅s⋅rad⁻¹
$(planckreduced(Conventional))

julia> lightspeed(Conventional) # m⋅s⁻¹
$(lightspeed(Conventional))

julia> permeability(Conventional) # H⋅m⁻¹
$(permeability(Conventional))

julia> electronmass(Conventional) # kg
$(electronmass(Conventional))
```
""" Conventional

@doc """
    IAU::UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001/Jₛ,ħ*day/Jₛ,day*𝘤/au,4π*1e-7*day^2/Jₛ,mₑ/mₛ}

Astronomical (solar) `UnitSystem` defined by International Astronomical Union.

```Julia
julia> boltzmann(IAU) # M⊙⋅au²⋅D⁻²⋅K⁻¹
$(boltzmann(IAU))

julia> planckreduced(IAU) # M⊙⋅au²⋅D⁻¹⋅rad⁻¹
$(planckreduced(IAU))

julia> lightspeed(IAU) # au⋅D⁻¹
$(lightspeed(IAU))

julia> permeability(IAU) # M⊙⋅au²⋅C⁻²
$(permeability(IAU))

julia> electronmass(IAU) # M⊙
$(electronmass(IAU))
```
""" IAU

@doc """
    English::UnitSystem{kB*rankine/slug/ft^2,ħ/slug/ft^2,𝘤/ft,4π,mₑ/slug}

Engineering `UnitSystem` historically used by Britain and United States.

```Julia
julia> boltzmann(English) # ft⋅lb⋅°R⁻¹
$(boltzmann(English))

julia> planckreduced(English) # ft⋅lb⋅s⋅rad⁻¹
$(planckreduced(English))

julia> lightspeed(English) # ft⋅s⁻¹
$(lightspeed(English))

julia> permeability(English) # slug⋅ft²⋅?⁻²
$(permeability(English))

julia> electronmass(English) # slugs
$(electronmass(English))
```
""" English

@doc """
    EnglishUS::UnitSystem{1000Rᵤ*mₑ/μₑᵤ*rankine/slug/ftUS^2,ħ/slug/ftUS^2,𝘤/ftUS,4π,mₑ/slug}

Engineering `UnitSystem` based on the geophysical US survey foot (1200/3937).

```Julia
julia> boltzmann(EnglishUS) # ftUS⋅lb⋅°R⁻¹
$(boltzmann(EnglishUS))

julia> planckreduced(EnglishUS) # ftUS⋅lb⋅s⋅rad⁻¹
$(planckreduced(EnglishUS))

julia> lightspeed(EnglishUS) # ftUS⋅s⁻¹
$(lightspeed(EnglishUS))

julia> permeability(EnglishUS) # slug⋅ftUS²⋅?⁻²
$(permeability(EnglishUS))

julia> electronmass(EnglishUS) # slugs
$(electronmass(EnglishUS))
```
""" EnglishUS

@doc """
    FFF::UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001*rankine/Jf,ħ/14day/Jf,14day*𝘤/201.168,0,mₑ/mf}

Furlong–firkin–fortnight `FFF` is a humorous `UnitSystem` based on unusal impractical units.

```Julia
julia> boltzmann(FFF) # fir⋅fur²⋅ftn⁻²⋅F⁻¹
$(boltzmann(FFF))

julia> planckreduced(FFF) # fir⋅fur²⋅ftn⁻¹⋅rad⁻¹
$(planckreduced(FFF))

julia> lightspeed(FFF) # fur⋅ftn⁻¹
$(lightspeed(FFF))

julia> permeability(FFF) # fir⋅fur²⋅Inf⁻²
$(permeability(FFF))

julia> electronmass(FFF) # fir
$(electronmass(FFF))
```
""" FFF

@doc """
    Kennelly::UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,ħ,𝘤,1e-7,mₑ,4π}

Kennelly ? variant `UnitSystem` of the standard `Metric` units ???

```Julia
julia> boltzmann(Kennelly) # J⋅K⁻¹
$(boltzmann(Kennelly))

julia> planckreduced(Kennelly) # J⋅s⋅rad⁻¹
$(planckreduced(Kennelly))

julia> lightspeed(Kennelly) # m⋅s⁻¹
$(lightspeed(Kennelly))

julia> permeability(Kennelly) # H⋅m⁻¹
$(permeability(Kennelly))

julia> electronmass(Kennelly) # kg
$(electronmass(Kennelly))

julia> rationalization(Kennelly)
$(rationalization(Kennelly))
```
""" Kennelly

# natural unit system docs

textunits(U,S) = """
```Julia
julia> boltzmann($S)
$(boltzmann(U))

julia> planckreduced($S)
$(planckreduced(U))

julia> lightspeed($S)
$(lightspeed(U))

julia> permeability($S)
$(permeability(U))

julia> electronmass($S)
$(electronmass(U))
```
"""

@doc """
    Planck::UnitSystem{1,1,1,1,√(4π*αG)}

Planck `UnitSystem` with the `electronmass` value `√(4π*αG)` using gravitational coupling.

$(textunits(Planck,:Planck))
""" Planck

@doc """
    PlanckGauss::UnitSystem{1,1,1,4π,√αG}

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
    Stoney::UnitSystem{1,αinv,1,4π,√(αG*αinv)}

Stoney `UnitSystem` with `permeability` of `4π` and `electronmass` coupling `√(αG*αinv)`.

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
    Hartree::UnitSystem{1,1,αinv,4π/αinv^2,1}

Hartree atomic `UnitSystem` with `lightspeed` of `αinv` and `permeability` of `4π/αinv^2`.

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
    Rydberg::UnitSystem{1,1,2αinv,π/αinv^2,1/2}

Rydberg `UnitSystem` with `lightspeed` of `2αinv` and `permeability` of `π/αinv^2`.

$(textunits(Rydberg,:Rydberg))
""" Rydberg

@doc """
    Schrodinger::UnitSystem{1,1,αinv,4π/αinv^2,√(αG*αinv)}

Schrodinger `UnitSystem` with `permeability` of `4π/αinv^2` and `electronmass` of `√(αG*αinv)`.

$(textunits(Schrodinger,:Schrodinger))
""" Schrodinger

@doc """
    Electronic::UnitSystem{1,αinv,1,4π,1}

Electronic `UnitSystem` with `planckreduced` of `αinv` and `permeability` of `4π`.

$(textunits(Electronic,:Electronic))
""" Electronic

@doc """
    Natural::UnitSystem{1,1,1,1,1}

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
    NaturalGauss::UnitSystem{1,1,1,4π,1}

Natural (Gauss) `UnitSystem` with the Gaussian `permeability` value of `4π`.

$(textunits(NaturalGauss,:NaturalGauss))
""" NaturalGauss

@doc """
    QCD::UnitSystem{1,1,1,1,1/μₚₑ}

Qunatum chromodynamics `UnitSystem` with `electronmass` of `1/μₚₑ` or `1/$μₚₑ`.

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
    QCDGauss::UnitSystem{1,1,1,4π,1/μₚₑ}

Qunatum chromodynamics (Gauss) `UnitSystem` with `electronmass` of `1/μₚₑ`.

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
    QCDoriginal::UnitSystem{1,1,1,4π/αinv,1/μₚₑ}

Qunatum chromodynamics (original) `UnitSystem` with `permeability` of `4π/αinv`.

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
