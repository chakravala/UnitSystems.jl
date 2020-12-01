module UnitSystems

#   This file is part of UnitSystems.jl. It is licensed under the MIT license
#   UnitSystems Copyright (C) 2020 Michael Reed

import Base: @pure

export slug, ft, KJ1990, KJ2014, RK1990, RK2014, mₑ2014
export mass, slugs, kilograms, poundal, meters, feet, rankine, kelvin, moles, molecules
export UnitSystem, CGS, CGS2019, Metric, SI2019, CODATA, Conventional, English, IAU
export Planck, PlanckGauss, Stoney, Hartree, Rydberg, Schrodinger, Electronic, Natural, NaturalGauss, QCD, QCDGauss, QCDoriginal

# unit systems

"""
    UnitSystem{kB,ħ,𝘤,μ₀,mₑ}

Standardized for engineering based on fundamental constants: `kB` Boltzmann's constant, `ħ` reduced Planck's constant, `𝘤` speed of light, `μ₀` vacuum permeability, and `mₑ` electron rest mass.
Primarily the `Metric` SI unit system is used in addition to the historic `English` engineering unit system.
These constants induce derived values for `avogadro`, `boltzmann`, `universal`, `planck`, `planckreduced`, `lightspeed`, `planckmass`, `atomicmass`, `protonmass`, `electronmass`, `newton`, `einstein`, `permeability`, `permittivity`, `coulomb`, and
additional constants `stefan`, `radiationintensity`, `impedance`, `charge`, `magneton`, `conductance`, `faraday`, `magneticflux`, `josephson`, `klitzing`, `hartree`, `rydberg`, `bohr`, `bohrreduced`, and `molarmass`.

Additional reference `UnitSystem` variants `CGS`, `CGS2019`, `SI2019`, `CODATA`, `Conventional`, `IAU`; along with several natural atomic units based on the fine structure constant `1/αinv` and the gravitational coupling constant `αG` (`Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal`).
""" #`Rᵤ,mᵤ,σ,ħ,μ₀,ε₀,kₑ,𝘦,𝔉,RK,Z₀,G₀`
struct UnitSystem{kB,ħ,𝘤,μ,mₑ} end
@pure boltzmann(::UnitSystem{k}) where k = k
@pure planckreduced(::UnitSystem{k,h}) where {k,h} = h
@pure lightspeed(::UnitSystem{k,h,c}) where {k,h,c} = c
@pure permeability(::UnitSystem{k,h,c,μ}) where {k,h,c,μ} = μ
@pure electronmass(::UnitSystem{k,h,c,μ,m}) where {k,h,c,μ,m} = m
# ΔνCs:s⁻¹, c:m⋅s⁻¹, h:kg⋅m²⋅s⁻¹, kB:kg⋅m²⋅s⁻²⋅K⁻¹, NA:mol⁻¹, Kcd: cd⋅sr⋅s³⋅kg⁻¹⋅m⁻²

@pure electronmass(R∞::Float64=R∞,𝘩::Float64=𝘩) = αinv^2*R∞*2𝘩/𝘤
@pure mass(U::UnitSystem,S::UnitSystem) = electronmass(S)/electronmass(U)
@pure newton(U::UnitSystem) = lightspeed(U)*planckreduced(U)/planckmass(U)^2
@pure planckmass(U::UnitSystem) = mass(mP,U)
@pure planck(U::UnitSystem) = 2π*planckreduced(U)
@pure impedance(U::UnitSystem) = permeability(U)*lightspeed(U)
@pure charge(U::UnitSystem) = sqrt(2planck(U)/impedance(U)/αinv) # fine structure
@pure charge(U::UnitSystem,S::UnitSystem) = charge(S)/charge(U)

for unit ∈ (:mass,:length,:time,:temperature,:charge)
    @eval @pure $unit(v::Real,U::UnitSystem,S::UnitSystem=Metric) = v/$unit(U,S)
    unit ≠ :charge && @eval @pure $unit(U::UnitSystem) = $unit(U,Metric)
end

Base.display(U::UnitSystem) = println("UnitSystem{kB=$(boltzmann(U)),ħ=$(planckreduced(U)),𝘤=$(lightspeed(U)),μ₀=$(permeability(U)),mᵤ=$(electronmass(U))}")

# common conversion factors

const atm,g₀,lbm = 101325.0,9.80665,32.17404856 # lb-f to pdl
const slug,ft,ftUS,rankine,kelvin = 0.45359237lbm,0.3048,1200/3937,9/5,5/9
const kcalₜₕ,kcal₄,kcal₁₀,kcal₂₀,kcalₘ,kcalᵢₜ = 4184,4204,4185.5,4182,4190,4186.8
const calₜₕ,cal₄,cal₁₀,cal₂₀,calₘ,calᵢₜ = (kcalₜₕ,kcal₄,kcal₁₀,kcal₂₀,kcalₘ,kcalᵢₜ)./1e3
const kcal = kcalₜₕ; const cal = kcal/1000 # calₜₕ thermal calorie

# fundamental constants, αinv = (34259-1/4366.8123)/250 # 137.036 exactly?

const ΔνCs,Kcd,mP = 9192631770.0,683.0,2.176434e-8 # planck mass (kg)
const NA,kB,𝘩,𝘤,𝘦 = 6.02214076e23,1.380649e-23,6.62607015e-34,299792458.,1.602176634e-19
const μₑᵤ,μₚᵤ,αinv,R∞ = 1/1822.888486209,1.007276466621,137.035999084,10973731.5681601
const μ₀ = 2𝘩/𝘤/αinv/𝘦^2 # ≈ 4π*(1e-7+5.5e-17), exact charge
const ħ,δμ₀,μₚₑ,Rᵤ,mₑ = 𝘩/2π,μ₀-4π*1e-7,μₚᵤ/μₑᵤ,NA*kB,electronmass(R∞,𝘩)
const RK1990,RK2014,KJ1990,KJ2014 = 25812.807,25812.8074555,4.835979e14,4.835978525e14
const ħ2014 = 2/RK2014/KJ2014^2/π; const mₑ2014 = electronmass(10973731.568508,ħ2014)

# engineering units

const CGS = UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,4π,1000mₑ}()
const CGS2019 = UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e7*μ₀,1000mₑ}()
const Metric = UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,ħ,𝘤,4π*1e-7,mₑ}()
const SI1976 = UnitSystem{8.31432mₑ/μₑᵤ/0.001,ħ,𝘤,4π*1e-7,mₑ}()
const SI2019 = UnitSystem{kB,ħ,𝘤,μ₀,mₑ}()
const CODATA = UnitSystem{Rᵤ*mₑ2014/μₑᵤ/0.001,ħ2014,𝘤,2RK2014/𝘤/αinv,mₑ2014}()
const Conventional = UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,2/RK1990/KJ1990^2/π,𝘤,2RK1990/𝘤/αinv,mₑ}()
const English = UnitSystem{5.657302466e-24,ħ/slug/ft^2,𝘤/ft,4π,mₑ/slug}()
const EnglishNew = UnitSystem{1000rankine/slug/ft*Rᵤ*mₑ/μₑᵤ,ħ/slug/ft^2,𝘤/ft,4π,mₑ/slug}()
const Sudgy = UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,ħ,𝘤,μ₀,mₑ}()

# astronomical units

const GMsun,GMearth,GMjupiter =  1.32712442099e20,398600441.8e6,1.26686534e17
const au,LD,day = 149597870.7e3,384402e3,60^2*24
const pc,ly,GG = au*648000/π,365.25𝘤*day,newton(SI2019)
const mₛ = GMsun/GG; const Jₛ = mₛ*au^2/day^2; export mₛ,Jₛ,au,day
const IAU = UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001/Jₛ,ħ*day/Jₛ,day*𝘤/au,4π*1e-7*day^2/Jₛ,mₑ/mₛ}()

# natural units

const αG = (mₑ/mP)^2
const Planck = UnitSystem{1,1,1,1,√(4π*αG)}()
const PlanckGauss = UnitSystem{1,1,1,4π,√αG}()
const Stoney = UnitSystem{1,αinv,1,4π,√(αG*αinv)}()
const Hartree = UnitSystem{1,1,αinv,4π/αinv^2,1}()
const Rydberg = UnitSystem{1,1,2αinv,π/αinv^2,1/2}()
const Schrodinger = UnitSystem{1,1,αinv,4π/αinv^2,√(αG*αinv)}()
const Electronic = UnitSystem{1,αinv,1,4π,1}()
const Natural = UnitSystem{1,1,1,1,1}()
const NaturalGauss = UnitSystem{1,1,1,4π,1}()
const QCD = UnitSystem{1,1,1,1,1/μₚₑ}()
const QCDGauss = UnitSystem{1,1,1,4π,1/μₚₑ}()
const QCDoriginal = UnitSystem{1,1,1,4π/αinv,1/μₚₑ}()

# physical constants

@pure molarmass(U::UnitSystem{1}) = 1
@pure molarmass(U::UnitSystem{boltzmann(CGS)}) = molarmass(Natural)
@pure molarmass(U::UnitSystem{kB}) where kB = molarmass(CGS)/1000
@pure molarmass(U::UnitSystem{1e7*kB}) = 1000molarmass(SI2019)
@pure molarmass(U::UnitSystem{kB}) = electronmass(U)*NA/μₑᵤ
@pure molarmass(U::UnitSystem{boltzmann(IAU)}) = 1/1000mₛ
@pure molarmass(U::UnitSystem{boltzmann(English)}) = Rᵤ*mₑ/μₑᵤ/slug/ft*rankine/boltzmann(English)

include("physics.jl")
include("convert.jl")

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
const ℓP = plancklength(SI2019)
const RK = klitzing(SI2019) #
const KJ = josephson(SI2019) #
const RH,Ry = R∞*mₚ/(mₑ+mₚ),𝘩*𝘤*R∞
const Mu,Ru,SB,hh,cc,m0,e0,ke,me,mp,mu,ee,FF,Z0,G0,Eh,a0,re,g0 = Mᵤ,Rᵤ,σ,𝘩,𝘤,μ₀,ε₀,kₑ,mₑ,mₚ,mᵤ,𝘦,𝔉,Z₀,G₀,Eₕ,a₀,rₑ,g₀
export κ, GG, NA, kB, Rᵤ, σ, 𝘩, ħ, 𝘤, μ₀, ε₀, kₑ, mₑ, mₚ, mᵤ, 𝘦, 𝔉, Φ₀, Z₀, G₀, Eₕ, R∞, a₀, rₑ, KJ, RK, Ru, SB, hh, cc, m0, e0, ke, me, mp, mu, ee, FF, Z0, G0, Eh, a0, re
export αG, αinv, μₚₑ, μₑᵤ, μₚᵤ, mpe, meu, mpu, mP, δμ₀, Mᵤ, Mu, RH, Ry, ΔνCs, Kcd, lbm
export cal, kcal, calₜₕ, kcalₜₕ, calᵢₜ, kcalᵢₜ, SI, SI1976, ℓP, plancklength, g₀, g0, atm
const mpe, mea, mpu, SI = μₚₑ, μₑᵤ, μₚᵤ, SI2019

const Constants = (:newton,:avogadro,:boltzmann,:planck,:planckreduced,:lightspeed,:universal,:permeability,:permittivity,:coulomb)
const Physics = (:electronmass,:protonmass,:atomicmass,:planckmass,:stefan,:radiationdensity,:einstein,:impedance,:charge,:faraday,:josephson,:klitzing,:hartree,:rydberg,:bohr,:bohrreduced,:electronradius,:conductance,:magneticflux,:magneton,:molarmass)

for op ∈ (Constants...,Physics...)
    @eval export $op
end

# engineering unit systems docs

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
    English::UnitSystem{5.657302466e-24,ħ/slug/ft^2,𝘤/ft,4π,mₑ/slug}

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
    CGS::UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,4π,1000mₑ}

Centimetre-gram-second `UnitSystem` variant of `Metric` system based on factors of `1e2,1e3`.

```Julia
julia> boltzmann(CGS) # erg⋅K⁻¹
$(boltzmann(CGS))

julia> planckreduced(CGS) # erg⋅s⋅rad⁻¹
$(planckreduced(CGS))

julia> lightspeed(CGS) # cm⋅s⁻¹
$(lightspeed(Metric))

julia> permeability(CGS) # erg⋅A⁻²⋅cm⁻¹
$(permeability(CGS))

julia> electronmass(CGS) # g
$(electronmass(CGS))
```
""" CGS

@doc """
    CGS2019::UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e7*μ₀,1000mₑ}

Centimetre-gram-second `UnitSystem` variant of the tuned `SI2019` unit specification.

```Julia
julia> boltzmann(CGS2019) # erg⋅K⁻¹
$(boltzmann(CGS2019))

julia> planckreduced(CGS2019) # erg⋅s⋅rad⁻¹
$(planckreduced(CGS2019))

julia> lightspeed(CGS2019) # cm⋅s⁻¹
$(lightspeed(CGS2019))

julia> permeability(CGS2019) # erg⋅A⁻²⋅cm⁻¹
$(permeability(CGS2019))

julia> electronmass(CGS2019 # g
$(electronmass(CGS2019))
```
""" CGS2019

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
$(permeability(CODATA))

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
    Conventional::UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,2/RK1990/KJ1990^2/π,𝘤,2RK1990/𝘤/αinv,mₑ}

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

# other unit system docs

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
""" PlanckGauss

@doc """
    Stoney::UnitSystem{1,αinv,1,4π,√(αG*αinv)}

Stoney `UnitSystem` with `permeability` of `4π` and `electronmass` coupling `√(αG*αinv)`.

$(textunits(Stoney,:Stoney))
""" Stoney

@doc """
    Hartree::UnitSystem{1,1,αinv,4π/αinv^2,1}

Hartree atomic `UnitSystem` with `lightspeed` of `αinv` and `permeability` of `4π/αinv^2`.

$(textunits(Hartree,:Hartree))
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

const Electronic = UnitSystem{1,αinv,1,4π,1}()
@doc """
    Electronic::UnitSystem{1,αinv,1,4π,1}

Electronic `UnitSystem` with `planckreduced` of `αinv` and `permeability` of `4π`.

$(textunits(Electronic,:Electronic))
""" Electronic

@doc """
    Natural::UnitSystem{1,1,1,1,1}

Natural `UnitSystem` with all primary constants having unit value.

$(textunits(Natural,:Natural))
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
""" QCD

@doc """
    QCDGauss::UnitSystem{1,1,1,4π,1/μₚₑ}

Qunatum chromodynamics (Gauss) `UnitSystem` with `electronmass` of `1/μₚₑ`.

$(textunits(QCDGauss,:QCDGauss))
""" QCDGauss

@doc """
    QCDoriginal::UnitSystem{1,1,1,4π/αinv,1/μₚₑ}

Qunatum chromodynamics (original) `UnitSystem` with `permeability` of `4π/αinv`.

$(textunits(QCDoriginal,:QCDoriginal))
""" QCDoriginal

end # module
