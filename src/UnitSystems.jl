module UnitSystems

#   This file is part of UnitSystems.jl. It is licensed under the MIT license
#   UnitSystems Copyright (C) 2020 Michael Reed

import Base: @pure, length, time

export slug, ft, KJ1990, KJ2014, RK1990, RK2014, mₑ1990, mₑ2014, temp, units
export slugs, kilograms, lbm, meters, feet, rankine, kelvin, moles, molecules
export UnitSystem, US, SI, CGS, CGS2019, CGSm, CGSe, HLU, FFF

const Systems = (:Metric,:SI2019,:CODATA,:Conventional,:MTS,:English,:EnglishUS,:IAU,:SI1976,:Mixed,:ESU2019,:EMU2019,:EMU,:ESU,:Gauss,:LorentzHeaviside,:Thomson,:Kennelly,:Planck,:PlanckGauss,:Stoney,:Hartree,:Rydberg,:Schrodinger,:Electronic,:Natural,:NaturalGauss,:QCD,:QCDGauss,:QCDoriginal)
const Constants = (:hyperfine,:lightspeed,:planck,:planckreduced,:electronmass,:molarmass,:boltzmann,:permeability,:rationalization,:lorentz,:luminousefficacy)
const Physics = (:atomicmass,:protonmass,:planckmass,:newton,:einstein,:hartree,:rydberg,:bohr,:bohrreduced,:electronradius,:avogadro,:universal,:stefan,:radiationdensity,:permittivity,:coulomb,:ampere,:biotsavart,:charge,:faraday,:impedance,:conductance,:klitzing,:josephson,:magneticflux,:magneton)
const Kinematic = (:time,:length,:area,:volume,:wavenumber,:fuelefficiency,:frequency,:frequencydrift,:speed,:acceleration,:jerk,:snap,:volumeflow)
const Mechanical = (:mass,:massflow,:lineardensity,:areadensity,:density,:specificvolume,:force,:stiffness,:pressure,:compressibility,:viscosity,:diffusivity,:rotationalinertia,:momentum,:angularmomentum,:yank,:energy,:specificenergy,:action,:fluence,:power,:powerdensity,:intensity,:spectralflux,:soundexposure,:impedance,:specificimpedance,:admittance,:compliance,:inertance)
const Electromagnetic = (:charge,:chargedensity,:linearchargedensity,:exposure,:mobility,:current,:currentdensity,:resistance,:conductance,:resistivity,:conductivity,:capacitance,:inductance,:reluctance,:permeance,:permittivity,:permeability,:susceptibility,:specificsusceptibility,:demagnetizingfactor,:vectorpotential,:electricpotential,:magneticpotential,:electricfield,:magneticfield,:electricflux,:magneticflux,:electricfluxdensity,:magneticfluxdensity,:electricdipolemoment,:magneticdipolemoment,:electricpolarizability,:magneticpolarizability,:magneticmoment,:magnetizability,:magnetization,:specificmagnetization,:rigidity,:polestrength)
const Thermodynamic = (:temperature,:entropy,:specificentropy,:volumeheatcapacity,:thermalconductivity,:thermalconductance,:thermalresistance,:thermalexpansion,:lapserate)
const Molar = (:molarmass,:molality,:mole,:molarity,:molarvolume,:molarentropy,:molarenergy,:molarconductivity,:molarsusceptibility,:catalysis,:specificity)
const Photometric = (:luminousflux,:luminance,:luminousenergy,:luminousexposure,:luminousefficacy)
const Mechanics = [Kinematic...,Mechanical...]
const Convert = [Mechanics...,Electromagnetic...,Thermodynamic...,Molar...,Photometric...]

listext(x) = join(x,"`, `")

# unit systems

"""
    UnitSystem{kB,ħ,𝘤,μ₀,mₑ,λ,αL}

Fundamental constants of physics are: `kB` Boltzmann's constant, `ħ` reduced Planck's constant, `𝘤` speed of light, `μ₀` vacuum permeability, `mₑ` electron rest mass, `λ` Gauss rationalization, and `αL` Lorentz's constant.
Primarily the `Metric` SI unit system is used in addition to the historic `English` engineering unit system.
These constants induce derived values for `avogadro`, `boltzmann`, `universal`, `planck`, `planckreduced`, `lightspeed`, `planckmass`, `atomicmass`, `protonmass`, `electronmass`, `newton`, `einstein`, `permeability`, `permittivity`, `coulomb`, and
additional constants `molarmass`, `hyperfine`, `luminousefficacy`, `stefan`, `radiationintensity`, `ampere`, `lorentz`, `biotsavart`, `rationalization`, `impedance`, `charge`, `magneton`, `conductance`, `faraday`, `magneticflux`, `josephson`, `klitzing`, `hartree`, `rydberg`, `bohr`, and `bohrreduced`.

Additional reference `UnitSystem` variants: `EMU`, `ESU`, `Gauss`, `LorentzHeaviside`, `MTS`, `SI2019`, `CODATA`, `Conventional`, `IAU`, `EnglishUS`; and natural atomic units based on gravitational coupling `αG` and the fine structure `1/αinv` constant (`Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal`).

**Derived unit conversions:**

Mechanics: `$(listext(Kinematic))`, `$(listext(Mechanical))`;
Electromagnetics: `$(listext(Electromagnetic))`;
Thermodynamics: `$(listext(Thermodynamic))`,
`$(listext(Molar))`, `$(listext(Photometric))`.
""" #`Rᵤ,mᵤ,σ,ħ,μ₀,ε₀,kₑ,𝘦,𝔉,RK,Z₀,G₀`
struct UnitSystem{kB,ħ,𝘤,μ₀,mₑ,λ,αL} end
@pure UnitSystem{k,ħ,𝘤,μ,m,λ}() where {k,ħ,𝘤,μ,m,λ} = UnitSystem{k,ħ,𝘤,μ,m,λ,1}()
@pure UnitSystem{k,ħ,𝘤,μ,m}() where {k,ħ,𝘤,μ,m} = UnitSystem{k,ħ,𝘤,μ,m,1}()
@pure boltzmann(::UnitSystem{k}) where k = k
@pure planckreduced(::UnitSystem{k,ħ}) where {k,ħ} = ħ
@pure lightspeed(::UnitSystem{k,ħ,𝘤}) where {k,ħ,𝘤} = 𝘤
@pure permeability(::UnitSystem{k,ħ,𝘤,μ}) where {k,ħ,𝘤,μ} = μ
@pure electronmass(::UnitSystem{k,ħ,𝘤,μ,m}) where {k,ħ,𝘤,μ,m} = m
@pure rationalization(::UnitSystem{k,ħ,𝘤,μ,m,λ}) where {k,ħ,𝘤,μ,m,λ} = λ
@pure lorentz(::UnitSystem{k,ħ,𝘤,μ,m,λ,α}) where {k,ħ,𝘤,μ,m,λ,α} = α
# ΔνCs:s⁻¹, c:m⋅s⁻¹, h:kg⋅m²⋅s⁻¹, kB:kg⋅m²⋅s⁻²⋅K⁻¹, NA:mol⁻¹, Kcd: cd⋅sr⋅s³⋅kg⁻¹⋅m⁻²

isrationalized(U::UnitSystem) = rationalization(U) ≠ 4π

Base.display(U::UnitSystem) = println("UnitSystem{kB=$(boltzmann(U)),ħ=$(planckreduced(U)),𝘤=$(lightspeed(U)),μ₀=$(permeability(U)),mᵤ=$(electronmass(U)),λ=$(isrationalized(U) ? rationalization(U) : "4π"),αL=$(lorentz(U))}")

@pure unit(x,y=1) = isapprox(y,x,rtol=eps()^0.9) ? y : x
@pure mass(U::UnitSystem,S::UnitSystem) = electronmass(U,S)
@pure electronmass(𝘩::Float64,R∞::Float64=R∞) = αinv^2*R∞*2𝘩/𝘤
@pure planckmass(U::UnitSystem) = mass(mP,U)
@pure planck(U::UnitSystem) = 2π*planckreduced(U)
@pure newton(U::UnitSystem) = lightspeed(U)*planckreduced(U)/planckmass(U)^2
@pure charge(U::UnitSystem) = sqrt(2planck(U)/impedance(U)/αinv) # fine structure
@pure impedance(U::UnitSystem) = permeability(U)*lightspeed(U)*rationalization(U)*lorentz(U)^2

for unit ∈ Constants
    @eval @pure $unit(U::UnitSystem,S::UnitSystem) = unit($unit(S)/$unit(U))
end
for unit ∈ Convert
    @eval begin
        @pure @inline $unit(v::Real,U::UnitSystem) = $unit(v,U,Metric)
        @pure @inline $unit(v::Real,U::UnitSystem,S::UnitSystem) = (u=$unit(U,S);isone(u) ? v : v/u)
        @pure @inline $unit(v::Real,U::UnitSystem{kB,ħ,𝘤,μ₀,mₑ},S::UnitSystem{kB,ħ,𝘤,μ₀,mₑ}) where {kB,ħ,𝘤,μ₀,mₑ} = v
        @pure @inline $unit(U::UnitSystem{kB,ħ,𝘤,μ₀,mₑ},S::UnitSystem{kB,ħ,𝘤,μ₀,mₑ}) where {kB,ħ,𝘤,μ₀,mₑ} = 1
    end
    if unit ∉ (Constants...,:permittivity,:charge,:magneticflux,:impedance,:conductance)
        @eval @pure @inline $unit(U::UnitSystem) = $unit(U,Metric)
    end
end
for unit ∈ (Systems...,Constants...,Physics...,Convert...)
    @eval export $unit
end

function (U::UnitSystem)(JK=1,Js=1,ms=1,Hm=1,kg=1)
    kB = boltzmann(U)*JK
    ħ = planckreduced(U)*Js
    𝘤 = lightspeed(U)*ms
    μ₀ = permeability(U)*Hm
    mₑ = electronmass(U)*kg
    λ = rationalization(U)
    αL = lorentz(U)
    UnitSystem{kB,ħ,𝘤,μ₀,mₑ,λ,isone(αL) ? αL : αL/ms}()
end

# common conversion factors

const atm,g₀,lbm = 101325.0,9.80665,32.17404856 # lb-f to pdl
const slug,ft,ftUS,rankine,kelvin = 0.45359237lbm,0.3048,1200/3937,5/9,9/5
const kcalₜₕ,kcal₄,kcal₁₀,kcal₂₀,kcalₘ,kcalᵢₜ = 4184,4204,4185.5,4182,4190,4186.8
const calₜₕ,cal₄,cal₁₀,cal₂₀,calₘ,calᵢₜ = (kcalₜₕ,kcal₄,kcal₁₀,kcal₂₀,kcalₘ,kcalᵢₜ)./1e3
const kcal = kcalₜₕ; const cal = kcal/1000 # calₜₕ thermal calorie

# fundamental constants, αinv = (34259-1/4366.8123)/250 # 137.036 exactly?

const ΔνCs,Kcd,mP = 9192631770.0,683.002,2.176434e-8 # planck mass (kg)
const NA,kB,𝘩,𝘤,𝘦 = 6.02214076e23,1.380649e-23,6.62607015e-34,299792458.,1.602176634e-19
const μₑᵤ,μₚᵤ,αinv,R∞ = 1/1822.888486209,1.007276466621,137.035999084,10973731.5681601
const αL,μ₀ = 0.01/𝘤,2𝘩/𝘤/αinv/𝘦^2 # ≈ 4π*(1e-7+5.5e-17), exact charge
const ħ,δμ₀,μₚₑ,Rᵤ,mₑ = 𝘩/2π,μ₀-4π*1e-7,μₚᵤ/μₑᵤ,NA*kB,electronmass(𝘩)
const RK1990,RK2014,KJ1990,KJ2014 = 25812.807,25812.8074555,4.835979e14,4.835978525e14
const ħ1990,ħ2014 = 2/RK1990/KJ1990^2/π,2/RK2014/KJ2014^2/π
const mₑ1990,mₑ2014 = electronmass(2π*ħ1990),electronmass(2π*ħ2014)

# engineering units # Thomson: αL = 1/2

const Gauss = UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,1,1000mₑ,4π,0.01/𝘤}()
const LorentzHeaviside = UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,1,1000mₑ,1,0.01/𝘤}()
const Thomson = UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,1,1000mₑ,4π,1/2}()
const Kennelly = UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,ħ,𝘤,1e-7,mₑ,4π}() # ?
const ESU = UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,(100𝘤)^-2,1000mₑ,4π}()
const ESU2019 = UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e3*μ₀/𝘤^2,1000mₑ}()
const EMU = UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,1,1000mₑ,4π}()
const EMU2019 = UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e7*μ₀,1000mₑ}()
const MTS = UnitSystem{1e6*Rᵤ*mₑ/μₑᵤ,1000ħ,𝘤,4π/1e4,mₑ/1000}()
const Mixed = UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,ħ,𝘤,μ₀,mₑ}()
const Metric = UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,ħ,𝘤,4π*1e-7,mₑ}()
const SI1976 = UnitSystem{8.31432mₑ/μₑᵤ/0.001,ħ,𝘤,4π*1e-7,mₑ}()
const SI2019 = UnitSystem{kB,ħ,𝘤,μ₀,mₑ}()
const CODATA = UnitSystem{Rᵤ*mₑ2014/μₑᵤ/0.001,ħ2014,𝘤,2RK2014/𝘤/αinv,mₑ2014}()
const Conventional = UnitSystem{Rᵤ*mₑ1990/μₑᵤ/0.001,ħ1990,𝘤,2RK1990/𝘤/αinv,mₑ1990}()
const English = UnitSystem{kB*rankine/slug/ft^2,ħ/slug/ft^2,𝘤/ft,4π,mₑ/slug}()
const EnglishUS = UnitSystem{1000Rᵤ*mₑ/μₑᵤ*rankine/slug/ftUS^2,ħ/slug/ftUS^2,𝘤/ftUS,4π,mₑ/slug}()

# astronomical units

const GMsun,GMearth,GMjupiter =  1.32712442099e20,398600441.8e6,1.26686534e17
const au,LD,day = 149597870.7e3,384402e3,60^2*24
const pc,ly,GG = au*648000/π,365.25𝘤*day,newton(SI2019)
const mₛ = GMsun/GG; const Jₛ = mₛ*au^2/day^2; export mₛ,Jₛ,au,day
const IAU = UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001/Jₛ,ħ/day/Jₛ,day*𝘤/au,4π*1e-7*day^2/Jₛ,mₑ/mₛ}()

# aliased & humorous units

const mf = mass(90/lbm,Metric,English); const Jf = mf*(201.168/14day)^2
const FFF = UnitSystem{1000Rᵤ*mₑ/μₑᵤ*rankine/Jf,ħ/14day/Jf,14day*𝘤/201.168,0,mₑ/mf}()
const units, US, SI, temp = UnitSystem, UnitSystem, SI2019, temperature
const CGS, CGS2019, CGSm, CGSe, HLU = Gauss, EMU2019, EMU, ESU, LorentzHeaviside

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
@pure molarmass(U::UnitSystem{kB}) = electronmass(U)*NA/μₑᵤ
@pure molarmass(U::UnitSystem{1e7*kB}) = 1000molarmass(SI2019)
@pure molarmass(U::UnitSystem{1e3*kB}) = molarmass(SI2019)/1000
@pure molarmass(U::UnitSystem{kB}) where kB = molarmass(CGS)/1000
@pure molarmass(U::UnitSystem{boltzmann(MTS)}) = molarmass(CGS)/1e6
@pure molarmass(U::UnitSystem{boltzmann(CGS)}) = molarmass(Natural)
@pure molarmass(U::UnitSystem{boltzmann(FFF)}) = molarmass(Natural)
@pure molarmass(U::UnitSystem{boltzmann(English)}) = 1000molarmass(SI2019)
@pure molarmass(U::UnitSystem{boltzmann(EnglishUS)}) = molarmass(Natural)
@pure molarmass(U::UnitSystem{boltzmann(IAU)}) = 1/1000mₛ

@pure luminousefficacy(U::UnitSystem{1}) = 1
@pure luminousefficacy(U::UnitSystem) = power(Kcd,SI2019,U)

include("kinematic.jl")
include("electromagnetic.jl")
include("thermodynamic.jl")
include("physics.jl")

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

end # module
