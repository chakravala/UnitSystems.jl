module UnitSystems

#   This file is part of UnitSystems.jl. It is licensed under the MIT license
#   UnitSystems Copyright (C) 2020 Michael Reed

import Base: @pure, length, time

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

@pure measure(x) = x

# universe

"""
    Coupling{αG,α,μₑᵤ,μₚᵤ}

Specification of `Universe` with the dimensionless `Coupling` constants `coupling`, `finestructure`, `electronunit`, `protonunit`, and `protonelectron`.
Alterations to these values can be facilitated and quantified using parametric polymorphism.
Due to the `Coupling` interoperability, the `MeasureSystems` package is made possible to support calculations with `Measurements` having error standard deviations.
"""
struct Coupling{αG,α,μₑᵤ,μₚᵤ} end
@pure coupling(U::Coupling{αG}) where αG = measure(αG)
@pure finestructure(U::Coupling{αG,α}) where {αG,α} = measure(α)
@pure electronunit(U::Coupling{αG,α,μₑᵤ}) where {αG,α,μₑᵤ} = measure(μₑᵤ)
@pure protonunit(U::Coupling{αG,α,μₑᵤ,μₚᵤ}) where {αG,α,μₑᵤ,μₚᵤ} = measure(μₚᵤ)
@pure protonelectron(U::Coupling) = protonunit(U)/electronunit(U)
Base.display(U::Coupling) = println("Coupling{αG=$(coupling(U)),α=$(finestructure(U)),μₑᵤ=$(electronunit(U)),μₚᵤ=$(protonunit(U))}")

# unit systems

"""
    UnitSystem{kB,ħ,𝘤,μ₀,mₑ,λ,αL}

Fundamental constants of physics are: `kB` Boltzmann's constant, `ħ` reduced Planck's constant, `𝘤` speed of light, `μ₀` vacuum permeability, `mₑ` electron rest mass, `λ` Gauss rationalization, and `αL` Lorentz's constant.
Primarily the `Metric` SI unit system is used in addition to the historic `English` engineering unit system.
These constants induce derived values for `avogadro`, `boltzmann`, `universal`, `planck`, `planckreduced`, `lightspeed`, `planckmass`, `atomicmass`, `protonmass`, `electronmass`, `newton`, `einstein`, `permeability`, `permittivity`, `coulomb`, and
additional constants `molarmass`, `hyperfine`, `luminousefficacy`, `stefan`, `radiationdensity`, `ampere`, `lorentz`, `biotsavart`, `rationalization`, `impedance`, `charge`, `magneton`, `conductance`, `faraday`, `magneticflux`, `josephson`, `klitzing`, `hartree`, `rydberg`, `bohr`, and `bohrreduced`.

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
@pure boltzmann(::UnitSystem{k}) where k = measure(k)
@pure planckreduced(::UnitSystem{k,ħ}) where {k,ħ} = measure(ħ)
@pure lightspeed(::UnitSystem{k,ħ,𝘤}) where {k,ħ,𝘤} = measure(𝘤)
@pure permeability(::UnitSystem{k,ħ,𝘤,μ}) where {k,ħ,𝘤,μ} = measure(μ)
@pure electronmass(::UnitSystem{k,ħ,𝘤,μ,m}) where {k,ħ,𝘤,μ,m} = measure(m)
@pure rationalization(::UnitSystem{k,ħ,𝘤,μ,m,λ}) where {k,ħ,𝘤,μ,m,λ} = measure(λ)
@pure lorentz(::UnitSystem{k,ħ,𝘤,μ,m,λ,α}) where {k,ħ,𝘤,μ,m,λ,α} = measure(α)
# ΔνCs:s⁻¹, c:m⋅s⁻¹, h:kg⋅m²⋅s⁻¹, kB:kg⋅m²⋅s⁻²⋅K⁻¹, NA:mol⁻¹, Kcd: cd⋅sr⋅s³⋅kg⁻¹⋅m⁻²

isrationalized(U::UnitSystem) = rationalization(U) ≠ 4π

Base.display(U::UnitSystem) = println("UnitSystem{kB=$(boltzmann(U)),ħ=$(planckreduced(U)),𝘤=$(lightspeed(U)),μ₀=$(permeability(U)),mᵤ=$(electronmass(U)),λ=$(isrationalized(U) ? rationalization(U) : "4π"),αL=$(lorentz(U))}")

@pure universe(::UnitSystem) = Universe
@pure unit(x,y=1) = isapprox(y,x,rtol=eps()^0.9) ? y : x
@pure mass(U::UnitSystem,S::UnitSystem) = electronmass(U,S)
@pure electronmass(𝘩::Float64) = αinv^2*R∞*2𝘩/𝘤
@pure electronmass(𝘩::Float64,C::Coupling) = inv(finestructure(C))^2*R∞*2𝘩/𝘤
@pure planckmass(U::UnitSystem,C::Coupling=universe(U)) = electronmass(U,C)/√coupling(C)
@pure planck(U::UnitSystem,C::Coupling=universe(U)) = 2π*planckreduced(U,C)
@pure newton(U::UnitSystem,C::Coupling=universe(U)) = lightspeed(U,C)*planckreduced(U,C)/planckmass(U,C)^2
@pure charge(U::UnitSystem,C::Coupling=universe(U)) = sqrt(2planck(U,C)*finestructure(C)/impedance(U,C))
@pure impedance(U::UnitSystem,C::Coupling=universe(U)) = permeability(U,C)*lightspeed(U,C)*rationalization(U)*lorentz(U)^2

for unit ∈ (:coupling,:finestructure,:electronunit,:protonunit,:protonelectron)
    @eval @pure $unit(U::UnitSystem) = $unit(universe(U))
end
for unit ∈ (:boltzmann,:planckreduced,:lightspeed,:permeability,:electronmass,:molarmass)
    @eval @pure $unit(U::UnitSystem,C::Coupling) = $unit(U)
end
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

const g₀,ft,ftUS = 9.80665,0.3048,1200/3937; const lbm,lbmUS = g₀/ft,g₀/ftUS
const slug,slugUS,rankine,kelvin,atm = 0.45359237lbm,0.45359237lbmUS,5/9,9/5,101325.0
const kcalₜₕ,kcal₄,kcal₁₀,kcal₂₀,kcalₘ,kcalᵢₜ = 4184,4204,4185.5,4182,4190,4186.8
const calₜₕ,cal₄,cal₁₀,cal₂₀,calₘ,calᵢₜ = (kcalₜₕ,kcal₄,kcal₁₀,kcal₂₀,kcalₘ,kcalᵢₜ)./1e3
const kcal = kcalₜₕ; const cal = kcal/1000 # calₜₕ thermal calorie

# fundamental constants, αinv = (34259-1/4366.8123)/250 # 137.036 exactly?

const ΔνCs,Kcd,mP = 9192631770.0,683*555.016/555,2.176434e-8 # planck mass (kg)
const NA,kB,𝘩,𝘤,𝘦 = 6.02214076e23,1.380649e-23,6.62607015e-34,299792458.,1.602176634e-19
const μₑᵤ,μₚᵤ,αinv,R∞ = 1/1822.888486209,1.007276466621,137.035999084,10973731.5681601
const mₑ,μ₀ = electronmass(𝘩),2𝘩/𝘤/αinv/𝘦^2 # ≈ 4π*(1e-7+5.5e-17), exact charge
const ħ,δμ₀,μₚₑ,Rᵤ,αL,αG = 𝘩/2π,μ₀-4π*1e-7,μₚᵤ/μₑᵤ,NA*kB,0.01/𝘤,(mₑ/mP)^2
const RK1990,RK2014,KJ1990,KJ2014 = 25812.807,25812.8074555,4.835979e14,4.835978525e14
const ħ1990,ħ2014 = 2/RK1990/KJ1990^2/π,2/RK2014/KJ2014^2/π
const mₑ1990,mₑ2014 = electronmass(2π*ħ1990),electronmass(2π*ħ2014)

# engineering units # Thomson: αL = 1/2

const Universe = Coupling{αG,1/αinv,μₑᵤ,μₚᵤ}()
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

const GMsun,GMearth,GMjupiter = 1.32712442099e20,398600441.8e6,1.26686534e17
const au,LD,day = 149597870.7e3,384402e3,60^2*24
const pc,ly,GG = au*648000/π,365.25𝘤*day,𝘤*ħ/mP^2
const mₛ = GMsun/GG; const Jₛ = mₛ*au^2/day^2; export mₛ,Jₛ,au,day
const IAU = UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001/Jₛ,ħ/day/Jₛ,day*𝘤/au,4π*1e-7*day^2/Jₛ,mₑ/mₛ}()

# aliased & humorous units

const mf = mass(90/lbm,Metric,English); const Jf = mf*(201.168/14day)^2
const FFF = UnitSystem{1000Rᵤ*mₑ/μₑᵤ*rankine/Jf,ħ/14day/Jf,14day*𝘤/201.168,0,mₑ/mf}()
const units, US, SI, MKS, temp = UnitSystem, UnitSystem, SI2019, Metric, temperature
const CGS, CGS2019, CGSm, CGSe, HLU = Gauss, EMU2019, EMU, ESU, LorentzHeaviside

# natural units

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

@pure electronmass(U::typeof(Planck),C::Coupling) = sqrt(4π*coupling(C))
@pure electronmass(U::typeof(PlanckGauss),C::Coupling) = sqrt(coupling(C))
@pure electronmass(U::UnitSystem{kB,ħ,𝘤,μ₀,√(αG*αinv)},C::Coupling) where {kB,ħ,𝘤,μ₀} = sqrt(coupling(C)/finestructure(C))
@pure electronmass(U::UnitSystem{kB,ħ,𝘤,μ₀,1/μₚₑ},C::Coupling) where {kB,ħ,𝘤,μ₀} = 1/protonelectron(C)
@pure permeability(U::UnitSystem{kB,ħ,𝘤,4π/αinv^2},C::Coupling) where {kB,ħ,𝘤} = 4π*finestructure(C)^2
@pure permeability(U::UnitSystem{kB,ħ,𝘤,π/αinv^2},C::Coupling) where {kB,ħ,𝘤} = π*finestructure(C)^2
@pure lightspeed(U::UnitSystem{kB,ħ,αinv},C::Coupling) where {kB,ħ} = 1/finestructure(C)
@pure lightspeed(U::UnitSystem{kB,ħ,2αinv},C::Coupling) where {kB,ħ} = 2/finestructure(C)
@pure planckreduced(U::UnitSystem{kB,αinv},C::Coupling) where kB = 1/finestructure(C)

@pure electronmass(U::UnitSystem{kB,ħ,𝘤,μ₀,mₑ},C::Coupling) where {kB,ħ,μ₀} = electronmass(planck(U),C)
@pure electronmass(U::UnitSystem{kB,ħ,100𝘤,μ₀,1000mₑ},C::Coupling) where {kB,ħ,μ₀} = 1000electronmass(SI,C)
@pure electronmass(U::UnitSystem{kB,ħ,𝘤,μ₀,mₑ/1000},C::Coupling) where {kB,ħ,μ₀} = electronmass(SI,C)/1000
@pure electronmass(U::UnitSystem{kB,ħ,𝘤,μ₀,mₑ2014},C::Coupling) where {kB,ħ,μ₀} = electronmass(planck(U),C)
@pure electronmass(U::UnitSystem{kB,ħ,𝘤,μ₀,mₑ1990},C::Coupling) where {kB,ħ,μ₀} = electronmass(planck(U),C)
@pure electronmass(U::UnitSystem{kB,ħ,𝘤/ftUS,μ₀,mₑ/slug},C::Coupling) where {kB,ħ,μ₀} = electronmass(SI,C)/slug
@pure permeability(U::UnitSystem{kB,ħ,𝘤,μ₀},C::Coupling) where {kB,ħ,𝘤} = finestructure(C)*2𝘩/𝘤/𝘦^2
@pure permeability(U::typeof(ESU2019),C::Coupling) = 1e3*permeability(SI,C)/𝘤^2
@pure permeability(U::typeof(EMU2019),C::Coupling) = 1e7*permeability(SI,C)
@pure permeability(U::typeof(CODATA),C::Coupling) = 2RK2014*finestructure(C)/𝘤
@pure permeability(U::typeof(Conventional),C::Coupling) = 2RK1990*finestructure(C)/𝘤

@pure molarmass(U::UnitSystem{1}) = 1
@pure molarmass(U::UnitSystem{kB},C::Coupling=universe(U)) = NA*electronmass(U,C)/electronunit(C)
@pure molarmass(U::UnitSystem{1e7*kB},C::Coupling=universe(U)) = 1000molarmass(SI2019,C)
@pure molarmass(U::UnitSystem{1e3*kB},C::Coupling=universe(U)) = molarmass(SI2019,C)/1000
@pure molarmass(U::UnitSystem{kB}) where kB = molarmass(CGS)/1000
@pure molarmass(U::UnitSystem{boltzmann(MTS)}) = molarmass(CGS)/1e6
@pure molarmass(U::UnitSystem{boltzmann(CGS)}) = molarmass(Natural)
@pure molarmass(U::UnitSystem{boltzmann(FFF)}) = molarmass(Natural)
@pure molarmass(U::UnitSystem{boltzmann(English)},C::Coupling=universe(U)) = 1000molarmass(SI2019,C)
@pure molarmass(U::UnitSystem{boltzmann(EnglishUS)}) = molarmass(Natural)
@pure molarmass(U::UnitSystem{boltzmann(IAU)}) = 1/1000mₛ

@pure luminousefficacy(U::UnitSystem{1}) = 1
@pure luminousefficacy(U::UnitSystem) = power(Kcd,SI2019,U)

include("kinematic.jl")
include("electromagnetic.jl")
include("thermodynamic.jl")
include("physics.jl")
include("systems.jl")

end # module
