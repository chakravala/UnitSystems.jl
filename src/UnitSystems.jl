module UnitSystems

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
#  _  _ _  _ _ ___ ____ _   _ ____ ___ ____ _  _ ____
#  |  | |\ | |  |  [__   \_/  [__   |  |___ |\/| [__
#  |__| | \| |  |  ___]   |   ___]  |  |___ |  | ___]

import Base: @pure, length, time, angle, rem

const Systems = (:Metric,:SI2019,:SI1976,:CODATA,:Conventional,:International,:InternationalMean,:MetricTurn,:MetricDegree,:MetricArcminute,:MetricArcsecond,:MetricGradian,:MetricEngineering,:GravitationalMetric,:MTS,:EMU,:ESU,:Gauss,:LorentzHeaviside,:Kennelly,:FPS,:IPS,:British,:English,:Survey,:FFF,:MPH,:KKH,:Nautical,:Meridian,:IAU☉,:IAUE,:IAUJ,:Hubble,:Cosmological,:CosmologicalQuantum,:Planck,:PlanckGauss,:Stoney,:Hartree,:Rydberg,:Schrodinger,:Electronic,:Natural,:NaturalGauss,:QCD,:QCDGauss,:QCDoriginal)
const Dimensionless = (:coupling,:finestructure,:electronunit,:protonunit,:protonelectron,:darkenergydensity)
const Constants = (:lightspeed,:planck,:planckreduced,:electronmass,:molarmass,:boltzmann,:permeability,:rationalization,:lorentz,:luminousefficacy,:gravity) #angle
const Physics = (:turn,:spat,:dalton,:protonmass,:planckmass,:gravitation,:gaussgravitation,:einstein,:hartree,:rydberg,:bohr,:electronradius,:avogadro,:molargas,:stefan,:radiationdensity,:vacuumpermeability,:vacuumpermittivity,:electrostatic,:magnetostatic,:biotsavart,:elementarycharge,:faraday,:vacuumimpedance,:conductancequantum,:klitzing,:josephson,:magneticfluxquantum,:magneton)
const Derived = (:hyperfine,:loschmidt,:wienwavelength,:wienfrequency,:mechanicalheat,:eddington,:solarmass,:jupitermass,:earthmass,:lunarmass,:earthradius,:greatcircle,:radarmile,:hubble,:cosmological,
    :radian,:steradian,:degree,:squaredegree,:gradian,:bradian,:arcminute,:arcsecond,
    :second,:minute,:hour,:day,:gaussianmonth,:siderealmonth,:synodicmonth,:year,:gaussianyear,:siderealyear,:jovianyear,
    :angstrom,:inch,:foot,:surveyfoot,:yard,:meter,:earthmeter,:mile,:statutemile,:meridianmile,:admiraltymile,:nauticalmile,:lunardistance,:astronomicalunit,:jupiterdistance,:lightyear,:parsec,
    :barn,:hectare,:acre,:surveyacre,
    :liter,:gallon,:quart,:pint,:cup,:fluidounce,:teaspoon,:tablespoon,
    :bubnoff,:ips,:fps,:fpm,:ms,:kmh,:mph,:knot,:mps,
    :grain,:gram,:earthgram,:kilogram,:tonne,:ton,:pound,:ounce,:slug,:slinch,:hyl,
    :dyne,:newton,:poundal,:poundforce,:kilopond,
    :psi,:pascal,:bar,:barye,:technicalatmosphere,:atmosphere,:inchmercury,:torr,
    :electronvolt,:erg,:joule,:footpound,:calorie,:kilocalorie,:meancalorie,:earthcalorie,:thermalunit,:gasgallon,:tontnt,
    :watt,:horsepower,:horsepowerwatt,:horsepowermetric,:electricalhorsepower,:tonsrefrigeration,:boilerhorsepower,
    :coulomb,:earthcoulomb,:ampere,:volt,:henry,:ohm,:siemens,:farad,:weber,:tesla,
    :abcoulomb,:abampere,:abvolt,:abhenry,:abohm,:abmho,:abfarad,:maxwell,:gauss,:oersted,:gilbert,
    :statcoulomb,:statampere,:statvolt,:stathenry,:statohm,:statmho,:statfarad,:statweber,:stattesla,
    :kelvin,:rankine,:celsius,:fahrenheit,:sealevel,:boiling,:mole,:earthmole,:poundmole,:slugmole,:slinchmole,:katal,:amagat,
    :lumen,:candela,:lux,:phot,:footcandle,:nit,:apostilb,:stilb,:lambert,:footlambert,:bril,
    :neper,:bel,:decibel,:hertz,:apm,:rpm,
    :kayser,:diopter,:gforce,:galileo,:eotvos,:darcy,:poise,:reyn,:stokes,:rayl,
    :mpge,:langley,:jansky,:solarflux,:curie,:sievert,:roentgen,:rem)
const Kinematic = (:solidangle,:time,:length,:area,:volume,:wavenumber,:angularwavenumber,:fuelefficiency,:numberdensity,:frequency,:angularfrequency,:frequencydrift,:speed,:acceleration,:jerk,:snap,:crackle,:pop,:volumeflow) #angle
const Mechanical = (:inertia,:mass,:massflow,:lineardensity,:areadensity,:density,:specificweight,:specificvolume,:force,:specificforce,:gravityforce,:pressure,:compressibility,:viscosity,:diffusivity,:rotationalinertia,:impulse,:momentum,:angularmomentum,:yank,:energy,:specificenergy,:action,:fluence,:power,:powerdensity,:irradiance,:radiance,:radiantintensity,:spectralflux,:spectralexposure,:soundexposure,:impedance,:specificimpedance,:admittance,:compliance,:inertance)
const Electromagnetic = (:charge,:chargedensity,:linearchargedensity,:exposure,:mobility,:current,:currentdensity,:resistance,:conductance,:resistivity,:conductivity,:capacitance,:inductance,:reluctance,:permeance,:permittivity,:permeability,:susceptibility,:specificsusceptibility,:demagnetizingfactor,:vectorpotential,:electricpotential,:magneticpotential,:electricfield,:magneticfield,:electricflux,:magneticflux,:electricfluxdensity,:magneticfluxdensity,:electricdipolemoment,:magneticdipolemoment,:electricpolarizability,:magneticpolarizability,:magneticmoment,:specificmagnetization,:polestrength)
const Thermodynamic = (:temperature,:entropy,:specificentropy,:volumeheatcapacity,:thermalconductivity,:thermalconductance,:thermalresistivity,:thermalresistance,:thermalexpansion,:lapserate)
const Molar = (:molarmass,:molality,:molaramount,:molarity,:molarvolume,:molarentropy,:molarenergy,:molarconductivity,:molarsusceptibility,:catalysis,:specificity)
const Photometric = (:luminousflux,:luminousintensity,:luminance,:illuminance,:luminousenergy,:luminousexposure,:luminousefficacy)
const Mechanics = [Kinematic...,Mechanical...]
const Convert = [Mechanics...,Electromagnetic...,Thermodynamic...,Molar...,Photometric...]

listext(x) = join(x,"`, `")

"""
    UnitSystems.similitude() = haskey(ENV,"SIMILITUDE")

An optional environment variable `ENV["SIMILITUDE"]` induces `UnitSystems.similitude()` to return `true`, giving flexibility for building dependencies whenever it is desirable to toggle usage between `UnitSystems` (default) and `Similitude` (requires environment variable specification). For example, in `MeasureSystems` and `Geophysics` this option is used to increase flexibility with variety in local compilation workflow.
"""
similitude() = haskey(ENV,"SIMILITUDE")

include("constant.jl")
Constant(D::Number) = Constant{D}()
@pure Constant(D::Float64) = Constant{D}()
@pure Constant(D::Int) = Constant{D}()
@pure Constant(D::Constant) = D
@pure constant(::Constant{D}) where D = D
Base.Int(::Constant{D}) where D = Constant(Int(D))
Base.show(io::IO,x::Constant{D}) where D  = show(io,D)
unit(x,::Constant{y}) where y = unit(x,y)
unit(::Constant{x},::Constant{y}) where {x,y} = Constant{unit(x,y)}()
unit(::Constant{x},y=1) where x = Constant{unit(x,y)}()

@pure measure(x) = x
@pure cache(x) = x
Constant(x) = x
Quantity(x) = x
Quantity(D,U,x) = x
const 𝟙,F,M,L,T,Q,Θ,N,J,A,Λ,C = 1,1,1,1,1,1,1,1,1,1,1,1

logdb(x) = 10log10(x)
expdb(x) = exp10(0.1)^x
const dB, Db = logdb, expdb

# universe

"""
    Coupling{αG,α,μₑᵤ,μₚᵤ,ΩΛ}

Specification of `Universe` with the dimensionless `Coupling` constants `coupling`, `finestructure`, `electronunit`, `protonunit`, `protonelectron`, and `darkenergydensity`.
Alterations to these values can be facilitated and quantified using parametric polymorphism.
Due to the `Coupling` interoperability, the `MeasureSystems` package is made possible to support calculations with `Measurements` having error standard deviations.
"""
struct Coupling{αG,α,μₑᵤ,μₚᵤ,ΩΛ} end
Coupling{αG,α,μₑᵤ,μₚᵤ}() where {αG,α,μₑᵤ,μₚᵤ} = Coupling(αG,α,μₑᵤ,μₚᵤ,ΩΛ)
Coupling(αG,α,μₑᵤ,μₚᵤ,ΩΛ) = Coupling{cache(αG),cache(α),cache(μₑᵤ),cache(μₚᵤ),cache(ΩΛ)}()
@pure coupling(U::Coupling{αG}) where αG = measure(αG)
@pure finestructure(U::Coupling{αG,α}) where {αG,α} = measure(α)
@pure electronunit(U::Coupling{αG,α,μₑᵤ}) where {αG,α,μₑᵤ} = measure(μₑᵤ)
@pure protonunit(U::Coupling{αG,α,μₑᵤ,μₚᵤ}) where {αG,α,μₑᵤ,μₚᵤ} = measure(μₚᵤ)
@pure protonelectron(U::Coupling) = protonunit(U)/electronunit(U)
@pure darkenergydensity(U::Coupling{αG,α,μₑᵤ,μₚᵤ,ΩΛ}) where {αG,α,μₑᵤ,μₚᵤ,ΩΛ} = measure(ΩΛ)
Base.display(U::Coupling) = println("Coupling{αG = $(coupling(U)), α = $(finestructure(U)), μₑᵤ = $(electronunit(U)), μₚᵤ = $(protonunit(U)), ΩΛ = $(darkenergydensity(U))}")

# unit systems

"""
    UnitSystem{kB, ħ, 𝘤, μ₀, mₑ, Mᵤ, (Kcd, θ, λ, αL, g₀, ...)}

Fundamental constants of physics are: `kB` Boltzmann's constant, `ħ` reduced Planck's constant, `𝘤` speed of light, `μ₀` vacuum permeability, `mₑ` electron rest mass, `Mᵤ` molar mass, `Kcd` luminous efficacy, `θ` angle measure, `λ` Gauss rationalization, `αL` Lorentz's constant, and `g₀` gravitational force reference.
Primarily the `Metric` SI unit system is used in addition to the historic `English` engineering unit system.
These constants induce derived values for `avogadro`, `boltzmann`, `molargas`, `planck`, `planckreduced`, `lightspeed`, `planckmass`, `dalton`, `protonmass`, `electronmass`, `newton`, `einstein`, `vacuumpermeability`, `vacuumpermittivity`, `electrostatic`, and
additional constants `molarmass`, `luminousefficacy`, `gravity`, `angle`, `turn`, `spat`, `stefan`, `radiationdensity`, `magnetostatic`, `lorentz`, `biotsavart`, `rationalization`, `vacuumimpedance`, `elementarycharge`, `magneton`, `conductancequantum`, `faraday`, `magneticfluxquantum`, `josephson`, `klitzing`, `hartree`, `rydberg`, `bohr`.

Standardized unit/derived quantities are `$(listext(Derived))`.

Additional reference `UnitSystem` variants: `EMU`, `ESU`, `Gauss`, `LorentzHeaviside`, `SI2019`, `SI1976`, `CODATA`, `Conventional`, `International`, `InternationalMean`, `MetricEngineering`, `GravitationalMetric`, `IAU`, `IAUE`, `IAUJ`, `FPS`, `IPS`, `British`, `Survey`, `Hubble`, `Cosmological`, `CosmologicalQuantum`, `Meridian`, `Nautical`, `MPH`, `KKH`, `MTS`, `FFF`; and natural atomic units based on gravitational `coupling` and `finestructure` constant (`Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal`).

**Dimensional unit conversions:**

Mechanics: `angle`, `$(listext(Kinematic))`, `$(listext(Mechanical))`;
Electromagnetics: `$(listext(Electromagnetic))`;
Thermodynamics: `$(listext(Thermodynamic))`,
`$(listext(Molar))`, `$(listext(Photometric))`.
""" #`Rᵤ,Da,σ,ħ,μ₀,ε₀,kₑ,𝘦,𝔉,RK,Z₀,G₀`
struct UnitSystem{kB,ħ,𝘤,μ₀,mₑ,Mᵤ,extra}
    @pure UnitSystem{kB,ħ,𝘤,μ₀,mₑ,Mᵤ,extra}() where {kB,ħ,𝘤,μ₀,mₑ,Mᵤ,extra} = new{kB,ħ,𝘤,μ₀,mₑ,Mᵤ,extra}()
end # UnitSystem{kB,ħ,𝘤,μ₀,mₑ,Mᵤ,(Kcd,θ,λ,αL,g,C,τ,𝟐,𝟑,𝟓,𝟕,𝟏𝟏,𝟏𝟗,𝟒𝟑)}
function UnitSystem(kB,ħ,𝘤,μ₀,mₑ,Mᵤ=1,Kcd=1,θ=1,λ=1,αL=1,g=1,C=Universe,τ=2π,x=2,y=3,z=5,w=7,u=11,v=19,q=43)
    UnitSystem{cache(kB),cache(ħ),cache(𝘤),cache(μ₀),cache(mₑ),cache(Mᵤ),(cache(Kcd),cache(θ),cache(λ),cache(αL),cache(g),C,cache(τ),cache(x),cache(y),cache(z),cache(w),cache(u),cache(v),cache(q))}()
end
@pure boltzmann(::UnitSystem{k}) where k = measure(k)
@pure planckreduced(::UnitSystem{k,ħ}) where {k,ħ} = measure(ħ)
@pure lightspeed(::UnitSystem{k,ħ,𝘤}) where {k,ħ,𝘤} = measure(𝘤)
@pure permeability(::UnitSystem{k,ħ,𝘤,μ}) where {k,ħ,𝘤,μ} = measure(μ)
@pure vacuumpermeability(U::UnitSystem) = permeability(U)
@pure electronmass(::UnitSystem{k,ħ,𝘤,μ,m}) where {k,ħ,𝘤,μ,m} = measure(m)
@pure molarmass(::UnitSystem{k,ħ,𝘤,μ,m,M}) where {k,ħ,𝘤,μ,m,M} = measure(M)
@pure luminousefficacy(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[1])
@pure angle(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[2])
@pure rationalization(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[3])
@pure lorentz(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[4])
@pure gravity(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[5])
@pure universe(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[6])
@pure twopi(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[7])
@pure two(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[8])
@pure three(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[9])
@pure five(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[10])
@pure seven(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[11])
@pure eleven(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[12])
@pure nineteen(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[13])
@pure fourtythree(::UnitSystem{k,ħ,𝘤,μ,m,M,e}) where {k,ħ,𝘤,μ,m,M,e} = measure(e[14])
# ΔνCs:s⁻¹, c:m⋅s⁻¹, h:kg⋅m²⋅s⁻¹, kB:kg⋅m²⋅s⁻²⋅K⁻¹, NA:mol⁻¹, Kcd: cd⋅sr⋅s³⋅kg⁻¹⋅m⁻²

function evaldim end
@pure isquantity(U) = false
@pure isquantity(A,B) = isquantity(A) && isquantity(B)
@pure isquantity(U::UnitSystem) = isquantity(boltzmann(U))
@pure isrationalized(U::UnitSystem) = rationalization(U) ≠ sphere(U)

normal(x) = x
@pure normal(x::Float64) = x
@pure normal(x::Int) = x
@pure normal(::UnitSystem{kB,ħ,𝘤,μ0,me,Mu,e}) where {kB,ħ,𝘤,μ0,me,Mu,e} = UnitSystem(normal(kB),normal(ħ),normal(𝘤),normal(μ0),normal(me),normal(Mu),normal(e[1]),normal(e[2]),normal(e[3]),normal(e[4]),normal(e[5]),e[6],normal(e[7]),normal(e[8]),normal(e[9]),normal(e[10]),normal(e[11]),normal(e[12]),normal(e[13]),normal(e[14]))

@pure unitname(::UnitSystem) = "Unknown"
Base.show(io::IO,U::UnitSystem) = print(io,unitname(normal(U)))
function Base.display(U::UnitSystem)
    println("UnitSystem: ", unitname(normal(U)))
    println("  kB  = $(boltzmann(U))")
    println("  ħ   = $(planckreduced(U))")
    println("  𝘤   = $(lightspeed(U))")
    println("  μ₀  = $(vacuumpermeability(U))")
    println("  mₑ  = $(electronmass(U))")
    println("  Mᵤ  = $(molarmass(U))")
    println("  Kcd = $(luminousefficacy(U))")
    println("  θ   = $(angle(U))")
    println("  λ   = $(rationalization(U)≠4π ? rationalization(U) : "4π")")
    println("  αL  = $(lorentz(U))")
    println("  g₀  = $(gravity(U))")
end

function (U::UnitSystem)(JK=1,Js=1,ms=1,Hm=1,kg=1)
    kB = boltzmann(U)*JK
    ħ = planckreduced(U)*Js
    𝘤 = lightspeed(U)*ms
    μ₀ = vacuumpermeability(U)*Hm
    mₑ = electronmass(U)*kg
    Mᵤ = molarmass(U)#kg#
    Kcd = luminousefficacy(U)#
    A = angle(U)
    λ = rationalization(U)
    αL = lorentz(U)
    g₀ = gravity(U)
    C = universe(U)
    τ = twopi(U)
    (x,y,z,u,v,w,p) = (two(U),three(U),five(U),seven(U),eleven(U),nineteen(U),fourtythree(U))
    UnitSystem(kB,ħ,𝘤,μ₀,mₑ,Mᵤ,Kcd,A,λ,isone(αL) ? αL : αL/ms,g₀,C,τ,x,y,z,u,v,w,p)
end

function DimensionSystem(U::UnitSystem,L,M,T,I,Θ,N,J,A,Λ,G,C)
    kB = boltzmann(U)*(M*L*L/(T*T)/Θ/G)
    ħ = reducedplanck(U)*(M*(L*L)/T/G/A)
    c = lightpseed(U)*(L/T)
    μ0 = vacuumpermeability(U)*(L*M/((T*T)*(I*I)*(C*C)*(A*A)*Λ))
    me = electronmass(U)*M
    Mu = molarmass(U)*(M/N)
    Kcd = luminousefficacy(U)*(J*(T*T*T)/(L*L*M)*G)
    λ = rationalization(U)*(Λ*(A*A))
    αL = lorentz(U)*C
    g₀ = gravity(U)*G
    τ = twopi(U)
    (x,y,z,u,v,w,p) = (two(U),three(U),five(U),seven(U),eleven(U),fourtythree(U))
    UnitSystem(kB,ħ,c,μ0,me,Mu,Kcd,A,λ,αL,g₀,universe(U),τ,x,y,z,u,v,w,p)
end

"""
    ElectricSystem(U::UnitSystem,Ω,V) = EntropySystem(U,𝟏,𝟏,V^2/Ω,𝟏,vacuumpermeability(U)/Ω)

Constructs new `UnitSystem` from `U` with `mass` rescaled by `electricpotential` and `resistance`. In the `International` system, `Ωᵢₜ` and `Vᵢₜ` are used as definitions from the more recent United States results, while in `InternationalMean` an earlier estimate based on other nations was used.
"""
ElectricSystem(u,Ω,V) = EntropySystem(u,one(u),one(u),V*V/Ω,one(u),vacuumpermeability(u)/Ω)

"""
    GaussSystem(U::UnitSystem,μ0,λ,αL=𝟏,l=centi,m=milli,g0=gravity(U))

Constructs new `UnitSystem` from `U` rescaled for `CGS` with electromagnetic options. The first three options are to set the values for `vacuumpermeability`, `rationalization`, and `lorentz` constants. The following two parameters are scaling for `length` and `mass`, while the last is an option to change the `gravity` reference.

Examples include `EMU`, `ESU`, `Gauss`, `LorentzHeaviside`, and `Kennelly`.
"""
function GaussSystem(u,μ0,λ,αL=one(u),l=inv((two(u)*five(u))^2),m=inv((two(u)*five(u))^3),g0=gravity(u))
    EntropySystem(u,one(u),l,m,one(u),μ0,m==1/1000 ? one(u) : molarmass(u)/m,g0,m*(l*l),λ,αL)
end

"""
    EntropySystem(U::UnitSystem,t,l,m,θ=𝟏)
    EntropySystem(U::UnitSystem,t,l,m,θ,μ0,Mu=molarmass(U)/m,g0=gravity(U))

Constructs new `UnitSystem` from `U` rescaled along `time`, `length`, `mass`, and `temperature` by the first four parameters. Additional optional parameters allow for customization of the `vacuumpermeability`, `molarmass`, and `gravity` constants.

Examples of this type include `Nautical`, `Meridian`, `MeridianEngineering`, `GravitatonalMeridian`, `GravitationalMetric`, `GravitationalSI2019`, `MTS`, `KKH`, `MPH`, `IAU☉`, `IAUE`, `IAUJ`, `Hubble`, `Cosmological`, `CosmologicalQuantum`.
However, most other constructors for `UnitSystem` derivations are based on internally calling `EntropySystem`, such as `AstronomicalSystem`, `ElectricSystem`, `GaussSystem`, and `RankineSystem`.
This means `EntropySystem` also constructs the examples listed there.
"""
function EntropySystem(u,t,l,m,θ=one(u))
    EntropySystem(u,t,l,m,θ,permeability(u)/(m*l),molarmass(u)/m,gravity(u),m*l*l/(t*t))
end
function EntropySystem(u,t,l,m,θ,μ0,Mu=molarmass(u)/m,g0=gravity(u),e=m*l*l/(t*t),λ=one(u),αL=one(u),Kcd=luminousefficacy(u)*e/t*g0)
    normal(UnitSystem(
        boltzmann(u)*θ/e/g0,
        planckreduced(u)/t/e/g0,
        lightspeed(u)*t/l,
        μ0,
        electronmass(u)/m,
        Mu,
        Kcd,
        angle(u),λ,αL,g0,universe(u),twopi(u),two(u),three(u),five(u),seven(u),eleven(u),nineteen(u),fourtythree(u)))
end

"""
    AstronomicalSystem(U::UnitSystem,t,l,m)

Constructs new `UnitSystem` from `U` rescaled along `time`, `length`, `mass`, and dimensionless `boltzmann` and `molarmass` constants.
Examples are `Hubble`, `Cosmological`, `CosmologicalQuantum`.
"""
function AstronomicalSystem(u,t,l,m,e=m*lightspeed(u)^2)
    EntropySystem(u,t,l,m,e/boltzmann(u),spat(u),one(u),one(u),e,one(u),one(u),one(u))
end

@pure unit(x,y=1) = isapprox(y,x,rtol=eps()^0.9) ? y : x
@pure Base.one(U::UnitSystem) = unit(two(U)/two(U))
@pure Base.zero(U::UnitSystem) = one(U)-one(U)
@pure turn(U::UnitSystem) = twopi(U)/angle(U)
@pure solidangle(U::UnitSystem,S::UnitSystem) = unit(angle(U,S)^2)
@pure spat(U::UnitSystem) = two(U)*turn(U)/angle(U)*unit(turn(U)/normal(turn(U)))
@pure mass(U::UnitSystem,S::UnitSystem) = electronmass(U,S)
@pure electronmass(𝘩::Number,C::Coupling) = inv(finestructure(C))^2*R∞*2𝘩/𝘤
@pure planckmass(U::UnitSystem,C::Coupling=universe(U)) = electronmass(U,C)/√coupling(C)
@pure planck(U::UnitSystem,C::Coupling=universe(U)) = turn(U)*planckreduced(U,C)
@pure gravitation(U::UnitSystem,C::Coupling=universe(U)) = lightspeed(U,C)*planck(U,C)/normal(twopi(U))/planckmass(U,C)^2
@pure elementarycharge(U::UnitSystem,C::Coupling=universe(U)) = sqrt(two(U)*planck(U)/(vacuumpermeability(U)/finestructure(C))/(lightspeed(U)*rationalization(U)*lorentz(U)^2))

for unit ∈ (:coupling,:finestructure,:electronunit,:protonunit,:protonelectron,:darkenergydensity)
    @eval @pure $unit(U::UnitSystem) = $unit(universe(U))
end
for unit ∈ (:boltzmann,:planckreduced,:lightspeed,:vacuumpermeability,:permeability,:electronmass,:molarmass)
    @eval @pure $unit(U::UnitSystem,C::Coupling) = $unit(U)
end
for unit ∈ (Constants...,:angle,:vacuumpermeability)
    @eval @pure $unit(U::UnitSystem,S::UnitSystem) = isquantity(U,S) ? evaldim($unit)(U,S) : unit($unit(S)/$unit(U))
end
for unit ∈ (Convert...,:angle)
    @eval begin
        @pure @inline $unit(v::Real,U::UnitSystem) = isquantity(U) ? evaldim($unit)(v,U) : $unit(v,U,Metric)
        @pure @inline $unit(v::Real,U::UnitSystem,S::UnitSystem) = isquantity(U,S) ? evaldim($unit)(v,U,S) : (u=$unit(U,S);isone(u) ? v : v/u)
        @pure @inline $unit(v::Real,U::UnitSystem{kB,ħ,𝘤,μ₀,mₑ,Mᵤ,extra},S::UnitSystem{kB,ħ,𝘤,μ₀,mₑ,Mᵤ,extra}) where {kB,ħ,𝘤,μ₀,mₑ,Mᵤ,extra} = v
    end
    if unit ∉ (Constants...,:angle,:permeability)
        @eval @pure @inline $unit(U::UnitSystem) = isquantity(U) ? evaldim($unit)(U) : $unit(Natural,U)
    end
end
for unit ∈ (Systems...,Dimensionless...,Constants...,Physics...,Derived...,Convert...)
    @eval export $unit
end

# fundamental constants, αinv = (34259-1/4366.8123)/250 # 137.036 exactly?

const g₀,atm,T₀ = Constant(9.80665),Constant(101325.0),Constant(273.15)
const ft,ftUS,lb = Constant(0.3048),Constant(1200/3937),Constant(0.45359237)
const inHg,Ωᵢₜ,Vᵢₜ = Constant(1/3386.389),Constant(1.000495),Constant(1.00033)
const ΔνCs,Kcd,mP = Constant(9192631770.0),Constant(683*555.016/555),Constant(2.176434e-8)
const NA,kB,𝘩 = Constant(6.02214076e23),Constant(1.380649e-23),Constant(6.62607015e-34)
const 𝘤,𝘦 = Constant(299792458.),Constant(1.602176634e-19)
const μₑᵤ,μₚᵤ,μE☾ = Constant(1/1822.888486209),Constant(1.007276466621),Constant(81.300568)
const αinv,R∞ = Constant(137.035999084),Constant(10973731.5681601)
const RK1990,KJ1990,Rᵤ2014 = Constant(25812.807),Constant(4.835979e14),Constant(8.3144598)
const RK2014,KJ2014 = Constant(25812.8074555),Constant(4.835978525e14)
const GME,GMJ = Constant(398600441.8e6),Constant(1.26686534e17)
const kG,H0,ΩΛ = Constant(3548.18761),Constant(67.66),Constant(0.6889)
const aⱼ,au = Constant(365.25),Constant(149597870.7e3)
const LD,JD = Constant(384399e3),Constant(778479e6)
const zetta,zepto = Constant(1e21),Constant(1e-21)
const yotta,yocto = Constant(1e24),Constant(1e-24)
const 𝟏,𝟐,𝟑,𝟓,𝟕,𝟏𝟎,𝟏𝟏,𝟏𝟗,𝟒𝟑,τ,α = Constant(1),Constant(2),Constant(3),Constant(5),Constant(7),Constant(10),Constant(11),Constant(19),Constant(43),Constant(2π),inv(αinv)

include("initdata.jl")

const slug,lbm,lbmUS,rankine,kelvin = lb*g₀/ft,g₀/ft,g₀/ftUS,°R,K
const ħ1990,ħ2014 = planckreduced(Conventional),planckreduced(CODATA)
const mₑ1990,mₑ2014 = electronmass(Conventional),electronmass(CODATA)
const δμ₀,ly,mₛ,GG = μ₀-4π*1e-7,aⱼ*𝘤*DAY,GM☉/G,G
const units, temp = UnitSystem, temperature

# common conversion factors

const kcalₜₕ,kcal₄,kcal₁₀,kcal₂₀,kcalₘ,kcalᵢₜ = 4184,4204,4185.5,4182,4190,4186.8
const calₜₕ,cal₄,cal₁₀,cal₂₀,calₘ,calᵢₜ = (kcalₜₕ,kcal₄,kcal₁₀,kcal₂₀,kcalₘ,kcalᵢₜ)./1e3

@pure deka(U::UnitSystem) = two(U)*five(U)
@pure hecto(U::UnitSystem) = deka(U)^2
@pure kilo(U::UnitSystem) = deka(U)^3
@pure mega(U::UnitSystem) = (Constant(1.0)*kilo(U))^2
@pure giga(U::UnitSystem) = (Constant(1.0)*kilo(U))^3
@pure tera(U::UnitSystem) = (Constant(1.0)*kilo(U))^4
@pure peta(U::UnitSystem) = (Constant(1.0)*kilo(U))^5
@pure exa(U::UnitSystem) = (Constant(1.0)*kilo(U))^6
@pure zetta(U::UnitSystem) = (Constant(1.0)*kilo(U))^7
@pure yotta(U::UnitSystem) = (Constant(1.0)*kilo(U))^8
@pure deci(U::UnitSystem) = inv(deka(U))
@pure centi(U::UnitSystem) = inv(hecto(U))
@pure milli(U::UnitSystem) = inv(kilo(U))
@pure micro(U::UnitSystem) = inv(mega(U))
@pure nano(U::UnitSystem) = inv(giga(U))
@pure pico(U::UnitSystem) = inv(tera(U))
@pure femto(U::UnitSystem) = inv(peta(U))
@pure atto(U::UnitSystem) = inv(exa(U))
@pure zepto(U::UnitSystem) = inv(zetta(U))
@pure yocto(U::UnitSystem) = inv(yotta(U))
@pure kibi(U::UnitSystem) = two(U)^10
@pure mebi(U::UnitSystem) = two(U)^20
@pure gibi(U::UnitSystem) = two(U)^30
@pure tebi(U::UnitSystem) = two(U)^40
@pure pebi(U::UnitSystem) = two(U)^50
@pure exbi(U::UnitSystem) = two(U)^60
@pure zebi(U::UnitSystem) = (Constant(1.0)*two(U))^70
@pure yobi(U::UnitSystem) = (Constant(1.0)*two(U))^80

# physical constants

@pure electronmass(U::typeof(Planck),C::Coupling) = sqrt(spat(U)*coupling(C))
@pure electronmass(U::typeof(PlanckGauss),C::Coupling) = sqrt(coupling(C))
@pure electronmass(U::UnitSystem{kB,ħ,𝘤,μ₀,√(αG*αinv)},C::Coupling) where {kB,ħ,𝘤,μ₀} = sqrt(coupling(C)/finestructure(C))
@pure electronmass(U::UnitSystem{kB,ħ,𝘤,μ₀,1/μₚₑ},C::Coupling) where {kB,ħ,𝘤,μ₀} = 1/protonelectron(C)
@pure vacuumpermeability(U::UnitSystem{kB,ħ,𝘤,4π/αinv^2},C::Coupling) where {kB,ħ,𝘤} = spat(U)*finestructure(C)^2
@pure vacuumpermeability(U::UnitSystem{kB,ħ,𝘤,π/αinv^2},C::Coupling) where {kB,ħ,𝘤} = spat(U)/two(U)^2*finestructure(C)^2
@pure lightspeed(U::UnitSystem{kB,ħ,αinv},C::Coupling) where {kB,ħ} = inv(finestructure(C))
@pure lightspeed(U::UnitSystem{kB,ħ,2αinv},C::Coupling) where {kB,ħ} = two(U)/finestructure(C)
@pure planckreduced(U::UnitSystem{kB,αinv},C::Coupling) where kB = inv(finestructure(C))

@pure electronmass(U::UnitSystem{kB,ħ,𝘤,μ₀,mₑ},C::Coupling) where {kB,ħ,μ₀} = electronmass(U)#electronmass(planck(U),C)
@pure electronmass(U::UnitSystem{kB,ħ,100𝘤,μ₀,1000mₑ},C::Coupling) where {kB,ħ,μ₀} = electronmass(SI,C)*(two(U)*five(U))^3
@pure electronmass(U::UnitSystem{kB,ħ,𝘤,μ₀,mₑ/1000},C::Coupling) where {kB,ħ,μ₀} = electronmass(SI,C)/(two(U)*five(U))^3
@pure electronmass(U::UnitSystem{kB,ħ,𝘤,μ₀,electronmass(CODATA)},C::Coupling) where {kB,ħ,μ₀} = electronmass(planck(U),C)
@pure electronmass(U::UnitSystem{kB,ħ,𝘤,μ₀,electronmass(Conventional)},C::Coupling) where {kB,ħ,μ₀} = electronmass(planck(U),C)
@pure electronmass(U::UnitSystem{kB,ħ,𝘤/ftUS,μ₀,mₑ*ft/lb/g₀},C::Coupling) where {kB,ħ,μ₀} = electronmass(SI,C)*ft/lb/g₀
@pure vacuumpermeability(U::UnitSystem{kB,ħ,𝘤,μ₀},C::Coupling) where {kB,ħ,𝘤} = finestructure(C)*2𝘩/𝘤/𝘦^2
#@pure vacuumpermeability(U::typeof(ESU2019),C::Coupling) = 1e3*vacuumpermeability(SI,C)/𝘤^2
#@pure vacuumpermeability(U::typeof(EMU2019),C::Coupling) = 1e7*vacuumpermeability(SI,C)
@pure vacuumpermeability(U::typeof(CODATA),C::Coupling) = 2RK2014*finestructure(C)/𝘤
@pure vacuumpermeability(U::typeof(Conventional),C::Coupling) = 2RK1990*finestructure(C)/𝘤

convertext(unit,fun) = """
```Julia
$unit(U::UnitSystem,S::UnitSystem) = $fun
$unit(v::Real,U::UnitSystem,S::UnitSystem) = v/$unit(U,S)
```
"""

unitext(unit,text) = """
```Julia
$unit(U::UnitSystem) = $text
```
"""

include("kinematic.jl")
include("electromagnetic.jl")
include("thermodynamic.jl")
include("physics.jl")
include("systems.jl")

const RK = klitzing(SI2019) #
const KJ = josephson(SI2019) #

end # module
