
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

const deka,byte = 𝟐*𝟓,𝟐^3
const hecto, kilo =  deka^2, deka^3
const mega,giga,tera,peta,exa = kilo^2, kilo^3, kilo^4, kilo^5, kilo^6
const deci,centi,milli,micro,nano,pico,femto,atto = inv(deka),inv(hecto),inv(kilo),inv(mega),inv(giga),inv(tera),inv(peta),inv(exa)
const kibi,mebi,gibi,tebi,pebi,exbi,zebi,yobi = 𝟐^10,𝟐^20,𝟐^30,𝟐^40,𝟐^50,𝟐^60,(Constant(1.0)*𝟐)^70,(Constant(1.0)*𝟐)^80

const fur,°R,K,HOUR,k = 𝟐^2*𝟑*𝟓*𝟏𝟏*ft,𝟓/𝟑^2,𝟑^2/𝟓,𝟐^4*𝟑^2*𝟓^2,kG*τ/(𝟐^7*𝟑^4*𝟓^3)
const mₑ,μ₀ = αinv^2*R∞*𝟐*𝘩/𝘤,𝟐*𝘩/𝘤*α/𝘦^2 # ≈ 4π*(1e-7+5.5e-17), exact charge
const ħ,μₚₑ,μₑₚ,Rᵤ,αL,αG,Mᵤ = 𝘩/τ,μₚᵤ/μₑᵤ,μₑᵤ/μₚᵤ,NA*kB,centi/𝘤,(mₑ/mP)^2,NA*mₑ/μₑᵤ
const pc,G,DAY,nm = au*𝟐^7*𝟑^4*𝟓^3/τ,𝘤*ħ/mP^2,𝟐^7*𝟑^3*𝟓^2,sqrt(GME/g₀)*τ/𝟐^5/𝟑^3/𝟓^2
const GM☉ =au^3*k^2/DAY^2; const th = 𝟏𝟎^3*pc/H0; const ΛC = 𝟑*ΩΛ*(th*𝘤)^-2
const lc,mc,ρΛ,𝘦ₙ = 𝟐*sqrt(τ/ΛC),𝘤^2/(𝟐*G*sqrt(τ*ΛC)),ΛC*𝘤^4/(𝟐^2*τ)/G,𝘦/√α
const lcq,mcq = sqrt.(sqrt.((𝘤*ħ/ρΛ,ρΛ*ħ^3/𝘤^5))); const 𝘦ᵣ = 𝘦ₙ/√(𝟐*τ)
const tcq,em,mi = lcq*sqrt(mcq/sqrt(sqrt(ρΛ*(𝘤*ħ)^3))),sqrt(GME/g₀)*τ/𝟐^9/𝟓^7,𝟐^5*𝟑*𝟓*𝟏𝟏

@pure sackurtetrode(U::UnitSystem,P=atm,T=𝟏,m=dalton(U)) = normal(log((Constant(exp(5/2))*boltzmann(U)*sqrt(boltzmann(U)/gravity(U)/turn(U)/planckreduced(U)^2)^3)*(T/P*sqrt(m*T)^3)))

const Universe = Coupling(αG,α,μₑᵤ,μₚᵤ,ΩΛ)

Constant(::UnitSystem{a,b,c,d,e,f,g}) where {a,b,c,d,e,f,g} = UnitSystem(Constant(a),Constant(b),Constant(c),Constant(d),Constant(e),Constant(f),Constant(g[1]),Constant(g[2]),Constant(g[3]),Constant(g[4]),Constant(g[5]),Universe,Constant(g[7]),𝟐,𝟑,𝟓,𝟕,𝟏𝟏,𝟏𝟗,𝟒𝟑)

export MetricSystem, ConventionalSystem, RankineSystem
export AstronomicalSystem, ElectricSystem, GaussSystem, EntropySystem

"""
    MetricSystem(Mu=Mᵤ,μ0=μ₀,Ru=Rᵤ,g0=𝟏,θ=𝟏,h=𝘩)

Constructs new `UnitSystem` from `molarmass` constant, `vacuumpermeability`, `molargas` constant, `gravity` force reference, `angle` scale, and `planck` constant.

Examples include `SI2019`, `Metric`, `SI2019Engineering`, `MetricEngineering`, `SI1976`, `MetricDegree`, `MetricGradian`. In addition, the `ConventionalSystem` constructor further builds on `MetricSystem`, resulting in variations.
"""
MetricSystem(Mu=Mᵤ,μ0=μ₀,Ru=Rᵤ,g0=𝟏,θ=𝟏,h=𝘩,me=αinv^2*R∞*𝟐*h/𝘤) = UnitSystem(Ru*me/Mu/μₑᵤ/g0,h/τ/g0,𝘤,μ0,me,Mu,Kcd*(mₑ/me)^2*(h/𝘩)*g0,θ,𝟏,𝟏,g0,Universe,τ,𝟐,𝟑,𝟓,𝟕,𝟏𝟏,𝟏𝟗,𝟒𝟑)

"""
    ConventionalSystem(RK,KJ,Ru=Rᵤ,g0=𝟏) = MetricSystem(milli,𝟐*RK/𝘤*α,Ru,g0,𝟐^2/RK/KJ^2)

Constructs new `UnitSystem` from von `klitzing` constant and `josephson` constant, with an optional specification of `universal` gas constant and `gravity` reference constant.

Examples include `Conventional` (based on 1990) and `CODATA` (based on 2014).
"""
ConventionalSystem(klitz,joseph,Ru=Rᵤ,g0=𝟏,θ=𝟏) = MetricSystem(milli,𝟐*klitz/𝘤*α,Ru,g0,θ,(𝟐*𝟐)/klitz/(joseph*joseph))

"""
    RankineSystem(U::UnitSystem,l,m,g0=𝟏) = EntropySystem(U,𝟏,l,m,°R,vacuumpermeability(U)/m/l/g0,kilo*molarmass(U),g0)

Constructs new `UnitSystem` from `U` rescaled along `length` and `mass` with optional `gravity` reference constant used to define technical and engineering units.

Examples: `FPS`, `British`, `IPS`, `English`, `Survey`.
"""
RankineSystem(u,l,m,g0=𝟏) = EntropySystem(u,𝟏,l,m,°R,UnitSystems.vacuumpermeability(u)/(m*l)/g0,UnitSystems.unit(kilo*UnitSystems.molarmass(u)),g0)

# historical units

const SI2019 = Quantity(MetricSystem())
const Metric = Quantity(MetricSystem(milli,τ/𝟐^6/𝟓^7))
#const SI2019Engineering = Quantity(MetricSystem(Mᵤ,μ₀/g₀,Rᵤ,g₀))
const MetricEngineering = Quantity(MetricSystem(milli,τ/𝟐^6/𝟓^7/g₀,Rᵤ,g₀))
const MetricTurn = Quantity(MetricSystem(milli,τ/𝟐^6/𝟓^7,Rᵤ,𝟏,𝟏/τ))
const MetricDegree = Quantity(MetricSystem(milli,τ/𝟐^6/𝟓^7,Rᵤ,𝟏,𝟐^3*𝟑^2*𝟓/τ))
const MetricArcminute = Quantity(MetricSystem(milli,τ/𝟐^6/𝟓^7,Rᵤ,𝟏,𝟐^5*𝟑^3*𝟓^2/τ))
const MetricArcsecond = Quantity(MetricSystem(milli,τ/𝟐^6/𝟓^7,Rᵤ,𝟏,𝟐^7*𝟑^4*𝟓^3/τ))
const MetricGradian = Quantity(MetricSystem(milli,τ/𝟐^6/𝟓^7,Rᵤ,𝟏,𝟐^4*𝟓^2/τ))
const SI1976 = Quantity(MetricSystem(milli,τ/𝟐^6/𝟓^7,Constant(8.31432)))
const CODATA = Quantity(ConventionalSystem(RK2014,KJ2014,Rᵤ2014))
const Conventional = Quantity(ConventionalSystem(RK1990,KJ1990))
const International = Quantity(ElectricSystem(Metric,Ωᵢₜ,Vᵢₜ))
const InternationalMean = Quantity(ElectricSystem(Metric,Constant(1.00049),Constant(1.00034)))

const EMU = Quantity(GaussSystem(Metric,𝟏,𝟐*τ))
const ESU = Quantity(GaussSystem(Metric,(hecto*𝘤)^-2,𝟐*τ))
const Gauss = Quantity(GaussSystem(Metric,𝟏,𝟐*τ,centi/𝘤))
const LorentzHeaviside = Quantity(GaussSystem(Metric,𝟏,𝟏,centi/𝘤))
#const Thomson = Quantity(GaussSystem(Metric,𝟏,𝟐*τ,inv(𝟐)))
const Kennelly = Quantity(GaussSystem(Metric,inv(𝟏𝟎^7),𝟐*τ,𝟏,𝟏,𝟏))

const British = Quantity(RankineSystem(Metric,ft,lb*g₀/ft))
#const British2019 = Quantity(RankineSystem(SI2019,ft,lb*g₀/ft))
const Survey = Quantity(RankineSystem(Metric,ftUS,lb,g₀/ftUS))
#const Survey2019 = Quantity(RankineSystem(SI2019,ftUS,lb,g₀/ftUS))
const English = Quantity(RankineSystem(Metric,ft,lb,g₀/ft))
#const English2019 = Quantity(RankineSystem(SI2019,ft,lb,g₀/ft))
const FPS = Quantity(RankineSystem(Metric,ft,lb))
#const FPS2019 = Quantity(RankineSystem(SI2019,ft,lb))
const IPS = Quantity(RankineSystem(Metric,ft/𝟐^2/𝟑,lb*g₀*𝟐^2*𝟑/ft))
#const IPS2019 = Quantity(RankineSystem(SI2019,ft/𝟐^2/𝟑,lb*g₀*𝟐^2*𝟑/ft))

#const Astronomical = Quantity(EntropySystem(Metric,𝟏,𝟏,𝟏/G))
const Hubble = Quantity(AstronomicalSystem(Metric,th,𝘤*th,mₑ))#ħ/th/𝘤^2/mP^2))
const Cosmological = Quantity(AstronomicalSystem(Metric,lc/𝘤,lc,mc))
const CosmologicalQuantum = Quantity(AstronomicalSystem(Metric,tcq,lcq,mcq))
#const EMU2019 = Quantity(EntropySystem(SI2019,𝟏,centi,milli))
#const ESU2019 = Quantity(EntropySystem(SI2019,𝟏,centi,milli,𝟏,kilo*μ₀/𝘤^2))
#const Mixed = Quantity(EntropySystem(Metric,𝟏,𝟏,𝟏,𝟏,μ₀))
const Nautical = Quantity(EntropySystem(Metric,HOUR,nm,em^3,𝟏,τ*𝟑^3/𝟐^10/𝟓^12,milli))
const Meridian = Quantity(EntropySystem(Metric,𝟏,em,em^3,𝟏,τ/𝟐^6/𝟓^7,milli))
#const MeridianEngineering = Quantity(EntropySystem(Metric,𝟏,em,em^3,𝟏,τ/𝟐^6/𝟓^7*em/g₀,milli,g₀/em))
#const GravitationalSI2019 = Quantity(EntropySystem(SI2019,𝟏,𝟏,g₀))
const GravitationalMetric = Quantity(EntropySystem(Metric,𝟏,𝟏,g₀))
#const GravitationalMeridian = Quantity(EntropySystem(Metric,𝟏,em,g₀*em^2,𝟏,τ/𝟐^6/𝟓^7*em/g₀,milli))
const IAU☉ = Quantity(EntropySystem(Metric,DAY,au,GM☉/G))
const IAUE = Quantity(EntropySystem(Metric,DAY,LD,GME/G))
const IAUJ = Quantity(EntropySystem(Metric,DAY,JD,GMJ/G))
const MTS = Quantity(EntropySystem(Metric,𝟏,𝟏,kilo))
const KKH = Quantity(EntropySystem(Metric,HOUR,kilo,𝟏))
const MPH = Quantity(EntropySystem(FPS,HOUR,mi,𝟏))
const FFF = Quantity(EntropySystem(Metric,𝟕*𝟐*DAY,fur,(𝟐*𝟑^2*𝟓)*lb,°R,Constant(0.),𝟏))

# natural units

const Planck = Quantity(UnitSystem(𝟏,𝟏,𝟏,𝟏,√(𝟐*τ*αG)))
const PlanckGauss = Quantity(UnitSystem(𝟏,𝟏,𝟏,𝟐*τ,√αG))
const Stoney = Quantity(UnitSystem(𝟏,αinv,𝟏,𝟐*τ,√(αG*αinv)))
const Hartree = Quantity(UnitSystem(𝟏,𝟏,αinv,𝟐*τ*α^2,𝟏))
const Rydberg = Quantity(UnitSystem(𝟏,𝟏,𝟐*αinv,τ/𝟐*α^2,inv(𝟐)))
const Schrodinger = Quantity(UnitSystem(𝟏,𝟏,αinv,𝟐*τ*α^2,√(αG*αinv)))
const Electronic = Quantity(UnitSystem(𝟏,αinv,𝟏,𝟐*τ,𝟏))
const Natural = Quantity(UnitSystem(𝟏,𝟏,𝟏,𝟏,𝟏))
const NaturalGauss = Quantity(UnitSystem(𝟏,𝟏,𝟏,𝟐*τ,𝟏))
const QCD = Quantity(UnitSystem(𝟏,𝟏,𝟏,𝟏,inv(μₚₑ)))
const QCDGauss = Quantity(UnitSystem(𝟏,𝟏,𝟏,𝟐*τ,inv(μₚₑ)))
const QCDoriginal = Quantity(UnitSystem(𝟏,𝟏,𝟏,𝟐*τ*α,inv(μₚₑ)))

export SI, MKS, ME, GM, CGS, CGS2019, CGSm, CGSe, HLU, FFF, AE, EE, BG # SIE, GSI, GSI2019
export EnglishEngineering, BritishGravitational, AbsoluteEnglish, EnglishUS, EE2019, IAU
#const SIE, GSI2019, GSI = SI2019Engineering, GravitationalSI2019, GravitationalSI2019
const SI, MKS, ME, GM, IAU = SI2019, Metric, MetricEngineering, GravitationalMetric, IAU☉
const CGS, CGSm, CGSe, HLU = Gauss, EMU, ESU, LorentzHeaviside
const EnglishEngineering, BritishGravitational, BG = English, British, British
const EnglishUS, AbsoluteEnglish, AE, EE = Survey, FPS, FPS, English
#const EnglishEngineering2019,BritishGravitational2019,BG2019 = English2019,British2019,British2019
#const EE2019,AbsoluteEnglish2019,AE2019 = English2019,FPS2019,FPS2019

for u ∈ Systems
    @eval unitname(::typeof(normal($u))) = $(string(u))
end
