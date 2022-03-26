
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

const centi,milli,HOUR = inv((𝟐*𝟓)^2),inv((𝟐*𝟓)^3),𝟐^4*𝟑^2*𝟓^2
const fur,slug,slugUS,°R,K = 𝟐^2*𝟑*𝟓*𝟏𝟏*ft,lb*g₀/ft,lb*g₀/ftUS,𝟓/𝟑^2,𝟑^2/𝟓
const mₑ,μ₀ = αinv^2*R∞*𝟐*𝘩/𝘤,𝟐*𝘩/𝘤*α/𝘦^2 # ≈ 4π*(1e-7+5.5e-17), exact charge
const ħ,μₚₑ,Rᵤ,αL,αG,Mᵤ = 𝘩/τ,μₚᵤ/μₑᵤ,NA*kB,centi/𝘤,(mₑ/mP)^2,NA*mₑ/μₑᵤ
const pc,G,DAY,nm = au*𝟐^7*𝟑^4*𝟓^3/τ,𝘤*ħ/mP^2,𝟐^7*𝟑^3*𝟓^2,𝟐^4*𝟓^5/𝟑^3 # nautical mile
const GM☉ =au^3*k^2/DAY^2; const th = 𝟏𝟎^3*pc/H0; const ΛC = 𝟑*ΩΛ*(th*𝘤)^-2
const lc,mc,ρΛ = 𝟐*sqrt(τ/ΛC),𝘤^2/(𝟐*sqrt(τ*ΛC*G)),ΛC*𝘤^4/(𝟐^2*τ)/G
const lcq,mcq = sqrt.(sqrt.((𝘤*ħ/ρΛ,ρΛ*ħ^3/𝘤^5)))
const tcq = lcq*sqrt(mcq/sqrt(sqrt(ρΛ*(𝘤*ħ)^3)))

const Universe = Coupling(αG,α,μₑᵤ,μₚᵤ,ΩΛ)

Constant(::UnitSystem{a,b,c,d,e,f,g}) where {a,b,c,d,e,f,g} = UnitSystem(Constant(a),Constant(b),Constant(c),Constant(d),Constant(e),Constant(f),Constant(g[1]),Constant(g[2]),Constant(g[3]),Constant(g[4]),Constant(g[5]),Universe,Constant(g[7]),𝟐,𝟑,𝟓,𝟕,𝟏𝟏,𝟏𝟗,𝟒𝟑)

export MetricSystem, ConventionalSystem, RankineSystem
export AstronomicalSystem, ElectricSystem, GaussSystem, EntropySystem

"""
    MetricSystem(Mu=Mᵤ,μ0=μ₀,Ru=Rᵤ,g0=𝟏,h=𝘩)

Constructs new `UnitSystem` from `molarmass` constant, `vacuumpermeability`, `universal` gas constant, `gravity` force reference, and `planck` constant.

Examples include `SI2019`, `Metric`, `SI2019Engineering`, `MetricEngineering`, `SI1976`. In addition, the `ConventionalSystem` constructor further builds on `MetricSystem`, resulting in variations.
"""
MetricSystem(Mu=Mᵤ,μ0=μ₀,Ru=Rᵤ,g0=𝟏,h=𝘩,me=αinv^2*R∞*𝟐*h/𝘤/g0) = UnitSystem(Ru*me/Mu/μₑᵤ/g0,h/τ/g0,𝘤,μ0,me,Mu,Kcd*(mₑ/me)^2*(h/𝘩)*g0,𝟏,𝟏,𝟏,g0,Universe,τ,𝟐,𝟑,𝟓,𝟕,𝟏𝟏,𝟏𝟗,𝟒𝟑)

"""
    ConventionalSystem(RK,KJ,Ru=Rᵤ,g0=𝟏) = MetricSystem(0.001,𝟐*RK/𝘤*α,Ru,g0,𝟒/RK/KJ^2)

Constructs new `UnitSystem` from von `klitzing` constant and `josephson` constant, with an optional specification of `universal` gas constant and `gravity` reference constant.

Examples include `Conventional` (based on 1990) and `CODATA` (based on 2014).
"""
ConventionalSystem(klitz,joseph,Ru=Rᵤ,g0=𝟏) = MetricSystem(milli,𝟐*klitz/𝘤*α,Ru,g0,(𝟐*𝟐)/klitz/(joseph*joseph))

"""
    RankineSystem(U::UnitSystem,l,m,g0=𝟏) = EntropySystem(U,𝟏,l,m,°R,𝟐*τ,1000molarmass(U),g0)

Constructs new `UnitSystem` from `U` rescaled along `length` and `mass` with optional `gravity` reference constant used to define technical and engineering units.

Examples: `British`, `British2019`, `Survey`, `Survey2019`, `English`, `English2019`, `FPS`, `FPS2019`.
"""
RankineSystem(u,l,m,g0=𝟏) = EntropySystem(u,𝟏,l,m,°R,𝟐*τ,UnitSystems.unit(𝟏𝟎^3*UnitSystems.molarmass(u)),g0)

# historical units

const SI2019 = Quantity(MetricSystem())
const Metric = Quantity(MetricSystem(milli,𝟐*τ/𝟏𝟎^7))
const SI2019Engineering = Quantity(MetricSystem(Mᵤ,μ₀,Rᵤ,g₀))
const MetricEngineering = Quantity(MetricSystem(milli,𝟐*τ/𝟏𝟎^7,Rᵤ,g₀))
const SI1976 = Quantity(MetricSystem(milli,𝟐*τ/𝟏𝟎^7,Constant(8.31432)))
const CODATA = Quantity(ConventionalSystem(RK2014,KJ2014,Rᵤ2014))
const Conventional = Quantity(ConventionalSystem(RK1990,KJ1990))
const International = Quantity(ElectricSystem(Metric,Ωᵢₜ,Vᵢₜ))
const InternationalMean = Quantity(ElectricSystem(Metric,Constant(1.00049),Constant(1.00034)))

const EMU = Quantity(GaussSystem(Metric,𝟏,𝟐*τ))
const ESU = Quantity(GaussSystem(Metric,(𝟏𝟎^2*𝘤)^-2,𝟐*τ))
const Gauss = Quantity(GaussSystem(Metric,𝟏,𝟐*τ,centi/𝘤))
const LorentzHeaviside = Quantity(GaussSystem(Metric,𝟏,𝟏,centi/𝘤))
const Thomson = Quantity(GaussSystem(Metric,𝟏,𝟐*τ,inv(𝟐)))
const Kennelly = Quantity(GaussSystem(Metric,inv(𝟏𝟎^7),𝟐*τ,𝟏,𝟏,𝟏))

const British = Quantity(RankineSystem(Metric,ft,slug))
const British2019 = Quantity(RankineSystem(SI2019,ft,slug))
const Survey = Quantity(RankineSystem(Metric,ftUS,lb,g₀/ftUS))
const Survey2019 = Quantity(RankineSystem(SI2019,ftUS,lb,g₀/ftUS))
const English = Quantity(RankineSystem(Metric,ft,lb,g₀/ft))
const English2019 = Quantity(RankineSystem(SI2019,ft,lb,g₀/ft))
const FPS = Quantity(RankineSystem(Metric,ft,lb))
const FPS2019 = Quantity(RankineSystem(SI2019,ft,lb))

const Astronomical = Quantity(AstronomicalSystem(Metric))
const Hubble = Quantity(EntropySystem(Metric,th,𝘤*th,𝟏))
const Cosmological = Quantity(EntropySystem(Metric,lc/𝘤,lc,mc))
const CosmologicalQuantum = Quantity(EntropySystem(Metric,tcq,lcq,mcq))
const EMU2019 = Quantity(EntropySystem(SI2019,𝟏,centi,milli))
const ESU2019 = Quantity(EntropySystem(SI2019,𝟏,centi,milli,𝟏,𝟏𝟎^3*μ₀/𝘤^2))
const Mixed = Quantity(EntropySystem(Metric,𝟏,𝟏,𝟏,𝟏,μ₀))
const GravitationalSI2019 = Quantity(EntropySystem(SI2019,𝟏,𝟏,g₀))
const GravitationalMetric = Quantity(EntropySystem(Metric,𝟏,𝟏,g₀))
const IAU☉ = Quantity(EntropySystem(Metric,DAY,au,GM☉/G))
const IAUE = Quantity(EntropySystem(Metric,DAY,au,GME/G))
const IAUJ = Quantity(EntropySystem(Metric,DAY,au,GMJ/G))
const MTS = Quantity(EntropySystem(Metric,𝟏,𝟏,𝟏𝟎^3))
const KKH = Quantity(EntropySystem(Metric,HOUR,𝟏𝟎^3,𝟏))
const MPH = Quantity(EntropySystem(English,HOUR,𝟐^5*𝟑*𝟓*𝟏𝟏,𝟏))
const Nautical = Quantity(EntropySystem(English,HOUR,𝟐^6*𝟓*𝟏𝟗,𝟏))
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

export SI, MKS, SIE, ME, GSI2019, GSI, GM, CGS, CGS2019, CGSm, CGSe, HLU, FFF, AE, EE, BG
export EnglishEngineering, BritishGravitational, AbsoluteEnglish, EnglishUS, EE2019, IAU
const SI, MKS, SIE, ME, IAU = SI2019, Metric, SI2019Engineering, MetricEngineering, IAU☉
const GSI2019, GSI, GM = GravitationalSI2019, GravitationalSI2019, GravitationalMetric
const CGS, CGS2019, CGSm, CGSe, HLU = Gauss, EMU2019, EMU, ESU, LorentzHeaviside
const EnglishEngineering, BritishGravitational, BG = English, British, British
const EnglishUS, AbsoluteEnglish, AE, EE, EE2019 = Survey, FPS, FPS, English, English2019
const EnglishEngineering2019,BritishGravitational2019,BG2019 = English2019,British2019,British2019
const AbsoluteEnglish2019,AE2019 = FPS2019,FPS2019

for u ∈ Systems
    @eval unitname(::typeof(normal($u))) = $(string(u))
end
