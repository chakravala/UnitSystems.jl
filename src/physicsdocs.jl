
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
    μₑᵤ, μₚᵤ, μₚₑ, αinv, αG, ΩΛ

Physical measured dimensionless `Coupling` values with uncertainty are the electron to proton mass ratio `μₑᵤ`, proton to atomic mass ratio `μₚᵤ`, proton to electron mass ratio `μₚₑ`, inverted fine structure constant `αinv`, and the gravitaional coupling constant `αG`.

```Julia
julia> μₑᵤ # electronunit(Universe)
$μₑᵤ

julia> μₚᵤ # protonunit(Universe)
$μₚᵤ

julia> μₚₑ # protonelectron(Universe)
$μₚₑ

julia> αinv # 1/finestructure(Universe)
$αinv

julia> αG # coupling(Universe)
$αG

julia> ΩΛ # darkenergydensity(Universe)
$ΩΛ
```
""" Universe, μₑᵤ, μₚᵤ, μₚₑ, μₑₚ, αinv, αG, meu, mpu, mpe, ainv, aG, ΩΛ, electronunit, protonunit, protonelectron, finestructure, coupling, darkenergydensity

@doc """
$(unitext(:turn,"2π/angle(U)"))

Complete rotation `angle` of revolution from a full circle.
```Julia
julia> turn(MetricEngineering) # rad
$(turn(MetricEngineering))

julia> turn(MetricDegree) # deg
$(turn(MetricDegree))

julia> turn(MetricArcminute) # amin
$(turn(MetricArcminute))

julia> turn(MetricArcsecond) # asec
$(turn(MetricArcsecond))

julia> turn(MetricGradian) # gon
$(turn(MetricGradian))
```
""" turn

@doc """
$(unitext(:spat,"4π/solidangle(U)"))

Complete spherical `solidangle` around point from a full sphere.
```Julia
julia> spat(MetricEngineering) # rad²
$(spat(MetricEngineering))

julia> spat(MetricDegree) # deg²
$(spat(MetricDegree))

julia> spat(MetricArcminute) # amin²
$(spat(MetricArcminute))

julia> spat(MetricArcsecond) # asec²
$(spat(MetricArcsecond))

julia> spat(MetricGradian) # gon²
$(spat(MetricGradian))
```
""" spat

@doc """
$(unitext(:luminousefficacy,"Kcd*power(U)\nluminousefficacy(U::UnitSystem{𝟏}) = 𝟏"))

Luminous efficacy of monochromatic radiation `Kcd` of frequency 540 THz (lm⋅W⁻¹).
```Julia
julia> luminousefficacy(Metric) # lm⋅W⁻¹
$(luminousefficacy(Metric))

julia> luminousefficacy(CODATA) # lm⋅W⁻¹
$(luminousefficacy(CODATA))

julia> luminousefficacy(Conventional) # lm⋅W⁻¹
$(luminousefficacy(Conventional))

julia> luminousefficacy(International) # lm⋅W⁻¹
$(luminousefficacy(International))

julia> luminousefficacy(British) # lm⋅s³⋅slug⋅ft⁻²
$(luminousefficacy(British))
```
""" luminousefficacy, Kcd

@doc """
$(unitext(:molarmass,"avogadro(U)*electronmass(U)/electronunit(U)"))

Molar mass constant `Mᵤ` is the ratio of the `molarmass` and `relativemass` of a chemical.
```Julia
julia> molarmass(CGS) # g⋅mol⁻¹
$(molarmass(CGS))

julia> molarmass(Metric) # kg⋅mol⁻¹
$(molarmass(Metric))

julia> molarmass(SI2019) # kg⋅mol⁻¹
$(molarmass(SI2019))

julia> molarmass(International) # kg⋅mol⁻¹
$(molarmass(International))
```
""" molarmass, Mᵤ, Mu

@doc """
$(unitext(:avogadro, "molargas(x)/boltzmann(x) # Mᵤ/dalton(x)"))

Avogadro `NA` is `molarmass(x)/dalton(x)` number of atoms in a 12 g sample of C₁₂.
```Julia
julia> avogadro(SI2019) # mol⁻¹
$(avogadro(SI2019))

julia> avogadro(Metric) # mol⁻¹
$(avogadro(Metric))

julia> avogadro(CODATA) # mol⁻¹
$(avogadro(CODATA))

julia> avogadro(Conventional) # mol⁻¹
$(avogadro(Conventional))

julia> avogadro(English) # lb-mol⁻¹
$(avogadro(English))

julia> avogadro(British) # slug-mol⁻¹
$(avogadro(British))
```
""" avogadro, NA

@doc """
$(unitext(:planckreduced,"planck(x)/turn(x)"))

Reduced Planck constant `ħ` is a Planck per radian (J⋅s⋅rad⁻¹ or ft⋅lb⋅s⋅rad⁻¹).

```Julia
julia> planckreduced(SI2019) # J⋅s⋅rad⁻¹
$(planckreduced(SI2019))

julia> planckreduced(SI2019)*lightspeed(SI2019) # J⋅m⋅rad⁻¹
$(planckreduced(SI2019)*lightspeed(SI2019))

julia> planckreduced(CODATA) # J⋅s⋅rad⁻¹
$(planckreduced(CODATA))

julia> planckreduced(Conventional) # J⋅s⋅rad⁻¹
$(planckreduced(Conventional))

julia> planckreduced(SI2019)/elementarycharge(SI2019) # eV⋅s⋅rad⁻¹
$(planckreduced(SI2019)/elementarycharge(SI2019))

julia> planckreduced(SI2019)*lightspeed(SI2019)/elementarycharge(SI2019) # eV⋅m⋅rad⁻¹
$(planckreduced(SI2019)*lightspeed(SI2019)/elementarycharge(SI2019))

julia> planckreduced(British) # ft⋅lb⋅s⋅rad⁻¹
$(planckreduced(British))
```
""" planckreduced, ħ

@doc """
$(unitext(:planck,"turn(x)*planckreduced(x)"))

Planck constant `𝘩` is energy per electromagnetic frequency (J⋅s or ft⋅lb⋅s).

```Julia
julia> planck(SI2019) # J⋅s
$(planck(SI2019))

julia> planck(SI2019)*lightspeed(SI2019) # J⋅m
$(planck(SI2019)*lightspeed(SI2019))

julia> planck(CODATA) # J⋅s
$(planck(CODATA))

julia> planck(Conventional) # J⋅s
$(planck(Conventional))

julia> planck(SI2019)/elementarycharge(SI2019) # eV⋅s
$(planck(SI2019)/elementarycharge(SI2019))

julia> planck(SI2019)*lightspeed(SI2019)/elementarycharge(SI2019) # eV⋅m
$(planck(SI2019)*lightspeed(SI2019)/elementarycharge(SI2019))

julia> planck(British) # ft⋅lb⋅s
$(planck(British))
```
""" planck, 𝘩, hh

@doc """
$(unitext(:boltzmann,"molargas(x)/avogadro(x)"))

Boltzmann constant `kB` is the entropy amount of a unit number microstate permutation.

```Julia
julia> boltzmann(SI2019) # J⋅K⁻¹
$(boltzmann(SI2019))

julia> boltzmann(Metric) # J⋅K⁻¹
$(boltzmann(Metric))

julia> boltzmann(SI2019)/elementarycharge(SI2019) # eV⋅K⁻¹
$(boltzmann(SI2019)/elementarycharge(SI2019))

julia> boltzmann(SI2019)/planck(SI2019) # Hz⋅K⁻¹
$(boltzmann(SI2019)/planck(SI2019))

julia> boltzmann(CGS) # erg⋅K⁻¹
$(boltzmann(CGS))

julia> boltzmann(SI2019)/calorie(SI2019) # calᵢₜ⋅K⁻¹
$(boltzmann(SI2019)/calorie(SI2019))

julia> boltzmann(SI2019)*°R/calorie(SI2019) # calᵢₜ⋅°R⁻¹
$(boltzmann(SI2019)*°R/calorie(SI2019))

julia> boltzmann(Brtish) # ft⋅lb⋅°R⁻¹
$(boltzmann(British))

julia> boltzmann(SI2019)/planck(SI2019)/lightspeed(SI2019) # m⁻¹⋅K⁻¹
$(boltzmann(SI2019)/planck(SI2019)/lightspeed(SI2019))

julia> avogadro(SI2019)*boltzmann(SI2019)/calorie(SI2019) # calᵢₜ⋅mol⁻¹⋅K⁻¹
$(avogadro(SI2019)*boltzmann(SI2019)/calorie(SI2019))

julia> dB(boltzmann(SI2019)) # dB(W⋅K⁻¹⋅Hz⁻¹)
$(dB(boltzmann(SI2019)))
```
""" boltzmann, kB

@doc """
$(unitext(:lightspeed,"𝟏/sqrt(vacuumpermeability(U)*vacuumpermittivity(U))/lorentz(U)"))

Speed of light in a vacuum `𝘤` for massless particles (m⋅s⁻¹ or ft⋅s⁻¹).

```Julia
julia> lightspeed(Metric) # m⋅s⁻¹
$(lightspeed(Metric))

julia> lightspeed(English) # ft⋅s⁻¹
$(lightspeed(English))

julia> lightspeed(IAU) # au⋅D⁻¹
$(lightspeed(IAU))
```
""" lightspeed, 𝘤, cc

@doc """
$(unitext(:vacuumpermeability,"𝟏/vacuumpermittivity(U)/(lightspeed(U)*lorentz(U))^2"))

Magnetic permeability in a classical vacuum defined as `μ₀` in SI units (H⋅m⁻¹, kg⋅m²⋅C⁻²).

```Julia
julia> vacuumpermeability(Metric) # H⋅m⁻¹
$(vacuumpermeability(Metric))

julia> vacuumpermeability(Conventional) # H⋅m⁻¹
$(vacuumpermeability(Conventional))

julia> vacuumpermeability(CODATA) # H⋅m⁻¹
$(vacuumpermeability(CODATA))

julia> vacuumpermeability(SI2019) # H⋅m⁻¹
$(vacuumpermeability(SI2019))

julia> vacuumpermeability(International) # H⋅m⁻¹
$(vacuumpermeability(International))

julia> vacuumpermeability(EMU) # abH⋅cm⁻¹
$(vacuumpermeability(EMU))

julia> vacuumpermeability(ESU) # statH⋅cm⁻¹
$(vacuumpermeability(ESU))
```
""" vacuumpermeability, μ₀, m0

@doc """
$(unitext(:lorentz,"spat(U)*biotsavart(U)/vacuumpermeability(U)/rationalization(U)"))

Electromagnetic proportionality constant `αL` for the Lorentz's law force (dimensionless).

```Julia
julia> lorentz(Metric)
$(lorentz(Metric))

julia> lorentz(LorentzHeaviside)
$(lorentz(LorentzHeaviside))

julia> lorentz(Gauss)
$(lorentz(Gauss))
```
""" lorentz, αL, aL, C

@doc """
$(unitext(:rationalization,"spat(U)*biotsavart(U)/vacuumpermeability(U)/lorentz(U)"))

Constant of magnetization and polarization density or `spat(U)*coulomb(U)*permittivity(U)`.

```Julia
julia> rationalization(Metric)
$(rationalization(Metric))

julia> rationalization(Gauss)
$(rationalization(Gauss))
```
""" rationalization, Λ

@doc """
$(unitext(:electronmass,"protonmass(U)/protonelectron(U) # αinv^2*R∞*2𝘩/𝘤"))

Electron rest mass `mₑ` of subatomic particle with `-𝘦` elementary charge  (kg or slugs).
```Julia
julia> electronmass(Metric) # kg
$(electronmass(Metric))

julia> electronmass(CODATA) # kg
$(electronmass(CODATA))

julia> electronmass(Conventional) # kg
$(electronmass(Conventional))

julia> electronmass(International) # kg
$(electronmass(International))

julia> electronmass(Metric)/dalton(Metric) # Da
$(electronmass(Metric)/dalton(Metric))

julia> electronmass(QCD) # mₚ
$(electronmass(QCD))

julia> electronmass(Hartree) # mₑ
$(electronmass(Hartree))

julia> electronmass(Metric)*lightspeed(Metric)^2 # J
$(electronmass(Metric)*lightspeed(Metric)^2)

julia> electronmass(SI2019)*lightspeed(SI2019)^2/elementarycharge(SI2019) # eV⋅𝘤⁻²
$(electronmass(SI2019)*lightspeed(SI2019)^2/elementarycharge(SI2019))

julia> electronmass(English) # lb
$(electronmass(English))
```
""" electronmass, mₑ, me

@doc """
$(unitext(:dalton,"molarmass(U)/avogadro(U)"))

Atomic mass unit `Da` of 1/12 of the C₁₂ carbon-12 atom's mass  (kg or slugs).
```Julia
julia> dalton(Metric) # kg
$(dalton(Metric))

julia> dalton(Hartree) # mₑ
$(dalton(Hartree))

julia> dalton(QCD) # mₚ
$(dalton(QCD))

julia> dalton(Metric)*lightspeed(Metric)^2 # J
$(dalton(Metric)*lightspeed(Metric)^2)

julia> dalton(SI2019)*lightspeed(SI2019)^2/elementarycharge(SI2019) # eV⋅𝘤⁻²
$(dalton(SI2019)*lightspeed(SI2019)^2/elementarycharge(SI2019))

julia> dalton(British) # lb
$(dalton(British))
```
""" dalton, Da, mu, mᵤ, atomicmass

@doc """
$(unitext(:protonmass,"protonunit(U)*dalton(U)"))

Proton mass `mₚ` of subatomic particle with `+𝘦` elementary charge  (kg or mass).
```Julia
julia> protonmass(Metric) # kg
$(protonmass(Metric))

julia> protonmass(SI2019)*lightspeed(SI2019)^2/elementarycharge(SI2019) # eV⋅𝘤⁻²
$(protonmass(SI2019)*lightspeed(SI2019)^2/elementarycharge(SI2019))

julia> protonmass(Metric)/dalton(Metric) # Da
$(protonmass(Metric)/dalton(Metric))

julia> protonmass(Hartree) # mₑ
$(protonmass(Hartree))

julia> protonmass(QCD) # mₚ
$(protonmass(QCD))
```
""" protonmass, mₚ, mp

@doc """
$(unitext(:planckmass,"electronmass(U)/sqrt(coupling(U))"))

Planck mass factor `mP` from the gravitational coupling constant `αG` (kg or slugs).
```Julia
juila> planckmass(Metric)*lightspeed(Metric)^2/elementarycharge(Metric) # eV⋅𝘤⁻²
$(planckmass(Metric)*lightspeed(Metric)^2/elementarycharge(Metric))

juila> planckmass(Metric) # kg
$(planckmass(Metric))

juila> planckmass(Metric)/dalton(Metric) # Da
$(planckmass(Metric)/dalton(Metric))

juila> planckmass(Metric)*lightspeed(Metric)^2/elementarycharge(Metric)/sqrt(𝟐^2*τ) # eV⋅𝘤⁻²
$(planckmass(Metric)*lightspeed(Metric)^2/elementarycharge(Metric)/sqrt(𝟐^2*τ))

julia> planckmass(PlanckGauss) # mP
$(planckmass(PlanckGauss))
```
""" planckmass, mP

@doc """
$(unitext(:gravitation,"lightspeed(U)*planckreduced(U)/planckmass(U)^2"))

Universal gravitational constant `G` of Newton's law (m³⋅kg⁻¹⋅s⁻² or ft³⋅slug⁻¹⋅s⁻²).
```Julia
juila> gravitation(Metric) # m³⋅kg⁻¹⋅s⁻²
$(gravitation(Metric))

julia> gravitation(English) # ft³⋅lbm⁻¹⋅s⁻²
$(gravitation(English))

julia> gravitation(PlanckGauss)
$(gravitation(PlanckGauss))
```
""" gravitation, G, GG

@doc """
$(unitext(:einstein,"𝟐^2*τ*gravitation(U)/lightspeed(U)^4"))

Einstein's gravitational constant from the Einstein field equations (s⋅²⋅m⁻¹⋅kg⁻¹).
```Julia
julia> einstein(Metric) # s²⋅m⁻¹⋅kg⁻¹
$(einstein(Metric))

julia> einstein(IAU) # day²⋅au⁻¹⋅M☉⁻¹
$(einstein(IAU))
```
""" einstein, κ

#=@doc """
$(unitext(:einstein2,"𝟐*spat(U)*gravitation(U)/lightspeed(U)^2"))

Einstein's gravitational constant from the Einstein field equations (m⋅kg⁻¹).
```Julia
julia> einstein2(Metric) # m⋅kg⁻¹
$(einstein2(Metric))

julia> einstein2(IAU) # au⋅M☉⁻¹
$(einstein2(IAU))
```
""" einstein2=#

@doc """
$(unitext(:gravity,"# mass*acceleration/force"))

Gravitational force reference used in technical engineering units (kg⋅m⋅N⁻¹⋅s⁻²).
```Julia
julia> gravity(Metric)
$(gravity(Metric))

julia> gravity(MetricEngineering) # m⋅kg⋅N⁻¹⋅s⁻²
$(gravity(MetricEngineering))

julia> gravity(English) # ft⋅lbm⋅lbf⁻¹⋅s⁻²
$(gravity(English))
```
""" gravity

@doc """
$(unitext(:molargas,"boltzmann(x)*avogadro(x)"))

Universal gas constant `Rᵤ` is factored into specific `gasconstant(x)*molarmass(x)` values.

```Julia
julia> molargas(SI2019) # J⋅K⁻¹⋅mol⁻¹
$(molargas(SI2019))

julia> molargas(English)/𝟐^4/𝟑^2 # psi⋅ft³⋅°R⁻¹⋅lb-mol⁻¹
$(molargas(English)/𝟐^4/𝟑^2)

julia> molargas(English)/atmosphere(English) # atm⋅ft³⋅R⁻¹⋅lb-mol⁻¹
$(molargas(English)/atmosphere(English))

julia> molargas(English)/thermalunit(English) # BTU⋅°R⁻¹⋅lb-mol⁻¹
$(molargas(English)/thermalunit(English))

julia> molargas(Metric)/atmosphere(Metric) # atm⋅m³⋅K⁻¹⋅mol⁻¹
$(molargas(Metric)/atmosphere(Metric))

julia> molargas(Metric)/torr(Metric) # m³⋅torr⋅K⁻¹⋅mol⁻¹
$(molargas(Metric)/torr(Metric))

julia> molargas(English)/torr(English) # ft³⋅torr⋅°R⁻¹⋅lb-mol⁻¹
$(molargas(English)/torr(English))

julia> molargas(CGS) # erg⋅K⁻¹⋅mol⁻¹
$(molargas(CGS))

julia> molargas(English) # ft⋅lb⋅°R⁻¹⋅lb-mol⁻¹
$(molargas(English))

julia> molargas(British) # ft⋅lb⋅°R⁻¹⋅slug-mol⁻¹
$(molargas(British))

julia> molargas(SI1976) # J⋅K⁻¹⋅mol⁻¹ (US1976 Standard Atmosphere)
$(molargas(SI1976))
```
""" molargas, Rᵤ, Ru

@doc """
$(unitext(:stefan,"τ^5/𝟐^4*boltzmann(U)^4/(𝟑*𝟓*planck(U)^3*lightspeed(U)^2)"))

Stefan-Boltzmann proportionality `σ` of black body radiation (W⋅m⁻²⋅K⁻⁴ or ?⋅ft⁻²⋅°R⁻⁴).

```Julia
julia> stefan(SI2019) # W⋅m⁻²⋅K⁻⁴
$(stefan(SI2019))

julia> stefan(Metric) # W⋅m⁻²⋅K⁻⁴
$(stefan(Metric))

julia> stefan(Conventional) # W⋅m⁻²⋅K⁻⁴
$(stefan(Conventional))

julia> stefan(CODATA) # W⋅m⁻²⋅K⁻⁴
$(stefan(CODATA))

julia> stefan(Metric)*day(Metric)/(calorie(Metric)*100^2) # cal⋅cm⁻²⋅day⁻¹⋅K⁻⁴
$(stefan(Metric)*day(Metric)/calorie(Metric))

julia> stefan(English) # lb⋅s⁻¹⋅ft⁻³⋅°R⁻⁴
$(stefan(English))
```
""" stefan, σ, SB

@doc """
$(unitext(:radiationdensity,"𝟐^2*stefan(U)/lightspeed(U)"))

Raditation density constant of black body radiation (J⋅m⁻³⋅K⁻⁴ or lb⋅ft⁻²⋅°R⁻⁴).

```Julia
julia> radiationdensity(Metric) # J⋅m⁻³⋅K⁻⁴
$(radiationdensity(Metric))

julia> radiationdensity(SI2019) # J⋅m⁻³⋅K⁻⁴
$(radiationdensity(SI2019))

julia> radiationdensity(Conventional) # J⋅m⁻³⋅K⁻⁴
$(radiationdensity(Conventional))

julia> radiationdensity(CODATA) # J⋅m⁻³⋅K⁻⁴
$(radiationdensity(CODATA))

julia> radiationdensity(International) # J⋅m⁻³⋅K⁻⁴
$(radiationdensity(International))
```
""" radiationdensity

@doc """
$(unitext(:vacuumpermittivity,"𝟏/vacuumpermeability(U)/(lightspeed(U)*lorentz(U))^2"))

Dielectric permittivity constant `ε₀` of a classical vacuum (C²⋅N⁻¹⋅m⁻²).

```Julia
julia> vacuumpermittivity(Metric) # F⋅m⁻¹
$(vacuumpermittivity(Metric))

julia> vacuumpermittivity(Conventional) # F⋅m⁻¹
$(vacuumpermittivity(Conventional))

julia> vacuumpermittivity(CODATA) # F⋅m⁻¹
$(vacuumpermittivity(CODATA))

julia> vacuumpermittivity(SI2019) # F⋅m⁻¹
$(vacuumpermittivity(SI2019))

julia> vacuumpermittivity(International) # F⋅m⁻¹
$(vacuumpermittivity(International))

julia> vacuumpermittivity(EMU) # abF⋅cm⁻¹
$(vacuumpermittivity(EMU))

julia> vacuumpermittivity(ESU) # statF⋅cm⁻¹
$(vacuumpermittivity(ESU))

julia> vacuumpermittivity(SI2019)/elementarycharge(SI2019) # 𝘦²⋅eV⁻¹⋅m⁻¹
$(vacuumpermittivity(SI2019)/elementarycharge(SI2019))
```
""" vacuumpermittivity, ε₀, ϵ₀, e0

@doc """
$(unitext(:electrostatic,"rationalization(U)/𝟐/τ/vacuumpermittivity(U)"))

Electrostatic proportionality constant `kₑ` for the Coulomb's law force (N⋅m²⋅C⁻²).

```Julia
julia> electrostatic(Metric) # N⋅m²⋅C⁻²
$(electrostatic(Metric))

julia> electrostatic(CODATA) # N·m²⋅C⁻²
$(electrostatic(CODATA))

julia> electrostatic(SI2019) # N·m²⋅C⁻²
$(electrostatic(SI2019))

julia> electrostatic(Conventional) # N·m²⋅C⁻²
$(electrostatic(Conventional))

julia> electrostatic(International) # N·m²⋅C⁻²
$(electrostatic(International))

julia> electrostatic(EMU) # dyn⋅cm²⋅abC⁻²
$(electrostatic(EMU))

julia> electrostatic(ESU) # dyn⋅cm²⋅statC⁻²
$(electrostatic(ESU))

julia> electrostatic(HLU) # dyn⋅cm²⋅hlC⁻²
$(electrostatic(HLU))
```
""" electrostatic, kₑ, ke

@doc """
$(unitext(:biotsavart,"vacuumpermeability(U)*lorentz(U)*rationalization(U)/𝟐/τ"))

Magnetostatic proportionality constant `αB` for the Biot-Savart's law (H/m).

```Julia
julia> biotsavart(Metric) # H⋅m⁻¹
$(biotsavart(Metric))

julia> biotsavart(CODATA) # H⋅m⁻¹
$(biotsavart(CODATA))

julia> biotsavart(SI2019) # H⋅m⁻¹
$(biotsavart(SI2019))

julia> biotsavart(Conventional) # H⋅m⁻¹
$(biotsavart(Conventional))

julia> biotsavart(International) # H⋅m⁻¹
$(biotsavart(International))

julia> biotsavart(InternationalMean) # H⋅m⁻¹
$(biotsavart(InternationalMean))

julia> biotsavart(EMU) # abH⋅cm⁻¹
$(biotsavart(EMU))

julia> biotsavart(ESU) # statH⋅cm⁻¹
$(biotsavart(ESU))

julia> biotsavart(Gauss) # abH⋅cm⁻¹
$(biotsavart(Gauss))

julia> biotsavart(HLU) # hlH⋅cm⁻¹
$(biotsavart(HLU))
```
""" biotsavart, αB, aB

@doc """
$(unitext(:magnetostatic,"lorentz(U)*biotsavart(U) # electrostatic(U)/lightspeed(U)^2"))

Magnetic proportionality constant `kₘ` for the Ampere's law force (N·s²⋅C⁻²).

```Julia
julia> magnetostatic(Metric) # H⋅m⁻¹
$(magnetostatic(Metric))

julia> magnetostatic(CODATA) # H⋅m⁻¹
$(magnetostatic(CODATA))

julia> magnetostatic(SI2019) # H⋅m⁻¹
$(magnetostatic(SI2019))

julia> magnetostatic(Conventional) # H⋅m⁻¹
$(magnetostatic(Conventional))

julia> magnetostatic(International) # H⋅m⁻¹
$(magnetostatic(International))

julia> magnetostatic(EMU) # abH⋅m⁻¹
$(magnetostatic(EMU))

julia> magnetostatic(ESU) # statH⋅m⁻¹
$(magnetostatic(ESU))

julia> magnetostatic(HLU) # hlH⋅m⁻¹
$(magnetostatic(HLU))
```
""" magnetostatic, kₘ, km

@doc """
$(unitext(:vacuumimpedance,"vacuumpermeability(U)*lightspeed(U)*rationalization(U)*lorentz(U)^2"))

Vacuum impedance of free space `Z₀` is magnitude ratio of electric to magnetic field (Ω).
```Julia
julia> vacuumimpedance(Metric) # Ω
$(vacuumimpedance(Metric))

julia> vacuumimpedance(Conventional) # Ω
$(vacuumimpedance(Conventional))

julia> vacuumimpedance(CODATA) # Ω
$(vacuumimpedance(CODATA))

julia> vacuumimpedance(SI2019) # Ω
$(vacuumimpedance(SI2019))

julia> vacuumimpedance(International) # Ω
$(vacuumimpedance(International))

julia> vacuumimpedance(InternationalMean) # Ω
$(vacuumimpedance(InternationalMean))

julia> 120π # 3e8*μ₀ # Ω
$(120π)

julia> vacuumimpedance(EMU) # abΩ
$(vacuumimpedance(EMU))

julia> vacuumimpedance(ESU) # statΩ
$(vacuumimpedance(ESU))

julia> vacuumimpedance(HLU) # hlΩ
$(vacuumimpedance(HLU))

julia> vacuumimpedance(IPS) # in⋅lb⋅s⋅C⁻²
$(vacuumimpedance(IPS))
```
""" vacuumimpedance, Z₀, Z0

@doc """
$(unitext(:elementarycharge,"√(𝟐*planck(U)*finestructure(U)/vacuumimpedance(U))"))

Quantized elementary charge `𝘦` of a proton or electron `2/(klitzing(U)*josephson(U))` (C).
```Julia
julia> elementarycharge(SI2019) # C
$(elementarycharge(SI2019))

julia> elementarycharge(Metric) # C
$(elementarycharge(Metric))

julia> elementarycharge(CODATA) # C
$(elementarycharge(CODATA))

julia> elementarycharge(Conventional) # C
$(elementarycharge(Conventional))

julia> elementarycharge(International) # C
$(elementarycharge(International))

julia> elementarycharge(EMU) # abC
$(elementarycharge(EMU))

julia> elementarycharge(ESU) # statC
$(elementarycharge(ESU))

julia> elementarycharge(Hartree) # 𝘦
$(elementarycharge(Hartree))
```
""" elementarycharge, 𝘦, ee

@doc """
$(unitext(:faraday,"elementarycharge(U)*avogadro(U)"))

Electric charge per mole of electrons `𝔉` based on elementary charge (C⋅mol⁻¹).
```Julia
julia> faraday(SI2019) # C⋅mol⁻¹
$(faraday(SI2019))

julia> faraday(Metric) # C⋅mol⁻¹
$(faraday(Metric))

julia> faraday(CODATA) # C⋅mol⁻¹
$(faraday(CODATA))

julia> faraday(Conventional) # C⋅mol⁻¹
$(faraday(Conventional))

julia> faraday(International) # C⋅mol⁻¹
$(faraday(International))

julia> faraday(InternationalMean) # C⋅mol⁻¹
$(faraday(InternationalMean))

julia> faraday(EMU) # abC⋅mol⁻¹
$(faraday(EMU))

julia> faraday(ESU) # statC⋅mol⁻¹
$(faraday(ESU))

julia> faraday(Metric)/kilocalorie(Metric) # kcal⋅(V-g-e)⁻¹
$(faraday(Metric)/kilocalorie(Metric))

julia> faraday(Metric)/3600 # A⋅h⋅mol⁻¹
$(faraday(Metric)/HOUR)
```
""" faraday, 𝔉, FF

@doc """
$(unitext(:josephson,"𝟐*elementarycharge(U)*lorentz(U)/planck(U)"))

Josephson constant `KJ` relating potential difference to irradiation frequency (Hz⋅V⁻¹).
```Julia
julia> josephson(SI2019) # Hz⋅V⁻¹
$(josephson(SI2019))

julia> josephson(Metric) # Hz⋅V⁻¹
$(josephson(Metric))

julia> josephson(Conventional) # Hz⋅V⁻¹
$(josephson(Conventional))

julia> josephson(CODATA) # Hz⋅V⁻¹
$(josephson(CODATA))

julia> josephson(International) # Hz⋅V⁻¹
$(josephson(International))

julia> josephson(EMU) # Hz⋅abV⁻¹
$(josephson(EMU))

julia> josephson(ESU) # Hz⋅statV⁻¹
$(josephson(ESU))
```
""" josephson, KJ

@doc """
$(unitext(:magneticfluxquantum,"planck(U)/𝟐/elementarycharge(U)/lorentz(U)"))

Magnetic flux quantum `Φ₀` is `𝟏/josephson(U)` (Wb).
```Julia
julia> magneticfluxquantum(SI2019) # Wb
$(magneticfluxquantum(SI2019))

julia> magneticfluxquantum(Metric) # Wb
$(magneticfluxquantum(Metric))

julia> magneticfluxquantum(Conventional) # Wb
$(magneticfluxquantum(Conventional))

julia> magneticfluxquantum(International) # Wb
$(magneticfluxquantum(International))

julia> magneticfluxquantum(InternationalMean) # Wb
$(magneticfluxquantum(InternationalMean))

julia> magneticfluxquantum(EMU) # Mx
$(magneticfluxquantum(EMU))

julia> magneticfluxquantum(ESU) # statWb
$(magneticfluxquantum(ESU))
```
""" magneticfluxquantum, Φ₀

@doc """
$(unitext(:klitzing,"planck(U)/elementarycharge(U)^2"))

Quantized Hall resistance `RK` (Ω).
```Julia
julia> klitzing(SI2019) # Ω
$(klitzing(SI2019))

julia> klitzing(Metric) # Ω
$(klitzing(Metric))

julia> klitzing(Conventional) # Ω
$(klitzing(Conventional))

julia> klitzing(International) # Ω
$(klitzing(International))

julia> klitzing(CODATA) # Ω
$(klitzing(CODATA))

julia> klitzing(EMU) # abΩ
$(klitzing(EMU))

julia> klitzing(ESU) # statΩ
$(klitzing(ESU))
```
""" klitzing, RK

@doc """
$(unitext(:conductancequantum,"𝟐*elementarycharge(U)^2/planck(U) # 2/klitzing(U)"))

Conductance quantum `G₀` is a quantized unit of electrical conductance (S).
```Julia
julia> conductancequantum(SI2019) # S
$(conductancequantum(SI2019))

julia> conductancequantum(Metric) # S
$(conductancequantum(Metric))

julia> conductancequantum(Conventional) # S
$(conductancequantum(Conventional))

julia> conductancequantum(CODATA) # S
$(conductancequantum(CODATA))

julia> conductancequantum(International) # S
$(conductancequantum(International))

julia> conductancequantum(InternationalMean) # S
$(conductancequantum(InternationalMean))

julia> conductancequantum(EMU) # abS
$(conductancequantum(EMU))

julia> conductancequantum(ESU) # statS
$(conductancequantum(ESU))
```
""" conductancequantum, G₀, G0

@doc """
$(unitext(:hartree,"electronmass(U)/gravity(U)*(lightspeed(U)*finestructure(U))^2"))

Hartree electric potential energy `Eₕ` of the hydrogen atom at ground state is `2R∞*𝘩*𝘤` (J).
```Julia
julia> hartree(SI2019)/elementarycharge(SI2019) # eV
$(hartree(SI2019)/elementarycharge(SI2019))

julia> hartree(Metric) # J
$(hartree(Metric))

julia> hartree(CGS) # erg
$(hartree(CGS))

julia> hartree(Metric)*avogadro(Metric)/kilo # kJ⋅mol⁻¹
$(hartree(Metric)*avogadro(Metric)/kilo)

julia> hartree(Metric)*avogadro(Metric)/kilocalorie(Metric) # kcal⋅mol⁻¹
$(hartree(Metric)*avogadro(Metric)/kilocalorie(Metric))

julia> 𝟐*rydberg(CGS) # Eₕ/𝘩/𝘤/100 cm⁻¹
$(𝟐*rydberg(CGS))

julia> hartree(Metric)/planck(Metric) # Hz
$(hartree(Metric)/planck(Metric))

julia> hartree(Metric)/boltzmann(Metric) # K
$(hartree(Metric)/boltzmann(Metric))
```
In a Gaussian unit system where `4π*ε₀ == 1` the Hartree energy is `𝘦^2/a₀`.
""" hartree, Eₕ, Eh

@doc """
$(unitext(:rydberg,"hartree(U)/2planck(U)/lightspeed(U) # Eₕ/2𝘩/𝘤"))

Rydberg constant `R∞` is lowest energy photon capable of ionizing atom at ground state (m⁻¹).
```Julia
julia> rydberg(Metric) # m⁻¹
$(rydberg(Metric))
```
The Rydberg constant for hydrogen `RH` is `R∞*mₚ/(mₑ+mₚ)` (m⁻¹).
```Julia
julia> rydberg(Metric)*protonmass(Metric)/(electronmass(Metric)+protonmass(Metric)) # m⁻¹
$(rydberg(Metric)*protonmass(Metric)/(electronmass(Metric)+protonmass(Metric)))
```
Rydberg unit of photon energy `Ry` is `𝘩*𝘤*R∞` or `Eₕ/2` (J).
```Julia
julia> hartree(Metric)/2 # J
$(hartree(Metric)/𝟐)

julia> hartree(SI2019)/𝟐/elementarycharge(SI2019) # eV
$(hartree(SI2019)/𝟐/elementarycharge(SI2019))
```
Rydberg photon frequency `𝘤*R∞` or `Eₕ/2𝘩` (Hz).
```Julia
julia> lightspeed(Metric)*rydberg(Metric) # Hz
$(lightspeed(Metric)*rydberg(Metric))
```
Rydberg wavelength `1/R∞` (m).
```Julia
julia> 𝟏/rydberg(Metric) # m
$(𝟏/rydberg(Metric))

julia> 𝟏/rydberg(Metric)/τ # m⋅rad⁻¹
$(𝟏/rydberg(Metric)/τ)
```
Precision measurements of the Rydberg constants are within a relative standard uncertainty of under 2 parts in 10¹², and is chosen to constrain values of other physical constants.
""" rydberg, R∞, RH, Ry

@doc """
$(unitext(:bohr,"planckreduced(U)*gravity(U)/electronmass(U)/lightspeed(U)/finestructure(U)"))

Bohr radius of the hydrogen atom in its ground state `a₀` (m).
```Julia
julia> bohr(Metric) # m
$(bohr(Metric))

julia> bohr(IPS) # in
$(bohr(IPS))

julia> bohr(Hartree) # a₀
$(bohr(Hartree))
```
""" bohr, a₀, a0
#julia> bohr(Metric)/length(PlanckGauss) # ℓP
#$(bohr(Metric)/length(PlanckGauss))

#=@doc """
$(unitext(:bohrreduced,"bohr(U)*(1+1/protonelectron(U))"))

Reduced Bohr radius including the effect of reduced mass in hydrogen atom (m).
```Julia
julia> bohrreduced(Metric) # m
$(bohrreduced(Metric))

julia> bohrreduced(Metric)/bohr(Metric) # a₀
$(bohrreduced(Metric)/bohr(Metric))
```
""" bohrreduced=#

@doc """
$(unitext(:electronradius,"finestructure(U)*planckreduced(U)*gravity(U)/electronmass(U)/lightspeed(U)"))

Classical electron radius or Lorentz radius or Thomson scattering length (m).
```Julia
julia> electronradius(Metric) # m
$(electronradius(Metric))

julia> electronradius(CODATA) # m
$(electronradius(CODATA))

julia> electronradius(Conventional) # m
$(electronradius(Conventional))

julia> electronradius(Hartree) # a₀
$(electronradius(Hartree))
```
""" electronradius, rₑ, re

@doc """
$(unitext(:magneton,"elementarycharge(U)*planckreduced(U)*lorentz(U)/2electronmass(U)"))

Bohr magneton `μB` natural unit for expressing magnetic moment of electron (J⋅T⁻¹).
```Julia
julia> magneton(SI2019) # J⋅T⁻¹
$(magneton(SI2019))

julia> magneton(Metric) # J⋅T⁻¹
$(magneton(Metric))

julia> magneton(CODATA) # J⋅T⁻¹
$(magneton(CODATA))

julia> magneton(Conventional) # J⋅T⁻¹
$(magneton(Conventional))

julia> magneton(International) # J⋅T⁻¹
$(magneton(International))

julia> magneton(ESU) # statA⋅cm²
$(magneton(ESU))

julia> magneton(SI2019)/elementarycharge(SI2019) # eV⋅T⁻¹
$(magneton(SI2019)/elementarycharge(SI2019))

julia> magneton(Hartree) # 𝘤⋅ħ⋅mₑ⁻¹
$(magneton(Hartree))
```
""" magneton, μB

@doc """
$(unitext(:wienwavelength,"planck(U)*lightspeed(U)/boltzmann(U)/(𝟓+W₀(-𝟓*exp(-𝟓)))"))

Wien wavelength displacement law constant based on Lambert `W₀` evaluation (m⋅K or ft⋅°R).
```Julia
julia> wienwavelength(Metric) # m⋅K
$(wienwavelength(Metric))

julia> wienwavelength(SI2019) # m⋅K
$(wienwavelength(SI2019))

julia> wienwavelength(Conventional) # m⋅K
$(wienwavelength(Conventional))

julia> wienwavelength(CODATA) # m⋅K
$(wienwavelength(CODATA))

julia> wienwavelength(English) # ft⋅°R
$(wienwavelength(English))
```
""" wienwavelength

@doc """
$(unitext(:wienfrequency,"(𝟑+W₀(-𝟑*exp(-𝟑)))*boltzmann(U)/planck(U)"))

Wien frequency radiation law constant based on Lambert `W₀` evaluation (Hz⋅K⁻¹).
```Julia
julia> wienfrequency(Metric) # Hz⋅K⁻¹
$(wienfrequency(Metric))

julia> wienfrequency(SI2019) # Hz⋅K⁻¹
$(wienfrequency(SI2019))

julia> wienfrequency(Conventional) # Hz⋅K⁻¹
$(wienfrequency(Conventional))

julia> wienfrequency(CODATA) # Hz⋅K⁻¹
$(wienfrequency(CODATA))

julia> wienfrequency(English) # Hz⋅°R⁻¹
$(wienfrequency(English))
```
""" wienfrequency

@doc """
$(unitext(:loschmidt,"atmosphere(U)/boltzmann(U)/temperature(T₀,SI2019,U)"))

Number of molecules (number density) of an ideal gas in a unit volume (m⁻³ or ft⁻³).
```Julia
julia> loschmidt(SI2019) # m⁻³
$(loschmidt(SI2019))

julia> loschmidt(Metric,atm,T₀) # m⁻³
$(loschmidt(Metric,atm,T₀))

julia> loschmidt(Conventional,atm,T₀) # m⁻³
$(loschmidt(Conventional,atm,T₀))

julia> loschmidt(CODATA,atm,T₀) # m⁻³
$(loschmidt(CODATA,atm,T₀))

julia> loschmidt(SI1976,atm,T₀) # m⁻³
$(loschmidt(SI1976,atm,T₀))

julia> loschmidt(English) # ft⁻³
$(loschmidt(English))

julia> loschmidt(IAU) # au⁻³
$(loschmidt(IAU))
```
""" loschmidt

@doc """
$(unitext(:amagat,"loschmidt(U)/avogadro(U)"))

Number of moles of an ideal gas in a unit volume (mol⋅m⁻³ or lb-mol⋅ft⁻³).
```Julia
julia> amagat(Metric) # mol⋅m⁻³
$(amagat(Metric))

julia> amagat(SI2019) # mol⋅m⁻³
$(amagat(SI2019))

julia> amagat(English) # slug-mol⋅ft⁻³
$(amagat(English))
```
""" amagat

@doc """
$(unitext(:hyperfine,"frequency($ΔνCs,U)"))

Unperturbed groundstate hyperfine transition frequency `ΔνCs` of caesium-133 atom (Hz).
```Julia
julia> hyperfine(Metric) # Hz
$(hyperfine(Metric))
```
""" hyperfine, ΔνCs

@doc """
$(unitext(:hubble,"time(U,Hubble)"))

Hubble universe expansion frequency parameter.
```Julia
julia> hubble(Metric)
$(hubble(Metric))

julia> hubble(Hubble)
$(hubble(Hubble))

julia> hubble(Cosmological)
$(hubble(Cosmological))

julia> 𝟏/hubble(Metric)/year(Metric)
$(𝟏/hubble(Metric)/year(Metric))
```
""" hubble

@doc """
$(unitext(:cosmological,"𝟑*darkenergydensity(U)*(hubble(U)/lightspeed(U))^2"))

Cosmological constant from Einstein's controversial theory expanded on by Hubble.
```Julia
julia> cosmological(Metric)
$(cosmological(Metric))

julia> cosmological(Hubble)
$(cosmological(Hubble))

julia> cosmological(Cosmological)
$(cosmological(Cosmological))
```
""" cosmological

@doc """
$(unitext(:solarmass,"mass($(GM☉/G),U)"))

Solar `mass` estimated from gravitational constant estimates (kg or slug).
```Julia
julia> solarmass(Metric) # kg
$(solarmass(Metric))

julia> solarmass(British) # slug
$(solarmass(British))

julia> solarmass(English) # lb
$(solarmass(English))

julia> solarmass(IAUE) # ME
$(solarmass(IAUE))

julia> solarmass(IAUJ) # MJ
$(solarmass(IAUJ))

julia> solarmass(QCD) # mₚ
$(solarmass(QCD))

julia> solarmass(Metric)/dalton(Metric) # Da
$(solarmass(Metric)/dalton(Metric))
```
""" solarmass, mₛ

@doc """
$(unitext(:earthmass,"mass($(GME/G),U)"))

Earth `mass` estimated from gravitational constant estimates (kg or slug).
```Julia
julia> earthmass(Metric) # kg
$(earthmass(Metric))

julia> earthmass(British) # slug
$(earthmass(British))

julia> earthmass(English) # lb
$(earthmass(English))

julia> earthmass(IAU) # M☉
$(earthmass(IAU))

julia> earthmass(IAUJ) # MJ
$(earthmass(IAUJ))

julia> earthmass(QCD) # mₚ
$(earthmass(QCD))

julia> earthmass(Metric)/dalton(Metric) # Da
$(earthmass(Metric)/dalton(Metric))
```
""" earthmass

@doc """
$(unitext(:jupitermass,"mass($(GMJ/G),U)"))

Jupiter `mass` estimated from gravitational constant estimates (kg or slug).
```Julia
julia> jupitermass(Metric) # kg
$(jupitermass(Metric))

julia> jupitermass(British) # slug
$(jupitermass(British))

julia> jupitermass(English) # lb
$(jupitermass(English))

julia> jupitermass(IAU) # M☉
$(jupitermass(IAU))

julia> jupitermass(IAUE) # ME
$(jupitermass(IAUE))

julia> jupitermass(QCD) # mₚ
$(jupitermass(QCD))

julia> jupitermass(Metric)/dalton(Metric) # Da
$(jupitermass(Metric)/dalton(Metric))
```
""" jupitermass

@doc """
$(unitext(:lunarmass,"earthmass(U)/μE☾"))

Lunar `mass` estimated from `μE☾` Earth-Moon mass ratio (kg or slug).
```Julia
julia> lunarmass(Metric) # kg
$(lunarmass(Metric))

julia> lunarmass(British) # slug
$(lunarmass(British))

julia> lunarmass(English) # lb
$(lunarmass(English))

julia> lunarmass(IAU) # M☉
$(lunarmass(IAU))

julia> lunarmass(IAUE) # ME
$(lunarmass(IAUE))

julia> lunarmass(IAUJ) # MJ
$(lunarmass(IAUJ))

julia> lunarmass(QCD) # mₚ
$(lunarmass(QCD))

julia> lunarmass(Metric)/dalton(Metric) # Da
$(lunarmass(Metric)/dalton(Metric))
```
""" lunarmass

@doc """
$(unitext(:gaussgravitation,"sqrt(gravitation(U)*solarmass(U)/astronomicalunit(U)^3)"))

Gaussian  gravitational constant `k` of Newton's laws (Hz or rad⋅D⁻¹).
```Julia
julia> gaussgravitation(MetricEngineering)
$(gaussgravitation(MetricEngineering))

julia> gaussgravitation(MetricGradian)
$(gaussgravitation(MetricGradian))

julia> gaussgravitation(MetricDegree)
$(gaussgravitation(MetricDegree))

julia> gaussgravitation(MetricArcminute)
$(gaussgravitation(MetricArcminute))

julia> gaussgravitation(MetricArcsecond)
$(gaussgravitation(MetricArcsecond))

juila> gaussgravitation(MPH)
$(gaussgravitation(MPH))

julia> gaussgravitation(IAU)
$(gaussgravitation(IAU))
```
""" gaussgravitation, k, kG

@doc """
$(unitext(:gaussianyear,"turn(U)/gaussgravitation(U)"))

Orbit `time` defined by `gaussgravitation` constant `kG` for neglible `mass` satellite (s).
```Julia
julia> gaussianyear(Metric) # s
$(gaussianyear(Metric))

julia> gaussianyear(MPH) # h
$(gaussianyear(MPH))

julia> gaussianyear(IAU) # D
$(gaussianyear(IAU))
```
""" gaussianyear

@doc """
$(unitext(:siderealyear,"gaussianyear(U)/√(𝟏+earthmass(IAU)+lunarmass(IAU))"))

Orbit `time` defined by `gaussgravitation` constant `kG` and Earth-Moon system `mass` (s).
```Julia
julia> siderealyear(Metric) # s
$(siderealyear(Metric))

julia> siderealyear(MPH) # h
$(siderealyear(MPH))

julia> siderealyear(IAU) # D
$(siderealyear(IAU))
```
""" siderealyear

@doc """
$(unitext(:jovianyear,"τ*day(U)*√(jupiterdistance(U)^3/solarmass(U)/gravitation(U))/√(𝟏+jupitermass(IAU))"))

Orbit `time` defined by `jupiterdistance` and the Sun-Jupiter system `mass` (s).
```Julia
julia> jovianyear(Metric) # s
$(jovianyear(Metric))

julia> jovianyear(MPH) # h
$(jovianyear(MPH))

julia> jovianyear(IAU) # D
$(jovianyear(IAU))
```
""" jovianyear

@doc """
$(unitext(:gaussianmonth,"τ*sqrt(lunardistance(U)^3/earthmass(U)/gravitation(U))"))

Orbit `time` defined by `lunardistance` and `earthmass` for neglible `mass` satellite (s).
```Julia
julia> gaussianmonth(Metric) # s
$(gaussianmonth(Metric))

julia> gaussianmonth(MPH) # h
$(gaussianmonth(MPH))

julia> gaussianmonth(IAU) # D
$(gaussianmonth(IAU))
```
""" gaussianmonth

@doc """
$(unitext(:siderealmonth,"gaussianmonth(U)/√(𝟏+lunarmass(IAU))"))

Orbit `time` defined by standard `lunardistance` and the Earth-Moon system `mass` (s).
```Julia
julia> siderealmonth(Metric) # s
$(siderealmonth(Metric))

julia> siderealmonth(MPH) # h
$(siderealmonth(MPH))

julia> siderealmonth(IAU) # D
$(siderealmonth(IAU))
```
""" siderealmonth

@doc """
$(unitext(:synodicmonth,"𝟏/(𝟏/siderealmonth(U)-𝟏/siderealyear(U))"))

Orbit `time` defined by `siderealmonth` and `siderealyear` of Sun-Earth-Moon system (s).
```Julia
julia> synodicmonth(Metric) # s
$(synodicmonth(Metric))

julia> synodicmonth(MPH) # h
$(synodicmonth(MPH))

julia> synodicmonth(IAU) # D
$(synodicmonth(IAU))
```
""" synodicmonth

@doc """
$(unitext(:earthradius,"sqrt(earthmass(U)*gravitation(U)/gforce(U))"))

Approximate `length` of standard Earth two-body radius consistent with units (m or ft).
```Julia
julia> earthradius(KKH) # km
$(earthradius(KKH))

julia> earthradius(Nautical) # nm
$(earthradius(Nautical))

julia> earthradius(IAU) # au
$(earthradius(IAU))
```
""" earthradius

@doc """
$(unitext(:greatcircle,"τ*earthradius(U)"))

Approximate `length` of standard Earth two-body circle consistent with units (m or ft).
```Julia
julia> greatcircle(KKH) # km
$(greatcircle(KKH))

julia> greatcircle(Nautical) # nm
$(greatcircle(Nautical))

julia> greatcircle(IAU) # au
$(greatcircle(IAU))
```
""" greatcircle

@doc """
    sackurtetrode(U::UnitSystem,P=atm,T=𝟏,m=Da) = log(kB*T/P*sqrt(m*kB*T/τ/ħ^2)^3)+5/2
    dimensionless : [𝟙], [𝟙], [𝟙], [𝟙], [𝟙]
    log(FL⁻²Θ⁻⁵ᐟ²A³ᐟ²⋅(μₑᵤ⁻³ᐟ²atm⁻¹τ⁻³ᐟ²exp(2⁻¹5) = 0.594141574194(26)))

Ideal gas entropy density for pressure `P`, temperature `T`, atomic mass `m` (dimensionless).
```Julia
julia> sackurtetrode(Metric)
$(sackurtetrode(Metric))

julia> sackurtetrode(SI2019)
$(sackurtetrode(SI2019))

julia> sackurtetrode(Conventional)
$(sackurtetrode(Conventional))

julia> sackurtetrode(CODATA)
$(sackurtetrode(CODATA))

julia> sackurtetrode(SI2019,𝟏𝟎^5)
$(sackurtetrode(SI2019,𝟏𝟎^5))
```
""" sackurtetrode

@doc """
    mechanicalheat(U::UnitSystem) = molargas(U)/molargas(Metric)*calorie(Metric)
    energy : [FL], [FL], [ML²T⁻²], [ML²T⁻²], [ML²T⁻²]

Heat to raise 1 `mass` unit of water by 1 `temperature` unit, or $(normal(molargas(Metric)/calorie(Metric))) `mechanicalheat` per `molaramount` per `temperature` units (J or ft⋅lb).
```Julia
julia> mechanicalheat(Metric) # J
$(mechanicalheat(Metric))

julia> mechanicalheat(English) # ft⋅lb
$(mechanicalheat(English))

julia> mechanicalheat(British) # ft⋅lb
$(mechanicalheat(British))
```
""" mechanicalheat

@doc """
$(unitext(:eddington,"mass(𝟏,U,Cosmological)"))

Approximate number of protons in the `Universe` as estimated by Eddington (kg or lb).
```Julia
julia> 𝟐^2^2^3/α # mₚ
$(𝟐^2^2^3/α)

julia> eddington(QCD) # mₚ
$(eddington(QCD))

julia> eddington(Metric) # kg
$(eddington(Metric))

julia> eddington(IAU) # M☉
$(eddington(IAU))

julia> eddington(Cosmological)
$(eddington(Cosmological))
```
""" eddington

include("derivedocs.jl")
