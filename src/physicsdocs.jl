
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
""" Universe, μₑᵤ, μₚᵤ, μₚₑ, αinv, αG, meu, mpu, mpe, ainv, aG, ΩΛ, electronunit, protonunit, protonelectron, finestructure, coupling, darkenergydensity

@doc """
    turn(U::UnitSystem) = 2π/angle(U)

Complete rotation `angle` of revolution from a full circle.
```Julia
julia> turn(MetricEngineering)
$(turn(MetricEngineering))
```
""" turn

@doc """
    sphere(U::UnitSystem) = turn(U)/angle(U)

Complete spherical `solidangle` of from a full sphere.
```Julia
julia> sphere(MetricEngineering)
$(sphere(MetricEngineering))
```
""" sphere

@doc """
```Julia
luminousefficacy(U::UnitSystem{1}) = 1
luminousefficacy(U::UnitSystem) = $(Kcd)power(U)
```

Luminous efficacy of monochromatic radiation `Kcd` of frequency 540 THz (lm⋅W⁻¹).
```Julia
julia> luminousefficacy(Metric) # lm⋅W⁻¹
$(luminousefficacy(Metric))

julia> luminousefficacy(CODATA) # lm⋅W⁻¹
$(luminousefficacy(CODATA))

julia> luminousefficacy(Conventional) # lm⋅W⁻¹
$(luminousefficacy(Conventional))

julia> luminousefficacy(CGS) # lm⋅s⋅erg⁻¹
$(luminousefficacy(CGS))

julia> luminousefficacy(British) # lm⋅s³⋅slug⋅ft⁻²
$(luminousefficacy(British))
```
""" luminousefficacy, Kcd

@doc """
    molarmass(U::UnitSystem) = avogadro(U)*electronmass(U)/electronunit(U)

Molar mass constant `Mᵤ` is the ratio of the `molarmass` and `relativemass` of a chemical.
```Julia
julia> molarmass(CGS) # g⋅mol⁻¹
$(molarmass(CGS))

julia> molarmass(CGS2019) # g⋅mol⁻¹
$(molarmass(CGS2019))

julia> molarmass(Metric) # kg⋅mol⁻¹
$(molarmass(Metric))

julia> molarmass(SI2019) # kg⋅mol⁻¹
$(molarmass(SI2019))
```
""" molarmass, Mᵤ, Mu

@doc """
    avogadro(x) = universalgas(x)/boltzmann(x) # Mᵤ/atomicmass(x), Mᵤ ≈ 0.001-3.5e-13

Avogadro `NA` is `molarmass(x)/atomicmass(x)` number of atoms in a 12 g sample of C₁₂.
```Julia
julia> avogadro(SI2019) # mol⁻¹
$(avogadro(SI2019))

julia> avogadro(Metric) # mol⁻¹
$(avogadro(Metric))

julia> avogadro(English) # lb-mol⁻¹
$(avogadro(English))

julia> avogadro(British) # slug-mol⁻¹
$(avogadro(British))
```
""" avogadro, NA

@doc """
    planckreduced(x) = planck(x)/turn(x)

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
    planck(x) = turn(x)*planckreduced(x)

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
$(planck(English))
```
""" planck, 𝘩, hh

@doc """
    boltzmann(x) = universalgas(x)/avogadro(x)

Boltzmann constant `kB` is the entropy amount of a unit number microstate permutation.
```Julia
pressure*molecularmass == density*boltzmann*temperature
```
It satisfies the ideal gas law.

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
    lightspeed(U::UnitSystem) = 1/sqrt(vacuumpermeability(U)*vacuumpermittivity(U))/lorentz(U)

Speed of light in a vacuum `𝘤` for massless particles (m⋅s⁻¹ or ft⋅s⁻¹).

```Julia
julia> lightspeed(Metric) # m⋅s⁻¹
$(lightspeed(Metric))

julia> lightspeed(English) # ft⋅s⁻¹
$(lightspeed(English))
```
""" lightspeed, 𝘤, cc

@doc """
    vacuumpermeability(U::UnitSystem) = 1/vacuumpermittivity(U)/(lightspeed(U)*lorentz(U))^2

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

julia> vacuumpermeability(EMU) # abH⋅cm⁻¹
$(vacuumpermeability(EMU))

julia> vacuumpermeability(ESU) # statH⋅cm⁻¹
$(vacuumpermeability(ESU))
```
""" vacuumpermeability, μ₀, m0

@doc """
    lorentz(U::UnitSystem) = sphere(U)*biotsavart(U)/vacuumpermeability(U)/rationalization(U)

Electromagnetic proportionality constant `αL` for the Lorentz's law force (dimensionless).

```Julia
julia> lorentz(Metric)
$(lorentz(Metric))

julia> lorentz(Thomson)
$(lorentz(Thomson))

julia> lorentz(Gauss)
$(lorentz(Gauss))
```
""" lorentz, αL, aL, C

@doc """
    rationalization(U::UnitSystem) = sphere(U)*biotsavart(U)/vacuumpermeability(U)/lorentz(U)

Constant of magnetization and polarization density or `sphere(U)*coulomb(U)*permittivity(U)`.

```Julia
julia> rationalization(Metric)
$(rationalization(Metric))

julia> rationalization(Gauss)
$(rationalization(Gauss))
```
""" rationalization, Λ

@doc """
    electronmass(U::UnitSystem) = protonmass(U)/protonelectron(U) # αinv^2*R∞*2𝘩/𝘤

Electron rest mass `mₑ` of subatomic particle with `-𝘦` elementary charge  (kg or slugs).
```Julia
julia> electronmass(Metric) # kg
$(electronmass(Metric))

julia> electronmass(Metric)/atomicmass(Metric) # Da
$μₑᵤ

julia> electronmass(Metric)*lightspeed(Metric)^2 # J
$(electronmass(Metric)*lightspeed(Metric)^2)

julia> electronmass(SI2019)*lightspeed(SI2019)^2/elementarycharge(SI2019) # eV⋅𝘤⁻²
$(electronmass(SI2019)*lightspeed(SI2019)^2/elementarycharge(SI2019))

julia> electronmass(English) # lb
$(electronmass(English))
```
""" electronmass, mₑ, me

@doc """
    atomicmass(U::UnitSystem) = Mᵤ/avogadro(U) # $(molarmass(SI2019)) ≈ 0.001-3.5e-13

Atomic mass unit `mᵤ` of 1/12 of the C₁₂ carbon-12 atom's mass  (kg or slugs).
```Julia
julia> atomicmass(Metric) # kg
$(atomicmass(Metric))

julia> atomicmass(Metric)/electronmass(Metric) # mₑ
$(atomicmass(Metric)/electronmass(Metric))

julia> atomicmass(Metric)*lightspeed(Metric)^2 # J
$(atomicmass(Metric)*lightspeed(Metric)^2)

julia> atomicmass(SI2019)*lightspeed(SI2019)^2/elementarycharge(SI2019) # eV⋅𝘤⁻²
$(atomicmass(SI2019)*lightspeed(SI2019)^2/elementarycharge(SI2019))

julia> atomicmass(British) # lb
$(atomicmass(British))
```
""" atomicmass, mᵤ, mu

@doc """
    protonmass(U::UnitSystem) = protonunit(U)*atomicmass(U)

Proton mass `mₚ` of subatomic particle with `+𝘦` elementary charge  (kg or mass).
```Julia
julia> protonmass(Metric) # kg
$(protonmass(Metric))

julia> protonmass(SI2019)*lightspeed(SI2019)^2/elementarycharge(SI2019) # eV⋅𝘤⁻²
$(protonmass(SI2019)*lightspeed(SI2019)^2/elementarycharge(SI2019))

julia> protonmass(Metric)/atomicmass(Metric) # mᵤ
$(protonmass(Metric)/atomicmass(Metric))

julia> protonmass(Metric)/electronmass(Metric) # mₑ
$(protonmass(Metric)/electronmass(Metric))
```
""" protonmass, mₚ, mp

@doc """
    planckmass(U::UnitSystem) = electronmass(U)/sqrt(coupling(U))

Planck mass factor `mP` from the gravitational coupling constant `αG` (kg or slugs).
```Julia
juila> planckmass(Metric)*lightspeed(Metric)^2/elementarycharge(Metric) # eV⋅𝘤⁻²
$(planckmass(Metric)*lightspeed(Metric)^2/elementarycharge(Metric))

juila> planckmass(Metric) # kg
$(planckmass(Metric))

juila> planckmass(Metric)/atomicmass(Metric) # mᵤ
$(planckmass(Metric)/atomicmass(Metric))

juila> planckmass(Metric)*lightspeed(Metric)^2/elementarycharge(Metric)/sqrt(𝟐^2*τ) # eV⋅𝘤⁻²
$(planckmass(Metric)*lightspeed(Metric)^2/elementarycharge(Metric)/sqrt(𝟐^2*τ))

juila> planckmass(Metric)/sqrt(𝟐^2*τ) # kg
$(planckmass(Metric)/sqrt(𝟐^2*τ))
```
""" planckmass, mP

@doc """
    newton(U::UnitSystem) = lightspeed(U)*planckreduced(U)/planckmass(U)^2

Universal gravitational constant `G` of Newton's law (m³⋅kg⁻¹⋅s⁻² or ft³⋅slug⁻¹⋅s⁻²).
```Julia
juila> newton(Metric) # m³⋅kg⁻¹⋅s⁻²
$(newton(Metric))

julia> newton(English) # ft³⋅lbm⁻¹⋅s⁻²
$(newton(English))

julia> newton(IAU) # au³⋅M☉⁻¹⋅day⁻²
$(newton(IAU))

julia> newton(Astronomical) # N⋅s⁴⋅m⁻⁴
$(newton(Astronomical))

julia> newton(PlanckGauss)
$(newton(PlanckGauss))
```
""" newton, G, GG

@doc """
    einstein(U::UnitSystem) = 2sphere(U)*newton(U)/lightspeed(U)^4

Einstein's gravitational constant from the Einstein field equations (s⋅²⋅m⁻¹⋅kg⁻¹).
```Julia
julia> einstein(Metric) # s²⋅m⁻¹⋅kg⁻¹
$(einstein(Metric))

julia> einstein(IAU) # day²⋅au⁻¹⋅M☉⁻¹
$(einstein(IAU))
```
""" einstein, κ

@doc """
    einstein2(U::UnitSystem) = 2sphere(U)*newton(U)/lightspeed(U)^2

Einstein's gravitational constant from the Einstein field equations (m⋅kg⁻¹).
```Julia
julia> einstein2(Metric) # m⋅kg⁻¹
$(einstein2(Metric))

julia> einstein2(IAU) # au⋅M☉⁻¹
$(einstein2(IAU))
```
""" einstein2

@doc """
    gravity(U::UnitSystem) # mass*acceleration/force

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
    universalgas(x) = boltzmann(x)*avogadro(x)

Universal gas constant `Rᵤ` is factored into specific `gasconstant(x)*molarmass(x)` values.
```Julia
pressure*molarmass == density*universal*temperature
```
It satisfies the ideal gas law.

```Julia
julia> universalgas(SI2019) # J⋅K⁻¹⋅mol⁻¹
$(universalgas(SI2019))

julia> universalgas(English)/𝟐^4/𝟑^2 # psi⋅ft³⋅°R⁻¹⋅lb-mol⁻¹
$(universalgas(English)/𝟐^4/𝟑^2)

julia> universalgas(English)/standardpressure(English) # atm⋅ft³⋅R⁻¹⋅lb-mol⁻¹
$(universalgas(English)/standardpressure(English))

julia> universalgas(English)/thermalunit(English) # BTU⋅°R⁻¹⋅lb-mol⁻¹
$(universalgas(English)/thermalunit(English))

julia> universalgas(Metric)/calorie(Metric) # cal⋅K⁻¹⋅mol⁻¹
$(universalgas(Metric)/calorie(Metric))

julia> universalgas(Metric)/standardpressure(Metric) # atm⋅m³⋅K⁻¹⋅mol⁻¹
$(universalgas(Metric)/standardpressure(Metric))

julia> universalgas(Metric)/torr(Metric) # m³⋅torr⋅K⁻¹⋅mol⁻¹
$(universalgas(Metric)/torr(Metric))

julia> universalgas(English)/torr(English) # ft³⋅torr⋅°R⁻¹⋅lb-mol⁻¹
$(universalgas(English)/torr(English))

julia> universalgas(CGS) # erg⋅K⁻¹⋅mol⁻¹
$(universalgas(CGS))

julia> universalgas(British) # ft⋅lb⋅°R⁻¹⋅slug-mol⁻¹
$(universalgas(British))
```
The 1976 United States Standard Atmosphere used R* = 8.31432 exactly.
""" universalgas, Rᵤ, Ru

@doc """
    stefan(U::UnitSystem) = π^4/2*sphere(U)*boltzmann(U)^4/(15planck(U)^3*lightspeed(U)^2)

Stefan-Boltzmann proportionality `σ` of black body radiation (W⋅m⁻²⋅K⁻⁴ or ?⋅ft⁻²⋅°R⁻⁴).

```Julia
julia> stefan(Metric) # W⋅m⁻²⋅K⁻⁴
$(stefan(Metric))

julia> stefan(CGS) # erg⋅cm⁻²⋅s⁻¹⋅K⁻⁴
$(stefan(CGS))

julia> stefan(Metric)*day(Metric)/(calorie(Metric)*100^2) # cal⋅cm⁻²⋅day⁻¹⋅K⁻⁴
$(stefan(Metric)*day(Metric)/calorie(Metric))

julia> stefan(English) # lb⋅s⁻¹⋅ft⁻³⋅°R⁻⁴
$(stefan(English))
```
""" stefan, σ, SB

@doc """
    radiationdensity(U::UnitSystem) = 4stefan(U)/lightspeed(U)

Raditation density constant of black body radiation (J⋅m⁻³⋅K⁻⁴ or lb⋅ft⁻²⋅°R⁻⁴).

```Julia
julia> radiationdensity(Metric) # J⋅m⁻³⋅K⁻⁴
$(radiationdensity(Metric))

julia> radiationdensity(CGS) # erg⋅cm⁻³⋅K⁻⁴
$(radiationdensity(CGS))
```
""" radiationdensity

@doc """
    vacuumpermittivity(U::UnitSystem) = 1/vacuumpermeability(U)/(lightspeed(U)*lorentz(U))^2

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

julia> vacuumpermittivity(EMU) # abF⋅cm⁻¹
$(vacuumpermittivity(EMU))

julia> vacuumpermittivity(ESU) # statF⋅cm⁻¹
$(vacuumpermittivity(ESU))

julia> vacuumpermittivity(SI2019)/elementarycharge(SI2019) # 𝘦²⋅eV⁻¹⋅m⁻¹
$(vacuumpermittivity(SI2019)/elementarycharge(SI2019))
```
""" vacuumpermittivity, ε₀, ϵ₀, e0

@doc """
    coulomb(U::UnitSystem) = rationalization(U)/sphere(U)/vacuumpermittivity(U)

Electrostatic proportionality constant `kₑ` for the Coulomb's law force (N⋅m²⋅C⁻²).

```Julia
julia> coulomb(Metric) # N⋅m²⋅C⁻²
$(coulomb(Metric))

julia> coulomb(CODATA) # N·m²⋅C⁻²
$(coulomb(CODATA))

julia> coulomb(SI2019) # N·m²⋅C⁻²
$(coulomb(SI2019))

julia> coulomb(Conventional) # N·m²⋅C⁻²
$(coulomb(Conventional))

julia> coulomb(EMU) # dyn⋅cm²⋅abC⁻²
$(coulomb(EMU))

julia> coulomb(ESU) # dyn⋅cm²⋅statC⁻²
$(coulomb(ESU))

julia> coulomb(HLU) # dyn⋅cm²⋅hlC⁻²
$(coulomb(HLU))
```
""" coulomb, kₑ, ke

@doc """
    biotsavart(U::UnitSystem) = vacuumpermeability(U)*lorentz(U)*rationalization(U)/sphere(U)

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
    ampere(U::UnitSystem) = lorentz(U)*biotsavart(U) # coulomb(U)/lightspeed(U)^2

Magnetic proportionality constant `kₘ` for the Ampere's law force (N·s²⋅C⁻²).

```Julia
julia> ampere(Metric) # H⋅m⁻¹
$(ampere(Metric))

julia> ampere(CODATA) # H⋅m⁻¹
$(ampere(CODATA))

julia> ampere(SI2019) # H⋅m⁻¹
$(ampere(SI2019))

julia> ampere(Conventional) # H⋅m⁻¹
$(ampere(Conventional))

julia> ampere(EMU) # abH⋅m⁻¹
$(ampere(EMU))

julia> ampere(ESU) # statH⋅m⁻¹
$(ampere(ESU))

julia> ampere(HLU) # hlH⋅m⁻¹
$(ampere(HLU))
```
""" ampere, kₘ, km

@doc """
    vacuumimpedance(U::UnitSystem) = vacuumpermeability(U)*lightspeed(U)*rationalization(U)*lorentz(U)^2

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

julia> 120π # 3e8*μ₀ # Ω
$(120π)

julia> vacuumimpedance(EMU) # abΩ
$(vacuumimpedance(EMU))

julia> vacuumimpedance(ESU) # statΩ
$(vacuumimpedance(ESU))

julia> vacuumimpedance(HLU) # hlΩ
$(vacuumimpedance(HLU))
```
""" vacuumimpedance, Z₀, Z0

@doc """
    elementarycharge(U::UnitSystem) = √(2planck(U)*finestructure(U)/vacuumimpedance(U))

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

julia> elementarycharge(EMU) # abC
$(elementarycharge(EMU))

julia> elementarycharge(ESU) # statC
$(elementarycharge(ESU))

julia> elementarycharge(Planck) # sqrt(4π/αinv)
$(elementarycharge(Planck))
```
""" elementarycharge, 𝘦, ee

@doc """
    faraday(U::UnitSystem) = elementarycharge(U)*avogadro(U)

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
    josephson(U::UnitSystem) = 2elementarycharge(U)*lorentz(U)/planck(U) # 1/magneticflux(U)

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

julia> josephson(EMU) # Hz⋅abV⁻¹
$(josephson(EMU))

julia> josephson(ESU) # Hz⋅statV⁻¹
$(josephson(ESU))
```
""" josephson, KJ

@doc """
    magneticfluxquantum(U::UnitSystem) = planck(U)/2elementarycharge(U)/lorentz(U)

Magnetic flux quantum `Φ₀` is `1/josephson(U)` (Wb).
```Julia
julia> magneticfluxquantum(SI2019) # Wb
$(magneticfluxquantum(SI2019))

julia> magneticfluxquantum(Metric) # Wb
$(magneticfluxquantum(Metric))

julia> magneticfluxquantum(Conventional) # Wb
$(magneticfluxquantum(Conventional))

julia> magneticfluxquantum(EMU) # Mx
$(magneticfluxquantum(EMU))

julia> magneticfluxquantum(ESU) # statWb
$(magneticfluxquantum(ESU))
```
""" magneticfluxquantum, Φ₀

@doc """
    klitzing(U::UnitSystem) = planck(U)/elementarycharge(U)^2

Quantized Hall resistance `RK` (Ω).
```Julia
julia> klitzing(SI2019) # Ω
$(klitzing(SI2019))

julia> klitzing(Metric) # Ω
$(klitzing(Metric))

julia> klitzing(Conventional) # Ω
$(klitzing(Conventional))

julia> klitzing(CODATA) # Ω
$(klitzing(CODATA))

julia> klitzing(EMU) # abΩ
$(klitzing(EMU))

julia> klitzing(ESU) # statΩ
$(klitzing(ESU))
```
""" klitzing, RK

@doc """
    conductancequantum(U::UnitSystem) = 2elementarycharge(U)^2/planck(U) # 2/klitzing(U)

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

julia> conductancequantum(EMU) # abS
$(conductancequantum(EMU))

julia> conductancequantum(ESU) # statS
$(conductancequantum(ESU))
```
""" conductancequantum, G₀, G0

@doc """
    hartree(U::UnitSystem) = electronmass(U)*(lightspeed(U)*finestructure(U))^2 # mₑ*(𝘤/αinv)^2

Hartree electric potential energy `Eₕ` of the hydrogen atom at ground state is `2R∞*𝘩*𝘤` (J).
```Julia
julia> hartree(SI2019)/elementarycharge(SI2019) # eV
$(hartree(SI2019)/elementarycharge(SI2019))

julia> hartree(Metric) # J
$(hartree(Metric))

julia> hartree(CGS) # erg
$(hartree(CGS))

julia> hartree(Metric)*avogadro(Metric)/1000 # kJ⋅mol⁻¹
$(hartree(Metric)*avogadro(Metric)/(𝟐*𝟑)^3)

julia> hartree(Metric)*avogadro(Metric)/kilocalorie(Metric) # kcal⋅mol⁻¹
$(hartree(Metric)*avogadro(Metric)/kilocalorie(Metric))

julia> 2rydberg(Metric)/100 # Eₕ/𝘩/𝘤/100 cm⁻¹
$(hartree(Metric)/planck(Metric)/lightspeed(Metric)/(𝟐*𝟓)^2)

julia> hartree(Metric)/planck(Metric)/10^12 # THz
$(hartree(Metric)/planck(Metric))

julia> hartree(Metric)/boltzmann(Metric) # K
$(hartree(Metric)/boltzmann(Metric))
```
In a Gaussian unit system where `4π*ε₀ == 1` the Hartree energy is `𝘦^2/a₀`.
""" hartree, Eₕ, Eh

@doc """
    rydberg(U::UnitSystem) = hartree(U)/2planck(U)/lightspeed(U) # Eₕ/2𝘩/𝘤

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

julia> hartree(SI2019)/𝟐/charge(SI2019) # eV
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
    bohr(U) = planckreduced(U)/electronmass(U)/lightspeed(U)/finestructure(U)

Bohr radius of the hydrogen atom in its ground state `a₀` (m).
```Julia
julia> bohr(Metric) # m
$(bohr(Metric))
```
""" bohr, a₀, a0
#julia> bohr(Metric)/length(PlanckGauss) # ℓP
#$(bohr(Metric)/length(PlanckGauss))

@doc """
    bohrreduced(U::UnitSystem) = bohr(U)*(1+1/protonelectron(U))

Reduced Bohr radius including the effect of reduced mass in hydrogen atom (m).
```Julia
julia> bohrreduced(Metric) # m
$(bohrreduced(Metric))

julia> bohrreduced(Metric)/bohr(Metric) # a₀
$(bohrreduced(Metric)/bohr(Metric))
```
""" bohrreduced

@doc """
    electronradius(U) = finestructure(U)*planckreduced(U)/electronmass(U)/lightspeed(U)

Classical electron radius or Lorentz radius or Thomson scattering length (m).
```Julia
julia> electronradius(Metric) # m
$(electronradius(Metric))

julia> electronradius(CODATA) # m
$(electronradius(CODATA))

julia> electronradius(Conventional) # m
$(electronradius(Conventional))
```
""" electronradius, rₑ, re

@doc """
    magneton(U::UnitSystem) = elementarycharge(U)*planckreduced(U)*lorentz(U)/2electronmass(U)

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

julia> magneton(EMU2019) # erg⋅G⁻¹
$(magneton(EMU2019))

julia> magneton(ESU2019) # statA⋅cm²
$(magneton(ESU2019))

julia> magneton(SI2019)/elementarycharge(SI2019) # eV⋅T⁻¹
$(magneton(SI2019)/elementarycharge(SI2019))

julia> magneton(Hartree) # 𝘤⋅ħ⋅mₑ⁻¹
$(magneton(Hartree))
```
""" magneton, μB

@doc """
    hyperfine(U::UnitSystem) = frequency($ΔνCs,U)

Unperturbed groundstate hyperfine transition frequency `ΔνCs` of caesium-133 atom (Hz).
```Julia
julia> hyperfine(Metric) # Hz
$(hyperfine(Metric))
```
""" hyperfine, ΔνCs

@doc """
    hubble(U::UnitSystem) = time(U,Hubble)

Hubble parameter.
```Julia
julia> hubble(SI2019)
$(hubble(SI2019))

julia> hubble(Hubble)
$(hubble(Hubble))
```
""" hubble

@doc """
    cosmological(U::UnitSystem) = 3darkenergydensity(U)*(hubble(U)/lightspeed(U))^2

Cosmological constant.
```Julia
julia> cosmological(SI2019)
$(cosmological(SI2019))
```
""" cosmological

@doc """
    standardgravity(U::UnitSystem) = acceleration($g₀,U)

Standard gravity `acceleration` `g₀` at geodetic reference latitude (m⋅s⁻² or ft⋅s⁻²).
```Julia
julia> standardgravity(Metric) # m⋅s⁻²
$(standardgravity(Metric))

julia> standardgravity(English) # ft⋅s⁻²
$(standardgravity(English))

julia> standardgravity(Survey) # ftUS⋅s⁻²
$(standardgravity(Survey))
```
""" standardgravity, g₀, g0, lbm

@doc """
    pressure(U::UnitSystem) = pressure($atm,U)

Standard `pressure` reference level of one atmosphere `atm` (Pa or lb⋅ft⁻²).
```Julia
julia> standardpressure(Metric) # Pa
$(standardpressure(Metric))

julia> standardpressure(English) # lbm⋅ft⁻¹⋅s⁻²
$(standardpressure(English))

julia> standardpressure(Survey) # lbm⋅ftUS⁻¹⋅s⁻²
$(standardpressure(Survey))
```
""" standardpressure, atm

@doc """
    temperature(U::UnitSystem) = temperature($atm,U)

Standard `temperature` reference level at sea level (K or °R).
```Julia
julia> standardtemperature(Metric) # K
$(standardtemperature(Metric))

julia> standardtemperature(SI2019) # K
$(standardtemperature(SI2019))

julia> standardtemperature(English) # °R
$(standardtemperature(English))

julia> standardtemperature(English2019) # °R
$(standardtemperature(English2019))
```
""" standardtemperature, Tₛ

@doc """
    solarmass(U::UnitSystem) = mass($(GM☉/G),U)

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
```
""" solarmass, mₛ

@doc """
    earthmass(U::UnitSystem) = mass($(GME/G),U)

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
```
""" earthmass

@doc """
    jupitermass(U::UnitSystem) = mass($(GMJ/G),U)

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
```
""" jupitermass

@doc """
    lunarmass(U::UnitSystem) = earthmass(U)/μE☾

Lunar `mass` estimated from `μE☾` Earth-Moon mass ratio (kg or slug).
```Julia
julia> lunarmass(Metric) # kg
$(lunarmass(Metric))

julia> jupitermass(British) # slug
$(lunarmass(British))

julia> lunarmass(English) # lb
$(lunarmass(English))

julia> lunarmass(IAU) # M☉
$(lunarmass(IAU))

julia> lunarmass(IAUE) # ME
$(lunarmass(IAUE))

julia> lunarmass(IAUJ) # MJ
$(lunarmass(IAUJ))
```
""" lunarmass

@doc """
    astronomicalunit(U::UnitSystem) = length($au,U)

Standard astronomical unit from the International Astronomical Union (m or ft).
```Julia
julia> astronomicalunit(Metric) # m
$(astronomicalunit(Metric))

julia> astronomicalunit(English) # ft
$(astronomicalunit(English))

julia> astronomicalunit(Survey) # ftUS
$(astronomicalunit(Survey))
```
""" astronomicalunit, au

@doc """
    lunardistance(U::UnitSystem) = length($LD,U)

Standard distance between the Earth and the Moon (m or ft).
```Julia
julia> lunardistance(Metric) # m
$(lunardistance(Metric))

julia> lunardistance(English) # ft
$(lunardistance(English))

julia> lunardistance(Survey) # ftUS
$(lunardistance(Survey))
```
""" lunardistance, LD

@doc """
    mile(U::UnitSystem) = length($(Constant(5280)*ft),U)

Statute mile (m or ft).
```Julia
julia> mile(Metric) # m
$(mile(Metric))

julia> mile(English) # ft
$(mile(English))

julia> mile(Survey) # ftUS
$(mile(Survey))
```
""" mile

@doc """
    clarkemile(U::UnitSystem) = length($nm,U)

Historic nautical mile as defined by the Clarke (m or ft).
```Julia
julia> clarkemile(Metric) # m
$(clarkemile(Metric))

julia> clarkemile(English) # ft
$(clarkemile(English))

julia> clarkemile(Survey) # ftUS
$(clarkemile(Survey))
```
""" clarkemile

@doc """
    nauticalmile(U::UnitSystem) = length($nm,U)

Historic nautical mile as defined by the French (m or ft).
```Julia
julia> nauticalmile(Metric) # m
$(nauticalmile(Metric))

julia> nauticalmile(English) # ft
$(nauticalmile(English))

julia> nauticalmile(Survey) # ftUS
$(nauticalmile(Survey))
```
""" nauticalmile, nm

@doc """
    kilocalorie(U::UnitSystem) = energy(𝟐^5*𝟓^4*𝟑^2/𝟒𝟑,U,International)

Heat energy required to raise 1 kg of water by 1 Kelvin (`kcal`).
```Julia
julia> kilocalorie(International)
$(kilocalorie(International))

julia> kilocalorie(Metric)
$(kilocalorie(Metric))
```
""" kilocalorie, kcal

@doc """
    calorie(U::UnitSystem) = kilocalorie(U)/𝟐^3/𝟓^3

Heat energy required to raise 1 g of water by 1 Kelvin (`kcal`) in `International` scale.
```Julia
julia> calorie(International)
$(calorie(International))

julia> calorie(Metric)
$(calorie(Metric))
```
""" calorie, cal

@doc """
    meancalorie(U::UnitSystem) = energy(𝟐^2*𝟓*𝟑^2/𝟒𝟑,U,InternationalMean)

Heat energy required to raise 1 g of water by 1 Kelvin (`kcal`) in `InternationalMean`.
```Julia
julia> meancalorie(InternationalMean)
$(meancalorie(InternationalMean))

julia> meancalorie(Metric)
$(meancalorie(Metric))
```
""" meancalorie

@doc """
    thermalunit(U::UnitSystem) = kilocalorie(U)*𝟑^2/𝟓/lb

Heat energy required to raise 1 lb of water by 1 Rankine (`BTU`) in `International` scale.
```Julia
julia> thermalunit(British)
$(thermalunit(British))

julia> thermalunit(International)
$(thermalunit(International))

julia> thermalunit(Metric)
$(thermalunit(Metric))
```
""" thermalunit, BTU, BTUJ, BTUftlb

@doc """
    tonsrefrigeration(U::UnitSystem) = frequency(𝟐*𝟓/𝟑,U,Metric)*thermalunit(U)

Unit of `power` derived from melting of 1 short ton of ice in 24 hours.
```Julia
julia> tonsrefrigeration(British)
$(tonsrefrigeration(British))

julia> tonsrefrigeration(Metric)
$(tonsrefrigeration(Metric))
```
""" tonsrefrigeration

@doc """
    boilerhorsepower(U::UnitSystem) = frequency(1339/𝟐^4/𝟑^2,U,Metric)*thermalunit(U)

Unit of `power` derived from evaporating 34.5 lb of boiling water in 1 hour.
```Julia
julia> boilerhorsepower(British)
$(boilerhorsepower(British))

julia> boilerhorsepower(Metric)
$(boilerhorsepower(Metric))
```
""" boilerhorsepower

@doc """
    horsepower(U::UnitSystem) = power(𝟐*𝟓^2*𝟏𝟏,U,British)

Unit of `power` derived from raising 550 lb by 1 ft in 1  in 1 s.
```Julia
julia> horsepower(British)
$(horsepower(British))

julia> horsepower(Metric)
$(horsepower(Metric))

julia> horsepower(MetricEngineering)
$(horsepower(MetricEngineering))
```
""" horsepower, HP

@doc """
    horsepowerwatt(U::UnitSystem) = power(𝟐^4*𝟑^3/𝟓*τ,U,British)

Unit of `power` derived from Watt's exact original horse power estimate.
```Julia
julia> horsepowerwatt(British)
$(horsepowerwatt(British))

julia> horsepowerwatt(Metric)
$(horsepowerwatt(Metric))

julia> horsepowerwatt(MetricEngineering)
$(horsepowerwatt(MetricEngineering))
```
""" horsepowerwatt

@doc """
    horsepowermetric(U::UnitSystem) = power(𝟑*𝟓^2,U,GravitationalMetric)

Unit of `power` derived from raising 75 kp by 1 m in 1  in 1 s.
```Julia
julia> horsepowermetric(British)
$(horsepowermetric(British))

julia> horsepowermetric(Metric)
$(horsepowermetric(Metric))

julia> horsepowermetric(MetricEngineering)
$(horsepowermetric(MetricEngineering))
```
""" horsepowermetric

@doc """
    electricalhorsepower(U::UnitSystem) = power(746,U,Metric)

Unit of `power` for electrical motors in the United States.
```Julia
julia> electricalhorsepower(British)
$(electricalhorsepower(British))

julia> electricalhorsepower(Metric)
$(electricalhorsepower(Metric))

julia> electricalhorsepower(MetricEngineering)
$(electricalhorsepower(MetricEngineering))
```
""" electricalhorsepower

@doc """
    gallon(U::UnitSystem) = volume(𝟕*𝟏𝟏/𝟐^2,U,English)

Unit of `volume` derived from the US liquid `gallon` in cubic inches.
```Julia
julia> gallon(English)
$(gallon(English))

julia> gallon(Metric)
$(gallon(Metric))
```
""" gallon, gal

@doc """
    litre(U::UnitSystem) = volume(𝟏𝟎^-3,U,Metric)

Unit of `volume` derived from 1 cubic decimetre.
```Julia
julia> litre(Metric)
$(gallon(Metric))

julia> gallon(English)
$(gallon(English))
```
""" litre

@doc """
    inchmercury(U::UnitSystem) = pressure(inHg,U,Metric)

Unit of `pressure` exerted by 1 inch of mercury at standard atmospheric conditions.
```Julia
juila> inchmercury(Metric)
$(inchmercury(Metric))

julia> inchmercury(English)
$(inchmercury(English))
```
""" inchmercury, inHg

@doc """
    torr(U::UnitSystem) = pressure(atm/𝟐^3/𝟓/𝟏𝟗,U,Metric)

Unit of `pressure` exerted by 1 mm of mercury at standard atmospheric conditions.
```Julia
juila> torr(Metric)
$(torr(Metric))

julia> torr(English)
$(torr(English))
```
""" torr

@doc """
    second(U::UnitSystem) = time(𝟏,U,Metric)

Unit of `time` defined by `hyperfine` transition frequency of Cs-133 atom.
```Julia
julia> second(Metric)
$(second(Metric))

julia> second(IAU)
$(second(IAU))
```
""" second

@doc """
    minute(U::UnitSystem) = 𝟐^2*𝟑*𝟓*second(U)

Unit of `time` defined by 60 `second` intervals.
```Julia
julia> minute(Metric)
$(minute(Metric))

julia> minute(IAU)
$(minute(IAU))
```
""" minute

@doc """
    hour(U::UnitSystem) = 𝟐^2*𝟑*𝟓*minute(U)

Unit of `time` defined by 60 `minute` intervals.
```Julia
julia> hour(Metric)
$(hour(Metric))

julia> hour(IAU)
$(hour(IAU))
```
""" hour, HOUR

@doc """
    day(U::UnitSystem) = 𝟐^3*𝟑*hour(U)

Unit of `time` defined by 24 `hour` intervals.
```Julia
julia> day(Metric)
$(day(Metric))

julia> day(IAU)
$(day(IAU))
```
""" day, DAY

@doc """
    year(U::UnitSystem) = aⱼ*day(U)

Unit of `time` defined by Julian calendar year interval.
```Julia
julia> year(Metric)
$(year(Metric))

julia> year(IAU)
$(year(IAU))
```
""" year, aⱼ

@doc """
    gaussianyear(U::UnitSystem) = (τ/k)*day(U)

Unit of `time` defined by Gaussian gravitational constant.
```Julia
julia> gaussianyear(Metric)
$(gaussianyear(Metric))

julia> gaussianyear(IAU)
$(gaussianyear(IAU))
```
""" gaussianyear

@doc """
    siderealyear(U::UnitSystem) = τ/k/√(𝟏+earthmass(IAU)+lunarmass(IAU))*day(U)

Unit of `time` defined by Gaussian gravitational constant and the Earth system mass.
```Julia
julia> siderealyear(Metric)
$(siderealyear(Metric))

julia> siderealyear(IAU)
$(siderealyear(IAU))
```
""" siderealyear

@doc """
    lightyear(U::UnitSystem) = year(U)*lightspeed(U)

Unit of `length` defined by distance traveled by light in 1 `year` unit.
```Julia
julia> lightyear(Metric)
$(lightyear(Metric))

julia> lightyear(IAU)
$(lightyear(IAU))
```
""" lightyear, ly

@doc """
    parsec(U::UnitSystem) = astronomicalunit(U)*𝟐^2*𝟑^4*𝟓^3/τ

Unit of `length` defined at which 1 `astronomicalunit` subtends an angle of 1 arcsecond.
```Julia
julia> parsec(Metric)
$(parsec(Metric))

julia> parsec(IAU)
$(parsec(IAU))
```
""" parsec, pc

#=@pure kilogram(U::UnitSystem) = mass(Metric,U)
@pure slug(U::UnitSystem) = mass(English,U)

@pure meter(U::UnitSystem) = length(Metric,U)
@pure foot(U::UnitSystem) = length(English,U)=#

#rankine, kelvin, moles/molecules
#add gravitional units of weight??
