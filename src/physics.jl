
#   This file is part of UnitSystems.jl. It is licensed under the MIT license
#   UnitSystems Copyright (C) 2020 Michael Reed

@doc """
    μₑᵤ, μₚᵤ, μₚₑ, αinv, αG

Physical measured dimensionless values with uncertainty are the electron to proton mass ratio `μₑᵤ`, proton to atomic mass ratio `μₚᵤ`, proton to electron mass ratio `μₚₑ`, inverted fine structure constant `αinv`, and the gravitaional coupling constant `αG`.

```Julia
julia> μₑᵤ
$μₑᵤ

julia> μₚᵤ
$μₚᵤ

julia> μₚₑ
$μₚₑ

julia> αinv
$αinv

julia> αG
$αG
```
""" μₑᵤ, μₚᵤ, μₚₑ, αinv, αG, meu, mpu, mpe, ainv, aG

@pure hyperfine(U::UnitSystem) = frequency(ΔνCs,U)
@doc """
    hyperfine(U::UnitSystem) = frequency($ΔνCs,U)

Unperturbed groundstate hyperfine transition frequency `ΔνCs` of caesium-133 atom (Hz).
```Julia
julia> hyperfine(Metric) # Hz
$(hyperfine(Metric))
```
""" hyperfine, ΔνCs

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

julia> luminousefficacy(English) # lm⋅s³⋅slug⋅ft⁻²
$(luminousefficacy(English))
```
""" luminousefficacy, Kcd

@doc """
    molarmass(U::UnitSystem) = avogadro(U)*electronmass(U)/μₑᵤ # 1/μₑᵤ = $(1/μₑᵤ-2e-13)

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

@pure avogadro(U::UnitSystem) = μₑᵤ*molarmass(U)/electronmass(U)
@doc """
    avogadro(x) = universal(x)/boltzmann(x) # Mᵤ/atomicmass(x), Mᵤ ≈ 0.001-3.5e-13

Avogadro `NA` is `molarmass(x)/atomicmass(x)` number of atoms in a 12 g sample of C₁₂.
```Julia
julia> avogadro(SI2019) # mol⁻¹
$(avogadro(SI2019))

julia> avogadro(Metric) # mol⁻¹
$(avogadro(Metric))

julia> avogadro(English) # slug-mol⁻¹
$(avogadro(English))
```
""" avogadro, NA

@doc """
    planckreduced(x) = planck(x)/2π

Reduced Planck constant `ħ` is a Planck per radian (J⋅s⋅rad⁻¹ or ft⋅lb⋅s⋅rad⁻¹).

```Julia
julia> planckreduced(SI2019) # J⋅s⋅rad⁻¹
$(planckreduced(SI2019))

julia> planckreduced(SI2019)*lightspeed(SI2019) # J⋅m⋅rad⁻¹
$(planckreduced(SI2019)*𝘤)

julia> planckreduced(CODATA) # J⋅s⋅rad⁻¹
$(planckreduced(CODATA))

julia> planckreduced(Conventional) # J⋅s⋅rad⁻¹
$(planckreduced(Conventional))

julia> planckreduced(SI2019)/electronmass(SI2019) # eV⋅s⋅rad⁻¹
$(planckreduced(SI2019)/charge(SI2019))

julia> planckreduced(SI2019)*lightspeed(SI2019)/charge(SI2019) # eV⋅m⋅rad⁻¹
$(planckreduced(SI2019)*lightspeed(SI2019)/charge(SI2019))

julia> planckreduced(English) # ft⋅lb⋅s⋅rad⁻¹
$(planckreduced(English))
```
""" planckreduced, ħ

@doc """
    planck(x) = 2π*planckreduced(x)

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

julia> planck(SI2019)/electronmass(SI2019) # eV⋅s⋅rad⁻¹
$(planck(SI2019)/𝘦)

julia> planck(SI2019)*lightspeed(SI2019)/charge(SI2019) # eV⋅m⋅rad⁻¹
$(planck(SI2019)*lightspeed(SI2019)/charge(SI2019))

julia> planck(English) # ft⋅lb⋅s
$(planck(English))
```
""" planck, 𝘩, hh

@doc """
    boltzmann(x) = universal(x)/avogadro(x)

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

julia> boltzmann(SI2019)/charge(SI2019) # eV⋅K⁻¹
$(boltzmann(SI2019)/charge(SI2019))

julia> boltzmann(SI2019)/planck(SI2019) # Hz⋅K⁻¹
$(boltzmann(SI2019)/planck(SI2019))

julia> boltzmann(CGS) # erg⋅K⁻¹
$(boltzmann(CGS))

julia> boltzmann(SI2019)/calᵢₜ # calᵢₜ⋅K⁻¹
$(boltzmann(SI2019)/calᵢₜ)

julia> boltzmann(SI2019)*rankine/calᵢₜ # calᵢₜ⋅°R⁻¹
$(boltzmann(SI2019)*rankine/calᵢₜ)

julia> boltzmann(English) # ft⋅lb⋅°R⁻¹
$(boltzmann(English))

julia> boltzmann(SI2019)/planck(SI2019)/lightspeed(SI2019) # m⁻¹⋅K⁻¹
$(boltzmann(SI2019)/planck(SI2019)/lightspeed(SI2019))

julia> avogadro(SI2019)*boltzmann(SI2019)/calᵢₜ # calᵢₜ⋅mol⁻¹⋅K⁻¹
$(avogadro(SI2019)*boltzmann(SI2019)/calᵢₜ)

julia> 10log10(boltzmann(SI2019)/1) # dB(W⋅K⁻¹⋅Hz⁻¹)
$(10log10(boltzmann(SI2019)))
```
""" boltzmann, kB

@doc """
    lightspeed(U::UnitSystem) = 1/sqrt(permeability(U)*permittivity(U)*lorentz(U)^2)

Speed of light in a vacuum `𝘤` for massless particles (m⋅s⁻¹ or ft⋅s⁻¹).

```Julia
julia> lightspeed(Metric) # m⋅s⁻¹
$(lightspeed(Metric))

julia> lightspeed(English) # ft⋅s⁻¹
$(lightspeed(English))
```
""" lightspeed, 𝘤, cc

@doc """
    permeability(U::UnitSystem) = 1/permittivity(U)/(lightspeed(U)*lorentz(U))^2

Magnetic permeability in a classical vacuum defined as `μ₀` in SI units (H⋅m⁻¹, kg⋅m²⋅C⁻²).

```Julia
julia> permeability(Metric) # H⋅m⁻¹
$(permeability(Metric))

julia> permeability(Conventional) # H⋅m⁻¹
$(permeability(Conventional))

julia> permeability(CODATA) # H⋅m⁻¹
$(permeability(CODATA))

julia> permeability(SI2019) # H⋅m⁻¹
$(permeability(SI2019))

julia> permeability(EMU) # abH⋅cm⁻¹
$(permeability(EMU))

julia> permeability(ESU) # statH⋅cm⁻¹
$(permeability(ESU))
```
""" permeability, μ₀, m0

@doc """
    lorentz(U::UnitSystem) = 4π*biotsavart(U)/permeability(U)/rationalization(U)

Electromagnetic proportionality constant `αL` for the Lorentz's law force (?).

```Julia
julia> lorentz(Metric)
$(lorentz(Metric))

julia> lorentz(Thomson)
$(lorentz(Thomson))

julia> lorentz(Gauss)
$(lorentz(Gauss))
```
""" lorentz, αL, aL

@doc """
    rationalization(U::UnitSystem) = 4π*biotsavart(U)/permeability(U)/lorentz(U)

Constant of magnetization and polarization density or `4π*coulomb(U)*permittivity(U)` (?).

```Julia
julia> rationalization(Metric)
$(rationalization(Metric))

julia> rationalization(Gauss)
$(rationalization(Gauss))
```
""" rationalization

@doc """
    electronmass(U::UnitSystem) = protonmass(U)/$μₚₑ # αinv^2*R∞*2𝘩/𝘤

Electron rest mass `mₑ` of subatomic particle with `-𝘦` elementary charge  (kg or slugs).
```Julia
julia> electronmass(Metric) # kg
$(electronmass(Metric))

julia> electronmass(Metric)/atomicmass(Metric) # Da
$μₑᵤ

julia> electronmass(Metric)*lightspeed(Metric)^2 # J
$(electronmass(Metric)*lightspeed(Metric)^2)

julia> electronmass(SI2019)*lightspeed(SI2019)^2/charge(SI2019) # eV⋅𝘤⁻²
$(electronmass(SI2019)*lightspeed(SI2019)^2/charge(SI2019))

julia> electronmass(English) # slugs
$(electronmass(English))
```
""" electronmass, mₑ, me

@pure atomicmass(U::UnitSystem) = electronmass(U)/μₑᵤ
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

julia> atomicmass(SI2019)*lightspeed(SI2019)^2/charge(SI2019) # eV⋅𝘤⁻²
$(atomicmass(SI2019)*lightspeed(SI2019)^2/charge(SI2019))

julia> atomicmass(English) # slugs
$(atomicmass(English))
```
""" atomicmass, mᵤ, mu

@pure protonmass(U::UnitSystem) =  μₚₑ*electronmass(U)
@doc """
    protonmass(U::UnitSystem) = $(μₚᵤ)atomicmass(U)

Proton mass `mₚ` of subatomic particle with `+𝘦` elementary charge  (kg or mass).
```Julia
julia> protonmass(Metric) # kg
$(protonmass(Metric))

julia> protonmass(SI2019)*lightspeed(SI2019)^2/charge(SI2019) # eV⋅𝘤⁻²
$(protonmass(SI2019)*lightspeed(SI2019)^2/charge(SI2019))

julia> protonmass(Metric)/atomicmass(Metric) # mᵤ
$(protonmass(Metric)/atomicmass(Metric))

julia> protonmass(Metric)/electronmass(Metric) # mₑ
$(protonmass(Metric)/electronmass(Metric))
```
""" protonmass, mₚ, mp

@doc """
    planckmass(U::UnitSystem) = sqrt(planckreduced(U)*lightspeed(U)/newton(U))

Planck mass factor `mP` from the gravitational coupling constant `αG` (kg or slugs).
```Julia
juila> planckmass(Metric)*lightspeed(Metric)^2/charge(Metric) # eV⋅𝘤⁻²
$(planckmass(Metric)*lightspeed(Metric)^2/charge(Metric))

juila> planckmass(Metric) # kg
$(planckmass(Metric))

juila> planckmass(Metric)/atomicmass(Metric) # mᵤ
$(planckmass(Metric)/atomicmass(Metric))

juila> planckmass(Metric)*lightspeed(Metric)^2/charge(Metric)/sqrt(8π) # eV⋅𝘤⁻²
$(planckmass(Metric)*lightspeed(Metric)^2/charge(Metric)/sqrt(8π))

juila> planckmass(Metric)/sqrt(8π) # kg
$(planckmass(Metric)/sqrt(8π))
```
""" planckmass, mP

@doc """
    newton(U::UnitSystem) = lightspeed(U)*planckreduced(U)/planckmass(U)^2

Universal gravitational constant `GG` of Newton's law (m³⋅kg⁻¹⋅s⁻² or ft³⋅slug⁻¹⋅s⁻²).
```Julia
juila> newton(Metric) # m³⋅kg⁻¹⋅s⁻²
$(newton(Metric))

julia> newton(English) # ft³⋅slug⁻¹⋅s⁻²
$(newton(English))
```
""" newton, GG

@pure einstein(U::UnitSystem) = 8π*newton(U)/lightspeed(U)^4
@doc """
    einstein(U::UnitSystem) = 8π*newton(U)/lightspeed(U)^4

Einstein's gravitational constant from the Einstein field equations (s⋅²⋅m⁻¹⋅kg⁻¹).
```Julia
julia> einstein(Metric) # s⋅²⋅m⁻¹⋅kg⁻¹
$(einstein(Metric))
```
""" einstein, κ

@pure universal(U::UnitSystem) = boltzmann(U)*avogadro(U)
@doc """
    universal(x) = boltzmann(x)*avogadro(x)

Universal gas constant `Rᵤ` is factored into specific `gasconstant(x)*molarmass(x)` values.
```Julia
pressure*molarmass == density*universal*temperature
```
It satisfies the ideal gas law.

```Julia
julia> universal(SI2019) # J⋅K⁻¹⋅mol⁻¹
$(universal(SI2019))

julia> universal(English)*lbm/2116.2 # atm⋅ft³⋅R⁻¹⋅lb-mol⁻¹
$(universal(English)*lbm/2116.2)

julia> universal(Metric)/cal # cal⋅K⁻¹⋅mol⁻¹
$(universal(Metric)/cal)

julia> universal(Metric)/pressure(Earth1959) # atm⋅m³⋅K⁻¹⋅mol⁻¹
$(universal(Metric)/atm)

julia> universal(CGS) # erg⋅K⁻¹⋅mol⁻¹
$(universal(CGS))

julia> universal(English) # ft⋅lb⋅°R⁻¹⋅slug-mol⁻¹
$(universal(English))
```
The 1976 United States Standard Atmosphere used R* = 8.31432 exactly.
""" universal, Rᵤ, Ru

@pure stefan(U::UnitSystem) = 2π^5*boltzmann(U)^4/(15planck(U)^3*lightspeed(U)^2)
@doc """
    stefan(U::UnitSystem) = 2π^5*boltzmann(U)^4/(15planck(U)^3*lightspeed(U)^2)

Stefan-Boltzmann proportionality `σ` of black body radiation (W⋅m⁻²⋅K⁻⁴ or ?⋅ft⁻²⋅°R⁻⁴).

```Julia
julia> stefan(Metric) # W⋅m⁻²⋅K⁻⁴
$(stefan(Metric))

julia> stefan(CGS) # erg⋅cm⁻²⋅s⁻¹⋅K⁻⁴
$(stefan(CGS))

julia> stefan(Metric)*24*60^2/(cal*100^2) # cal⋅cm⁻²⋅day⁻¹⋅K⁻⁴
$(stefan(Metric)*24*0.6^2/cal)

julia> stefan(English) # lb⋅s⁻¹⋅ft⁻³⋅°R⁻⁴
$(stefan(English))
```
""" stefan, σ, SB

"""
    radiationdensity(U::UnitSystem) = 4stefan(U)/lightspeed(U)

Raditation density constant of black body radiation (J⋅m⁻³⋅K⁻⁴ or lb⋅ft⁻²⋅°R⁻⁴).

```Julia
julia> radiationdensity(Metric) # J⋅m⁻³⋅K⁻⁴
$(radiationdensity(Metric))

julia> radiationdensity(CGS) # erg⋅cm⁻³⋅K⁻⁴
$(radiationdensity(CGS))

julia> radiationdensity(English) # lb⋅ft⁻²⋅°R⁻⁴
$(radiationdensity(English))
```
"""
@pure radiationdensity(U::UnitSystem) = 4stefan(U)/lightspeed(U)

@pure permittivity(U::UnitSystem) = inv(permeability(U)*(lightspeed(U)*lorentz(U))^2)
@doc """
    permittivity(U::UnitSystem) = 1/permeability(U)/(lightspeed(U)*lorentz(U))^2

Dielectric permittivity constant `ε₀` of a classical vacuum (C²⋅N⁻¹⋅m⁻²).

```Julia
julia> permittivity(Metric) # F⋅m⁻¹
$(permittivity(Metric))

julia> permittivity(Conventional) # F⋅m⁻¹
$(permittivity(Conventional))

julia> permittivity(CODATA) # F⋅m⁻¹
$(permittivity(CODATA))

julia> permittivity(SI2019) # F⋅m⁻¹
$(permittivity(SI2019))

julia> permittivity(EMU) # abF⋅cm⁻¹
$(permittivity(EMU))

julia> permittivity(ESU) # statF⋅cm⁻¹
$(permittivity(ESU))

julia> permittivity(SI2019)/charge(SI2019) # 𝘦²⋅eV⁻¹⋅m⁻¹
$(permittivity(SI2019)/charge(SI2019))
```
""" permittivity, ε₀, ϵ₀, e0

@pure coulomb(U::UnitSystem) = rationalization(U)/4π/permittivity(U)
@doc """
    coulomb(U::UnitSystem) = rationalization(U)/4π/permittivity(U)

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

@pure biotsavart(U::UnitSystem) = permeability(U)*lorentz(U)*(rationalization(U)/4π)
@doc """
    biotsavart(U::UnitSystem) = permeability(U)*lorentz(U)*rationalization(U)/4π

Matnetostatic proportionality constant `αB` for the Biot-Savart's law (?).

```Julia
julia> biotsavart(Metric)
$(biotsavart(Metric))

julia> biotsavart(CODATA)
$(biotsavart(CODATA))

julia> biotsavart(SI2019)
$(biotsavart(SI2019))

julia> biotsavart(Conventional)
$(biotsavart(Conventional))

julia> biotsavart(EMU)
$(biotsavart(EMU))

julia> biotsavart(ESU)
$(biotsavart(ESU))

julia> biotsavart(Gauss)
$(biotsavart(Gauss))

julia> biotsavart(HLU)
$(biotsavart(HLU))
```
""" biotsavart, αB, aB

@pure ampere(U::UnitSystem) = lorentz(U)*biotsavart(U)
@doc """
    ampere(U::UnitSystem) = lorentz(U)*biotsavart(U) # coulomb(U)/lightspeed(U)^2

Magnetic proportionality constant `kₘ` for the Ampere's law force (N·s²⋅C⁻²).

```Julia
julia> ampere(Metric) # N·s²⋅C⁻²
$(ampere(Metric))

julia> ampere(CODATA) # N·s²⋅C⁻²
$(ampere(CODATA))

julia> ampere(SI2019) # N·s²⋅C⁻²
$(ampere(SI2019))

julia> ampere(Conventional) # N·s²⋅C⁻²
$(ampere(Conventional))

julia> ampere(EMU) # dyn·s²⋅abC⁻²
$(ampere(EMU))

julia> ampere(ESU) # dyn·s²⋅statC⁻²
$(ampere(ESU))

julia> ampere(HLU) # dyn·s²⋅hlC⁻²
$(ampere(HLU))
```
""" ampere, kₘ, km

@doc """
    impedance(U::UnitSystem) = permeability(U)*lightspeed(U)*rationalization(U)*lorentz(U)^2

Vacuum impedance of free space `Z₀` is magnitude ratio of electric to magnetic field (Ω).
```Julia
julia> impedance(Metric) # Ω
$(impedance(Metric))

julia> impedance(Conventional) # Ω
$(impedance(Conventional))

julia> impedance(CODATA) # Ω
$(impedance(CODATA))

julia> impedance(SI2019) # Ω
$(impedance(SI2019))

julia> 120π # 3e8*μ₀ # Ω
$(120π)

julia> impedance(EMU) # abΩ
$(impedance(EMU))

julia> impedance(ESU) # statΩ
$(impedance(ESU))

julia> impedance(HLU) # hlΩ
$(impedance(HLU))
```
""" impedance, Z₀, Z0

@doc """
    charge(U::UnitSystem) = sqrt(2𝘩/$(αinv)impedance(U)) # faraday(U)/avogadro(U)

Quantized elementary charge `𝘦` of a proton or electron `2/(klitzing(U)*josephson(U))` (C).
```Julia
julia> charge(SI2019) # C
$(charge(SI2019))

julia> charge(Metric) # C
$(charge(Metric))

julia> charge(CODATA) # C
$(charge(CODATA))

julia> charge(Conventional) # C
$(charge(Conventional))

julia> charge(EMU) # abC
$(charge(EMU))

julia> charge(ESU) # statC
$(charge(ESU))

julia> charge(Planck) # sqrt(4π/αinv)
$(charge(Planck))
```
""" charge, 𝘦, ee

@pure faraday(U::UnitSystem) = charge(U)*avogadro(U)
@doc """
    faraday(U::UnitSystem) = charge(U)*avogadro(U)

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

julia> faraday(Metric)/kcal # kcal⋅(V-g-e)⁻¹
$(faraday(Metric)/kcal)

julia> faraday(Metric)/3600 # A⋅h⋅mol⁻¹
$(faraday(Metric)/3600)
```
""" faraday, 𝔉, FF

@pure josephson(U::UnitSystem) = 2charge(U)*lorentz(U)/planck(U)
@doc """
    josephson(U::UnitSystem) = 2charge(U)*lorentz(U)/planck(U) # 1/magneticflux(U)

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

@pure magneticflux(U::UnitSystem) = inv(josephson(U))
@doc """
    magneticflux(U::UnitSystem) = planck(U)/2charge(U)/lorentz(U)

Magnetic flux quantum `Φ₀` is `1/josephson(U)` (Wb).
```Julia
julia> magneticflux(SI2019) # Wb
$(magneticflux(SI2019))

julia> magneticflux(Metric) # Wb
$(magneticflux(Metric))

julia> magneticflux(Conventional) # Wb
$(magneticflux(Conventional))

julia> magneticflux(EMU) # Mx
$(magneticflux(EMU))

julia> magneticflux(ESU) # statWb
$(magneticflux(ESU))
```
""" magneticflux, Φ₀

@pure klitzing(U::UnitSystem) = planck(U)/charge(U)^2
@doc """
    klitzing(U::UnitSystem) = planck(U)/charge(U)^2

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

@pure conductance(U::UnitSystem) = 2charge(U)^2/planck(U)
@doc """
    conductance(U::UnitSystem) = 2charge(U)^2/planck(U) # 2/klitzing(U)

Conductance quantum `G₀` is a quantized unit of electrical conductance (S).
```Julia
julia> conductance(SI2019) # S
$(conductance(SI2019))

julia> conductance(Metric) # S
$(conductance(Metric))

julia> conductance(Conventional) # S
$(conductance(Conventional))

julia> conductance(CODATA) # S
$(conductance(CODATA))

julia> conductance(EMU) # abS
$(conductance(EMU))

julia> conductance(ESU) # statS
$(conductance(ESU))
```
""" conductance, G₀, G0

@pure hartree(U::UnitSystem) = electronmass(U)*(lightspeed(U)/αinv)^2
@doc """
    hartree(U::UnitSystem) = electronmass(U)*(lightspeed(U)/$αinv)^2 # mₑ*(𝘤/αinv)^2

Hartree electric potential energy `Eₕ` of the hydrogen atom at ground state is `2R∞*𝘩*𝘤` (J).
```Julia
julia> hartree(SI2019)/charge(SI2019) # eV
$(hartree(SI2019)/charge(SI2019))

julia> hartree(Metric) # J
$(hartree(Metric))

julia> hartree(CGS) # erg
$(hartree(CGS))

julia> hartree(Metric)*avogadro(Metric)/1000 # kJ⋅mol⁻¹
$(hartree(Metric)*avogadro(Metric)/1000)

julia> hartree(Metric)*avogadro(Metric)/kcal # kcal⋅mol⁻¹
$(hartree(Metric)*avogadro(Metric)/kcal)

julia> 2rydberg(Metric)/100 # Eₕ/𝘩/𝘤/100 cm⁻¹
$(hartree(Metric)/planck(Metric)/lightspeed(Metric)/100)

julia> hartree(Metric)/planck(Metric)/10^12 # THz
$(hartree(Metric)/planck(Metric))

julia> hartree(Metric)/boltzmann(Metric) # K
$(hartree(Metric)/boltzmann(Metric))
```
In a Gaussian unit system where `4π*ε₀ == 1` the Hartree energy is `𝘦^2/a₀`.
""" hartree, Eₕ, Eh

@pure rydberg(U::UnitSystem) = hartree(U)/2planck(U)/lightspeed(U)
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
$(hartree(Metric)/2)

julia> hartree(SI2019)/2charge(SI2019) # eV
$(hartree(SI2019)/2charge(SI2019))
```
Rydberg photon frequency `𝘤*R∞` or `Eₕ/2𝘩` (Hz).
```Julia
julia> lightspeed(Metric)*rydberg(Metric) # Hz
$(lightspeed(Metric)*rydberg(Metric))
```
Rydberg wavelength `1/R∞` (m).
```Julia
julia> 1/rydberg(Metric) # m
$(1/rydberg(Metric))

julia> 1/rydberg(Metric)/2π # m⋅rad⁻¹
$(1/rydberg(Metric)/2π)
```
Precision measurements of the Rydberg constants are within a relative standard uncertainty of under 2 parts in 10¹², and is chosen to constrain values of other physical constants.
""" rydberg, R∞, RH, Ry

@pure bohr(U::UnitSystem) = αinv*planckreduced(U)/electronmass(U)/lightspeed(U)
@doc """
    bohr(U) = $αinv*planckreduced(U)/electronmass(U)/lightspeed(U)

Bohr radius of the hydrogen atom in its ground state `a₀` (m).
```Julia
julia> bohr(Metric) # m
$(bohr(Metric))

julia> 12bohr(English) # in
$(12bohr(English))

julia> bohr(Metric)/length(PlanckGauss) # ℓP
$(bohr(Metric)/length(PlanckGauss))
```
""" bohr, a₀, a0

"""
    bohrreduced(U::UnitSystem) = bohr(U)*(1+1/$μₚₑ)

Reduced Bohr radius including the effect of reduced mass in hydrogen atom (m).
```Julia
julia> bohrreduced(Metric) # m
$(bohrreduced(Metric))

julia> bohrreduced(Metric) # a₀
$(bohrreduced(Metric)/bohr(Metric))
```
"""
@pure bohrreduced(U::UnitSystem) = bohr(U)*(1+1/μₚₑ)

@pure electronradius(U::UnitSystem) = planckreduced(U)/electronmass(U)/lightspeed(U)/αinv
@doc """
    electronradius(U) = planckreduced(U)/electronmass(U)/lightspeed(U)/$αinv

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

@pure magneton(U::UnitSystem) = charge(U)*planckreduced(U)*lorentz(U)/2electronmass(U)
"""
    magneton(U::UnitSystem) = charge(U)*planckreduced(U)*lorentz(U)/2electronmass(U)

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

julia> magneton(SI2019)/charge(SI2019) # eV⋅T⁻¹
$(magneton(SI2019)/charge(SI2019))

julia> magneton(Hartree) # 𝘤⋅ħ⋅mₑ⁻¹
$(magneton(Hartree))
```
""" magneton, μB
