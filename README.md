# UnitSystems.jl

*Physical unit systems (Metric, English, Natural, etc...)*

[![DOI](https://zenodo.org/badge/317419353.svg)](https://zenodo.org/badge/latestdoi/317419353)
[![Build Status](https://travis-ci.org/chakravala/UnitSystems.jl.svg?branch=master)](https://travis-ci.org/chakravala/UnitSystems.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/r4gmftclfm30ik9n?svg=true)](https://ci.appveyor.com/project/chakravala/unitsystems-jl)
[![Coverage Status](https://coveralls.io/repos/chakravala/UnitSystems.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/chakravala/UnitSystems.jl?branch=master)
[![codecov.io](https://codecov.io/github/chakravala/UnitSystems.jl/coverage.svg?branch=master)](https://codecov.io/github/chakravala/UnitSystems.jl?branch=master)

> In fact there is nothing transcendental about dimensions; the ultimate principle is precisely expressible (in Newton's terminology) as one of *similitude*, exact or approximate, to be tested by the rule that mere change in the magnitudes of the ordered scheme of units of measurement that is employed must not affect sensibly the forms of the equations that are the adequate expression of the underlying relations of the problem. (J.L., 1914)

Specifications for dimensional units are in the [UnitSystems.jl](https://github.com/chakravala/UnitSystems.jl) and [Similitude.jl](https://github.com/chakravala/Similitude.jl) and [MeasureSystems.jl](https://github.com/chakravala/MeasureSystems.jl) repositories.
The three packages are designed so that they can be interchanged with compatibility.
On its own `UnitSystems` is the fastest package, while `Similitude` (provides `Quantity` type) and `MeasureSystems` (introduces [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl) uncertainty) build additional features on top of `UnitSystems` base defintions.
Additionally, in the `UnitSystems` repository there is an equivalent [Wolfram language paclet](https://reference.wolfram.com/language/guide/Paclets) `Kernel` and also an unmaintained Rust `src` implementation.
Defaults are shared across the packages: `Metric`, `SI2019`, `CODATA`, `Conventional`, `International`, `InternationalMean`, `MetricEngineering`, `SI2019Engineering`, `GravitationalMetric`, `GravitationalSI2019`, `British`, `British2019`, `Survey`, `Survey2019`, `English`, `English2019`, `FPS`, `FPS2019`, `Gauss`, `LorentzHeaviside`, `Thomson`, `EMU`, `ESU`, `EMU2019`, `ESU2019`, `IAU`, `IAUE`, `IAUJ`, `Astronomical`, `Hubble`, `Cosmological`, `CosmologicalQuantum`, `Nautical`, `MPH`, `KKH`, `MTS`, `FFF`, `Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal`.

```Julia
julia> using UnitSystems # or Similitude or MeasureSystems
```

A `UnitSystem` is a consistent set of dimensional values selected to accomodate a particular use case or standardization.
It is possible to convert derived physical quantities from any `UnitSystem` specification into any other using accurate values.
Eleven fundamental constants `kB`, `ħ`, `𝘤`, `μ₀`, `mₑ`, `Mᵤ`, `Kcd`, `θ`, `λ`, `αL`, `g₀` are used to govern a specific unit system consistent scaling.
These are the constants `boltzmann`, `planckreduced`, `lightspeed`, `vacuumpermeability`, `electronmass`, `molarmass`, `luminousefficacy`, `angle`, `rationalization`, `lorentz`, and `gravity`.
Different choices of natural units or physical measurements result in a variety of unit systems for many purposes.

Main documentation is at https://geophysics.crucialflow.com/dev/units

Historically, older electromagnetic unit systems also relied on a `rationalization` constant `λ` and a `lorentz` force proportionality constant `αL`.
In most unit systems these extra constants have a value of `1` unless otherwise specified.

```Julia
    UnitSystem{kB, ħ, 𝘤, μ₀, mₑ, Mᵤ, (Kcd, θ, λ, αL, g₀, ...)}
```

Fundamental constants of physics are: `kB` Boltzmann's constant, `ħ` reduced Planck's constant, `𝘤` speed of light, `μ₀` vacuum permeability, `mₑ` electron rest mass, `Mᵤ` molar mass, `Kcd` luminous efficacy, `θ` angle measure, `λ` Gauss rationalization, `αL` Lorentz's constant, and `g₀` gravitational force reference.
Primarily the `Metric` SI unit system is used in addition to the historic `English` engineering unit system.
These constants induce derived values for `avogadro`, `boltzmann`, `universalgas`, `planck`, `planckreduced`, `lightspeed`, `planckmass`, `atomicmass`, `protonmass`, `electronmass`, `newton`, `einstein`, `vacuumpermeability`, `vacuumpermittivity`, `coulomb`, and
additional constants `molarmass`, `luminousefficacy`, `gravity`, `angle`, `turn`, `sphere`, `stefan`, `radiationdensity`, `ampere`, `lorentz`, `biotsavart`, `rationalization`, `vacuumimpedance`, `elementarycharge`, `magneton`, `conductancequantum`, `faraday`, `magneticfluxquantum`, `josephson`, `klitzing`, `hartree`, `rydberg`, `bohr`, `bohrreduced`.
Derived quantities are `second`, `minute`, `hour`, `day`, `year`, `gaussianyear`, `siderealyear`, `hyperfine`, `hubble`, `cosmological`, `solarmass`, `earthmass`, `jupitermass`, `lunarmass`, `astronomicalunit`, `lunardistance`, `mile`, `clarkemile`, `nauticalmile`, `parsec`, `lightyear`, `gallon`, `litre`, `standardgravity`, `standardtemperature`, `standardpressure`, `inchmercury`, `torr`, `kilocalorie`, `calorie`, `meancalorie`, `thermalunit`, `tonsrefrigeration`, `horsepower`, `horsepowerwatt`, `horsepowermetric`, `electricalhorsepower`, `boilerhorsepower`.

Physics constant documentation is at https://geophysics.crucialflow.com/dev/constants

Additional reference `UnitSystem` variants: `EMU`, `ESU`, `Gauss`, `LorentzHeaviside`, `SI2019`, `SI1976`, `CODATA`, `Conventional`, `International`, `InternationalMean`, `MetricEngineering`, `GravitationalMetric`, `Astronomical`, `Hubble`, `Cosmological`, `CosmologicalQuantum`, `IAU`, `IAUE`, `IAUJ`, `MTS`, `FPS`, `British`, `Survey`, `Nautical`, `MPH`, `KKH`, `FFF`; and natural atomic units based on gravitational `coupling` and `finestructure` constant (`Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal`).

Unit conversion documentation is at https://geophysics.crucialflow.com/dev/convert

**Derived Unit conversions:**

Mechanics: `angle`, `solidangle`, `time`, `length`, `area`, `volume`, `wavenumber`, `angularwavenumber`, `fuelefficiency`, `frequency`, `angularfrequency`, `frequencydrift`, `speed`, `acceleration`, `jerk`, `snap`, `crackle`, `pop`, `volumeflow`,
`inertia`, `mass`, `massflow`, `lineardensity`, `areadensity`, `density`, `specificweight`, `specificvolume`, `force`, `gforce`, `stiffness`, `pressure`, `compressibility`, `viscosity`, `diffusivity`, `rotationalinertia`, `impulse`, `momentum`, `angularmomentum`, `yank`, `energy`, `specificenergy`, `action`, `fluence`, `power`, `powerdensity`, `intensity`, `spectralflux`, `soundexposure`, `impedance`, `specificimpedance`, `admittance`, `compliance`, `inertance`;
Electromagnetics: `charge`, `chargedensity`, `linearchargedensity`, `exposure`, `mobility`, `current`, `currentdensity`, `resistance`, `conductance`, `resistivity`, `conductivity`, `capacitance`, `inductance`, `reluctance`, `permeance`, `permittivity`, `permeability`, `susceptibility`, `specificsusceptibility`, `demagnetizingfactor`, `vectorpotential`, `electricpotential`, `magneticpotential`, `electricfield`, `magneticfield`, `electricflux`, `magneticflux`, `electricfluxdensity`, `magneticfluxdensity`, `electricdipolemoment`, `magneticdipolemoment`, `electricpolarizability`, `magneticpolarizability`, `magneticmoment`, `magnetizability`, `magnetization`, `specificmagnetization`, `rigidity`, `polestrength`;
Thermodynamics: `temperature`, `entropy`, `specificentropy`, `volumeheatcapacity`, `thermalconductivity`, `thermalconductance`, `thermalresistance`, `thermalexpansion`, `lapserate`,
`molarmass`, `molality`, `mole`, `molarity`, `molarvolume`, `molarentropy`, `molarenergy`, `molarconductivity`, `molarsusceptibility`, `catalysis`, `specificity`,
`luminousflux`, `luminance`, `luminousenergy`, `luminousexposure`, `luminousefficacy`.

**Generalized dimensionless `Coupling`:**

```Julia
Coupling{αG,α,μₑᵤ,μₚᵤ,ΩΛ}
```
Specification of `Universe` with the dimensionless `Coupling` constants `coupling`, `finestructure`, `electronunit`, `protonunit`, `protonelectron`, and `darkenergydensity`. Alterations to these values can be facilitated and quantified using parametric polymorphism.
Due to the `Coupling` interoperability, the `MeasureSystems` package is made possible to support calculations with `Measurements` having error standard deviations.

Other similar packages include [Similitude.jl](https://github.com/chakravala/Similitude.jl), [MeasureSystems.jl](https://github.com/chakravala/MeasureSystems.jl), [PhysicalConstants.jl](https://github.com/JuliaPhysics/PhysicalConstants.jl), [MathPhysicalConstants.jl](https://github.com/LaGuer/MathPhysicalConstants.jl), [Unitful.jl](https://github.com/PainterQubits/Unitful.jl.git), [UnitfulUS.jl](https://github.com/PainterQubits/UnitfulUS.jl), [UnitfulAstro.jl](https://github.com/JuliaAstro/UnitfulAstro.jl), [UnitfulAtomic.jl](https://github.com/sostock/UnitfulAtomic.jl), [NaturallyUnitful.jl](https://github.com/MasonProtter/NaturallyUnitful.jl), and [UnitfulMoles.jl](https://github.com/rafaqz/UnitfulMoles.jl).
