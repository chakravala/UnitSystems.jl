# UnitSystems.jl

*Physical unit systems (Metric, English, Natural, etc...)*

[![DOI](https://zenodo.org/badge/317419353.svg)](https://zenodo.org/badge/latestdoi/317419353)
[![Build Status](https://travis-ci.org/chakravala/UnitSystems.jl.svg?branch=master)](https://travis-ci.org/chakravala/UnitSystems.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/r4gmftclfm30ik9n?svg=true)](https://ci.appveyor.com/project/chakravala/unitsystems-jl)
[![Coverage Status](https://coveralls.io/repos/chakravala/UnitSystems.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/chakravala/UnitSystems.jl?branch=master)
[![codecov.io](https://codecov.io/github/chakravala/UnitSystems.jl/coverage.svg?branch=master)](https://codecov.io/github/chakravala/UnitSystems.jl?branch=master)

Specifications for dimensional units are in the [UnitSystems.jl](https://github.com/chakravala/UnitSystems.jl) and [UnitfulSystems.jl](https://github.com/chakravala/UnitfulSystems.jl) repositories.
The two packages are designed so that they can be interchanged if compatibility with [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) is desired or not.
However, the `UnitfulSystems` package has fewer `UnitSystem` specifications available than the `UnitSystems` package due to limitations in combination with the `Unitful` package.
Specifically, `Metric`, `SI2019`, `CODATA`, `Conventional`, `MTS`, `EMU2019`, `English`, and `EnglishUS` can have `Unitful` values; while `Gauss`, `LorentzHeaviside`, `Thomson`, `EMU`, `ESU`, `ESU2019`, `IAU`, `FFF`, `Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal` are plain valued.

```Julia
pkg> add UnitSystems # or UnitfulSystems

julia> using UnitSystems
```

A `UnitSystem` is a consistent set of dimensional values selected to accomodate a particular use case or standardization.
It is possible to convert derived physical quantities from any `UnitSystem` specification into any other using accurate values.
In total, five fundamental constants `kB,ħ,𝘤,μ₀,mₑ` are used to specify a specific unit system.
These are the constants of `boltzmann`, `planckreduced`, `lightspeed`, `permeability`, and `electronmass`.
Different choices of natural units or physical measurements result in a variety of unit systems for many purposes.

Main documentation is at https://geophysics.crucialflow.com/dev/units

Another important additional definition is the `molarmass` constant `Mᵤ`, which is automatically selected based on the choice of `boltzmann` constant (but can also be customized if necessary).
Historically, older electromagnetic unit systems also relied on a `rationalization` constant `λ` and a `lorentz` force proportionality constant `αL`.
In most unit systems these extra constants have a value of `1` unless otherwise specified.

```Julia
    UnitSystem{kB,ħ,𝘤,μ₀,mₑ,λ,αL}
```

Fundamental constants of physics are: `kB` Boltzmann's constant, `ħ` reduced Planck's constant, `𝘤` speed of light, `μ₀` vacuum permeability, `mₑ` electron rest mass, `λ` Gauss rationalization, and `αL` Lorentz's constant.
Primarily the `Metric` SI unit system is used in addition to the historic `English` engineering unit system.
These constants induce derived values for `avogadro`, `boltzmann`, `universal`, `planck`, `planckreduced`, `lightspeed`, `planckmass`, `atomicmass`, `protonmass`, `electronmass`, `newton`, `einstein`, `permeability`, `permittivity`, `coulomb`, and
additional constants `molarmass`, `hyperfine`, `luminousefficacy`, `stefan`, `radiationintensity`, `ampere`, `lorentz`, `biotsavart`, `rationalization`, `impedance`, `charge`, `magneton`, `conductance`, `faraday`, `magneticflux`, `josephson`, `klitzing`, `hartree`, `rydberg`, `bohr`, and `bohrreduced`.

Physics constant documentation is at https://geophysics.crucialflow.com/dev/constants

Additional reference `UnitSystem` variants: `EMU`, `ESU`, `Gauss`, `LorentzHeaviside`, `MTS`, `SI2019`, `CODATA`, `Conventional`, `IAU`, `EnglishUS`; and natural atomic units based on gravitational coupling `αG` and the fine structure `1/αinv` constant (`Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal`).

Unt conversion documentation is at https://geophysics.crucialflow.com/dev/convert

**Derived Unit conversions:**

Mechanics: `time`, `length`, `area`, `volume`, `wavenumber`, `fuelefficiency`, `frequency`, `frequencydrift`, `speed`, `acceleration`, `jerk`, `snap`, `volumeflow`,
`mass`, `massflow`, `lineardensity`, `areadensity`, `density`, `specificvolume`, `force`, `stiffness`, `pressure`, `compressibility`, `viscosity`, `diffusivity`, `rotationalinertia`, `momentum`, `angularmomentum`, `yank`, `energy`, `specificenergy`, `action`, `fluence`, `power`, `powerdensity`, `intensity`, `spectralflux`, `soundexposure`, `impedance`, `specificimpedance`, `admittance`, `compliance`, `inertance`;
Electromagnetics: `charge`, `chargedensity`, `linearchargedensity`, `exposure`, `mobility`, `current`, `currentdensity`, `resistance`, `conductance`, `resistivity`, `conductivity`, `capacitance`, `inductance`, `reluctance`, `permeance`, `permittivity`, `permeability`, `susceptibility`, `specificsusceptibility`, `demagnetizingfactor`, `vectorpotential`, `electricpotential`, `magneticpotential`, `electricfield`, `magneticfield`, `electricflux`, `magneticflux`, `electricfluxdensity`, `magneticfluxdensity`, `electricdipolemoment`, `magneticdipolemoment`, `electricpolarizability`, `magneticpolarizability`, `magneticmoment`, `magnetizability`, `magnetization`, `specificmagnetization`, `rigidity`, `polestrength`;
Thermodynamics: `temperature`, `entropy`, `specificentropy`, `volumeheatcapacity`, `thermalconductivity`, `thermalconductance`, `thermalresistance`, `thermalexpansion`, `lapserate`,
`molarmass`, `molality`, `mole`, `molarity`, `molarvolume`, `molarentropy`, `molarenergy`, `molarconductivity`, `molarsusceptibility`, `catalysis`, `specificity`,
`luminousflux`, `luminance`, `luminousenergy`, `luminousexposure`, `luminousefficacy`.

Other similar packages include [PhysicalConstants.jl](https://github.com/JuliaPhysics/PhysicalConstants.jl), [MathPhysicalConstants.jl](https://github.com/LaGuer/MathPhysicalConstants.jl), [Unitful.jl](https://github.com/PainterQubits/Unitful.jl.git), [UnitfulSystems.jl](https://github.com/chakravala/UnitfulSystems.jl), [UnitfulUS.jl](https://github.com/PainterQubits/UnitfulUS.jl), [UnitfulAstro.jl](https://github.com/JuliaAstro/UnitfulAstro.jl), [UnitfulAtomic.jl](https://github.com/sostock/UnitfulAtomic.jl), [NaturallyUnitful.jl](https://github.com/MasonProtter/NaturallyUnitful.jl), and [UnitfulMoles.jl](https://github.com/rafaqz/UnitfulMoles.jl).
