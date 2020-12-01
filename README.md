# UnitSystems.jl

*Physical unit system constants (Metric, English, Natural, etc...)*

[![Build Status](https://travis-ci.org/chakravala/UnitSystems.jl.svg?branch=master)](https://travis-ci.org/chakravala/UnitSystems.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/r4gmftclfm30ik9n?svg=true)](https://ci.appveyor.com/project/chakravala/unitsystems-jl)
[![Coverage Status](https://coveralls.io/repos/chakravala/UnitSystems.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/chakravala/UnitSystems.jl?branch=master)
[![codecov.io](https://codecov.io/github/chakravala/UnitSystems.jl/coverage.svg?branch=master)](https://codecov.io/github/chakravala/UnitSystems.jl?branch=master)

A `UnitSystem` is a consistent set of dimensional values selected to accomodate a particular use case or standardization.
In total, five fundamental constants `kB,ħ,𝘤,μ₀,mₑ` are used to specify a specific unit system. These are the constants of `boltzmann`, `planckreduced`, `lightspeed`, `permeability`, and `electronmass`.
Different choices of natural units or physical measurements result in a variety of unit systems optimized for many purposes.

```Julia
    UnitSystem{kB,ħ,𝘤,μ₀,mₑ}
```

Standardized for engineering based on fundamental constants: `kB` Boltzmann's constant, `ħ` reduced Planck's constant, `𝘤` speed of light, `μ₀` vacuum permeability, and `mₑ` electron rest mass.
Primarily the `Metric` SI unit system is used in addition to the historic `English` engineering unit system.
These constants induce derived values for `avogadro`, `boltzmann`, `universal`, `planck`, `planckreduced`, `lightspeed`, `planckmass`, `atomicmass`, `protonmass`, `electronmass`, `newton`, `einstein`, `permeability`, `permittivity`, `coulomb`, and
additional constants `stefan`, `radiationintensity`, `impedance`, `charge`, `magneton`, `conductance`, `faraday`, `magneticflux`, `josephson`, `klitzing`, `hartree`, `rydberg`, `bohr`, `bohrreduced`, and `molarmass`.

Main documentation is at https://geophysics.crucialflow.com/dev/units

Additional reference `UnitSystem` variants `CGS`, `CGS2019`, `SI2019`, `CODATA`, `Conventional`, `IAU`; along with several natural atomic units based on the fine structure constant `1/αinv` and the gravitational coupling constant `αG` (`Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal`).

Other similar packages include [PhysicalConstants.jl](https://github.com/JuliaPhysics/PhysicalConstants.jl), [MathPhysicalConstants.jl](https://github.com/LaGuer/MathPhysicalConstants.jl), [Unitful.jl](https://github.com/PainterQubits/Unitful.jl.git), [UnitfulUS](https://github.com/PainterQubits/UnitfulUS.jl), [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl), [UnitfulAtomic](https://github.com/sostock/UnitfulAtomic.jl), [NaturallyUnitful](https://github.com/MasonProtter/NaturallyUnitful.jl), and [UnitfulMoles](https://github.com/rafaqz/UnitfulMoles.jl).
