# UnitSystems.jl

*Physical unit systems (Metric, English, Natural, etc...)*

[![DOI](https://zenodo.org/badge/317419353.svg)](https://zenodo.org/badge/latestdoi/317419353)
[![Build Status](https://travis-ci.org/chakravala/UnitSystems.jl.svg?branch=master)](https://travis-ci.org/chakravala/UnitSystems.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/r4gmftclfm30ik9n?svg=true)](https://ci.appveyor.com/project/chakravala/unitsystems-jl)
[![Coverage Status](https://coveralls.io/repos/chakravala/UnitSystems.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/chakravala/UnitSystems.jl?branch=master)
[![codecov.io](https://codecov.io/github/chakravala/UnitSystems.jl/coverage.svg?branch=master)](https://codecov.io/github/chakravala/UnitSystems.jl?branch=master)

In aggregate, the `UnitSystem` data generated here constitutes a new universal standardization for dimensional analysis, which generalizes upon previous historical systems up to the 2019 redefinition and unifies them in a common `Universe`.
This enables a more precise and generalized standardization than the 2019 redefinition, which was comparatively limited in scope.
Specified default `UnitSystem` values are to be taken as a newly defined mutually-compatible recommended standard, verified to be consistent and coherent.
A `UnitSystem` can only be useful as a measuring standard if it can be scientifically reproduced, so the data here has been implemented in several important scientific programming languages (initially in the Julia language but also Wolfram language and Rust langauge) as well as presented abstractly in terms of dimensional formulas.

> In fact there is nothing transcendental about dimensions; the ultimate principle is precisely expressible (in Newton's terminology) as one of *similitude*, exact or approximate, to be tested by the rule that mere change in the magnitudes of the ordered scheme of units of measurement that is employed must not affect sensibly the forms of the equations that are the adequate expression of the underlying relations of the problem. (J.L., 1914)

Specifications for dimensional units are in the [UnitSystems.jl](https://github.com/chakravala/UnitSystems.jl) and [Similitude.jl](https://github.com/chakravala/Similitude.jl) and [MeasureSystems.jl](https://github.com/chakravala/MeasureSystems.jl) repositories.
The three packages are designed so that they can be interchanged with compatibility.
On its own `UnitSystems` is the fastest package, while `Similitude` (provides `Quantity` type) and `MeasureSystems` (introduces [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl) uncertainty) build additional features on top of `UnitSystems` base defintions.
Additionally, in the `UnitSystems` repository there is an equivalent [Wolfram language paclet](https://reference.wolfram.com/language/guide/Paclets) `Kernel` and also an unmaintained Rust `src` implementation.
Defaults are shared: `Metric`, `SI2019`, `CODATA`, `Conventional`, `International`, `InternationalMean`, `MetricEngineering`, `SI2019Engineering`, `GravitationalMetric`, `GravitationalSI2019`, `FPS`, `IPS`, `British`, `English`, `Survey`, `Gauss`, `LorentzHeaviside`, `EMU`, `ESU`, `IAU`, `IAUE`, `IAUJ`, `Hubble`, `Cosmological`, `CosmologicalQuantum`, `Meridian`, `Nautical`, `MPH`, `KKH`, `MTS`, `FFF`, `Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, `QCDoriginal`.

```Julia
julia> using UnitSystems # or Similitude or MeasureSystems
```

A `UnitSystem` is a consistent set of dimensional values selected to accomodate a particular use case or standardization.
It is possible to convert derived physical quantities from any `UnitSystem` specification into any other using accurate values.
Eleven fundamental constants `kB`, `ħ`, `𝘤`, `μ₀`, `mₑ`, `Mᵤ`, `Kcd`, `θ`, `λ`, `αL`, `g₀` are used to govern a specific unit system consistent scaling.
These are the constants `boltzmann`, `planckreduced`, `lightspeed`, `vacuumpermeability`, `electronmass`, `molarmass`, `luminousefficacy`, `angle`, `rationalization`, `lorentz`, and `gravity`.
Different choices of natural units or physical measurements result in a variety of unit systems for many purposes.

Main documentation is at https://geophysics.crucialflow.com/dev/unitsystems

Historically, older electromagnetic unit systems also relied on a `rationalization` constant `λ` and a `lorentz` force proportionality constant `αL`.
In most unit systems these extra constants have a value of `1` unless specified.

```Julia
    UnitSystem{kB, ħ, 𝘤, μ₀, mₑ, Mᵤ, (Kcd, θ, λ, αL, g₀, ...)}
```

Fundamental constants of physics are: `kB` Boltzmann's constant, `ħ` reduced Planck's constant, `𝘤` speed of light, `μ₀` vacuum permeability, `mₑ` electron rest mass, `Mᵤ` molar mass, `Kcd` luminous efficacy, `θ` angle measure, `λ` Gauss rationalization, `αL` Lorentz's constant, and `g₀` gravitational force reference.
Primarily the `Metric` SI unit system is used in addition to the historic `English` engineering unit system.
These constants induce derived values for `avogadro`, `boltzmann`, `molargas`, `planck`, `planckreduced`, `lightspeed`, `planckmass`, `atomicmass`, `protonmass`, `electronmass`, `newton`, `einstein`, `vacuumpermeability`, `vacuumpermittivity`, `electrostatic`, and
additional constants `molarmass`, `luminousefficacy`, `gravity`, `angle`, `turn`, `spat`, `stefan`, `radiationdensity`, `magnetostatic`, `lorentz`, `biotsavart`, `rationalization`, `vacuumimpedance`, `elementarycharge`, `magneton`, `conductancequantum`, `faraday`, `magneticfluxquantum`, `josephson`, `klitzing`, `hartree`, `rydberg`, `bohr`.

Physics constant documentation is at https://geophysics.crucialflow.com/dev/constants

Standardized unit/derived quantities are `hyperfine`, `loschmidt`, `wienwavelength`, `wienfrequency`, `mechanicalheat`, `solarmass`, `jupitermass`, `earthmass`, `lunarmass`, `earthradius`, `greatcircle`, `radarmile`, `hubble`, `cosmological`, `radian`, `steradian`, `degree`, `gradian`, `arcminute`, `arcsecond`, `second`, `minute`, `hour`, `day`, `year`, `gaussianyear`, `siderealyear`, `angstrom`, `inch`, `foot`, `surveyfoot`, `yard`, `meter`, `earthmeter`, `mile`, `statutemile`, `meridianmile`, `admiraltymile`, `nauticalmile`, `lunardistance`, `astronomicalunit`, `lightyear`, `parsec`, `barn`, `hectare`, `acre`, `surveyacre`, `liter`, `gallon`, `quart`, `pint`, `cup`, `fluidounce`, `teaspoon`, `tablespoon`, `grain`, `gram`, `earthgram`, `kilogram`, `tonne`, `ton`, `pound`, `ounce`, `slug`, `slinch`, `hyl`, `dyne`, `newton`, `poundal`, `poundforce`, `kilopond`, `psi`, `pascal`, `bar`, `barye`, `technicalatmosphere`, `atmosphere`, `inchmercury`, `torr`, `electronvolt`, `erg`, `joule`, `footpound`, `calorie`, `kilocalorie`, `meancalorie`, `earthcalorie`, `thermalunit`, `gasgallon`, `tontnt`, `watt`, `horsepower`, `horsepowerwatt`, `horsepowermetric`, `electricalhorsepower`, `tonsrefrigeration`, `boilerhorsepower`, `coulomb`, `earthcoulomb`, `ampere`, `volt`, `henry`, `ohm`, `siemens`, `farad`, `weber`, `tesla`, `abcoulomb`, `abampere`, `abvolt`, `abhenry`, `abohm`, `abmho`, `abfarad`, `maxwell`, `gauss`, `oersted`, `gilbert`, `statcoulomb`, `statampere`, `statvolt`, `stathenry`, `statohm`, `statmho`, `statfarad`, `statweber`, `stattesla`, `kelvin`, `rankine`, `sealevel`, `mole`, `earthmole`, `poundmole`, `slugmole`, `slinchmole`, `katal`, `amagat`, `lumen`, `candela`, `lux`, `phot`, `footcandle`, `nit`, `apostilb`, `stilb`, `lambert`, `footlambert`, `bril`, `neper`, `bel`, `decibel`, `hertz`, `rpm`, `kayser`, `diopter`, `bubnoff`, `gforce`, `galileo`, `eotvos`, `darcy`, `poise`, `reyn`, `stokes`, `rayl`, `mpge`, `langley`, `jansky`, `solarflux`, `curie`, `sievert`, `roentgen`, `rem`.

Standard physics units are at https://geophysics.crucialflow.com/dev/units

Additional reference `UnitSystem` variants: `EMU`, `ESU`, `Gauss`, `LorentzHeaviside`, `SI2019`, `SI1976`, `CODATA`, `Conventional`, `International`, `InternationalMean`, `MetricEngineering`, `GravitationalMetric`, `IAU`, `IAUE`, `IAUJ`, `FPS`, `IPS`, `British`, `Survey`, `Hubble`, `Cosmological`, `CosmologicalQuantum`, `Meridian`, `Nautical`, `MPH`, `KKH`, `MTS`, `FFF`; and natural atomic units based on gravitational `coupling` and `finestructure` constant (`Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal`).

Unit conversion documentation is at https://geophysics.crucialflow.com/dev/convert

**Derived Unit conversions:**

Mechanics: `angle`, `solidangle`, `time`, `length`, `area`, `volume`, `wavenumber`, `angularwavenumber`, `fuelefficiency`, `numberdensity`, `frequency`, `angularfrequency`, `frequencydrift`, `speed`, `acceleration`, `jerk`, `snap`, `crackle`, `pop`, `volumeflow`,
`inertia`, `mass`, `massflow`, `lineardensity`, `areadensity`, `density`, `specificweight`, `specificvolume`, `force`, `specificforce`, `gravityforce`, `pressure`, `compressibility`, `viscosity`, `diffusivity`, `rotationalinertia`, `impulse`, `momentum`, `angularmomentum`, `yank`, `energy`, `specificenergy`, `action`, `fluence`, `power`, `powerdensity`, `intensity`, `spectralflux`, `soundexposure`, `impedance`, `specificimpedance`, `admittance`, `compliance`, `inertance`;
Electromagnetics: `charge`, `chargedensity`, `linearchargedensity`, `exposure`, `mobility`, `current`, `currentdensity`, `resistance`, `conductance`, `resistivity`, `conductivity`, `capacitance`, `inductance`, `reluctance`, `permeance`, `permittivity`, `permeability`, `susceptibility`, `specificsusceptibility`, `demagnetizingfactor`, `vectorpotential`, `electricpotential`, `magneticpotential`, `electricfield`, `magneticfield`, `electricflux`, `magneticflux`, `electricfluxdensity`, `magneticfluxdensity`, `electricdipolemoment`, `magneticdipolemoment`, `electricpolarizability`, `magneticpolarizability`, `magneticmoment`, `specificmagnetization`, `polestrength`;
Thermodynamics: `temperature`, `entropy`, `specificentropy`, `volumeheatcapacity`, `thermalconductivity`, `thermalconductance`, `thermalresistivity`, `thermalresistance`, `thermalexpansion`, `lapserate`,
`molarmass`, `molality`, `mole`, `molarity`, `molarvolume`, `molarentropy`, `molarenergy`, `molarconductivity`, `molarsusceptibility`, `catalysis`, `specificity`,
`luminousflux`, `luminousintensity`, `luminance`, `illuminance`, `luminousenergy`, `luminousexposure`, `luminousefficacy`.

**Generalized dimensionless `Coupling`:**

```Julia
Coupling{αG,α,μₑᵤ,μₚᵤ,ΩΛ}
```
Specification of `Universe` with the dimensionless `Coupling` constants `coupling`, `finestructure`, `electronunit`, `protonunit`, `protonelectron`, and `darkenergydensity`. Alterations to these values can be facilitated and quantified using parametric polymorphism.
Due to the `Coupling` interoperability, the `MeasureSystems` package is made possible to support calculations with `Measurements` having error standard deviations.

Other similar packages include [Similitude.jl](https://github.com/chakravala/Similitude.jl), [MeasureSystems.jl](https://github.com/chakravala/MeasureSystems.jl), [PhysicalConstants.jl](https://github.com/JuliaPhysics/PhysicalConstants.jl), [MathPhysicalConstants.jl](https://github.com/LaGuer/MathPhysicalConstants.jl), [Unitful.jl](https://github.com/PainterQubits/Unitful.jl.git), [UnitfulUS.jl](https://github.com/PainterQubits/UnitfulUS.jl), [UnitfulAstro.jl](https://github.com/JuliaAstro/UnitfulAstro.jl), [UnitfulAtomic.jl](https://github.com/sostock/UnitfulAtomic.jl), [NaturallyUnitful.jl](https://github.com/MasonProtter/NaturallyUnitful.jl), and [UnitfulMoles.jl](https://github.com/rafaqz/UnitfulMoles.jl).
