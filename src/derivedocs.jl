
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
```Julia
julia> deka
$deka

julia> hecto
$hecto

julia> kilo
$kilo

julia> mega
$mega

julia> giga
$giga

julia> tera
$tera

julia> peta
$peta

julia> exa
$exa

julia> zetta
$zetta

julia> yotta
$yotta
```
""" deka,hecto,kilo,mega,giga,tera,peta,exa,zetta,yotta

@doc """
```Julia
julia> deci
$deci

julia> centi
$centi

julia> milli
$milli

julia> micro
$micro

julia> nano
$nano

julia> pico
$pico

julia> femto
$femto

julia> atto
$atto

julia> zepto
$zepto

julia> yocto
$yocto
```
""" deci,centi,milli,micro,nano,pico,femto,atto,zepto,yocto

@doc """
```Julia
julia> byte
$byte

julia> kibi
$kibi

julia> mebi
$mebi

julia> gibi
$gibi

julia> tebi
$tebi

julia> pebi
$pebi

julia> exbi
$exbi

julia> zebi
$zebi

julia> yobi
$yobi
```
""" byte,kibi,mebi,gibi,tebi,pebi,exbi,zebi,yobi

# angle

@doc """
$(unitext(:radian,"angle(𝟏,U,Metric)"))

Unit of `angle` which is dimensionless (rad).
```Julia
julia> radian(MetricEngineering) # rad
$(radian(MetricEngineering))

julia> radian(MetricDegree) # deg
$(radian(MetricDegree))

julia> radian(MetricArcminute) # amin
$(radian(MetricArcminute))

julia> radian(MetricArcsecond) # asec
$(radian(MetricArcsecond))

julia> radian(MetricGradian) # gon
$(radian(MetricGradian))
```
""" radian

@doc """
$(unitext(:degree,"angle(𝟏,U,MetricDegree)"))

Unit of `angle` which divides a `turn` into `360` parts (rad).
```Julia
julia> degree(MetricEngineering) # rad
$(degree(MetricEngineering))

julia> degree(MetricDegree) # deg
$(degree(MetricDegree))

julia> degree(MetricArcminute) # amin
$(degree(MetricArcminute))

julia> degree(MetricArcsecond) # asec
$(degree(MetricArcsecond))

julia> degree(MetricGradian) # gon
$(degree(MetricGradian))
```
""" degree

@doc """
$(unitext(:gradian,"angle(𝟏,U,MetricGradian)"))

Unit of `angle` which divides a `turn` into `400` parts (rad).
```Julia
julia> gradian(MetricEngineering) # rad
$(gradian(MetricEngineering))

julia> gradian(MetricDegree) # deg
$(gradian(MetricDegree))

julia> gradian(MetricArcminute) # amin
$(gradian(MetricArcminute))

julia> gradian(MetricArcsecond) # asec
$(gradian(MetricArcsecond))

julia> gradian(MetricGradian) # gon
$(gradian(MetricGradian))
```
""" gradian

@doc """
$(unitext(:bradian,"angle(τ/𝟐^8,U,Metric)"))

Unit of `angle` which divides a `turn` into `𝟐^8` or `256` parts (rad).
```Julia
julia> bradian(MetricEngineering) # rad
$(bradian(MetricEngineering))

julia> bradian(MetricDegree) # deg
$(bradian(MetricDegree))

julia> bradian(MetricArcminute) # amin
$(bradian(MetricArcminute))

julia> bradian(MetricArcsecond) # asec
$(bradian(MetricArcsecond))

julia> bradian(MetricGradian) # gon
$(bradian(MetricGradian))
```
""" bradian

@doc """
$(unitext(:arcminute,"angle(𝟏,U,MetricArcminute)"))

Unit of `angle` which divides a `degree` into `60` parts (rad).
```Julia
julia> arcminute(MetricEngineering) # rad
$(arcminute(MetricEngineering))

julia> arcminute(MetricDegree) # deg
$(arcminute(MetricDegree))

julia> arcminute(MetricArcminute) # amin
$(arcminute(MetricArcminute))

julia> arcminute(MetricArcsecond) # asec
$(arcminute(MetricArcsecond))

julia> arcminute(MetricGradian) # gon
$(arcminute(MetricGradian))
```
""" arcminute

@doc """
$(unitext(:arcsecond,"angle(𝟏,U,MetricArcsecond)"))

Unit of `angle` which divides a `arcminute` into `60` parts (rad).
```Julia
julia> arcsecond(MetricEngineering) # rad
$(arcsecond(MetricEngineering))

julia> arcsecond(MetricDegree) # deg
$(arcsecond(MetricDegree))

julia> arcsecond(MetricArcminute) # amin
$(arcsecond(MetricArcminute))

julia> arcsecond(MetricArcsecond) # asec
$(arcsecond(MetricArcsecond))

julia> arcsecond(MetricGradian) # gon
$(arcsecond(MetricGradian))
```
""" arcsecond

# solidangle

@doc """
$(unitext(:steradian,"solidangle(𝟏,U,Metric)"))

Unit of `solidangle` which is dimensionless (rad²).
```Julia
julia> steradian(MetricEngineering) # rad²
$(steradian(MetricEngineering))

julia> steradian(MetricDegree) # deg²
$(steradian(MetricDegree))

julia> steradian(MetricArcminute) # amin²
$(steradian(MetricArcminute))

julia> steradian(MetricArcsecond) # asec²
$(steradian(MetricArcsecond))

julia> steradian(MetricGradian) # gon²
$(steradian(MetricGradian))
```
""" steradian

@doc """
$(unitext(:squaredegree,"solidangle(𝟏,U,MetricDegree)"))

Unit of `solidangle` which is a `degree` squared (rad²).
```Julia
julia> squaredegree(MetricEngineering) # rad²
$(squaredegree(MetricEngineering))

julia> squaredegree(MetricDegree) # deg²
$(squaredegree(MetricDegree))

julia> squaredegree(MetricArcminute) # amin²
$(squaredegree(MetricArcminute))

julia> squaredegree(MetricArcsecond) # asec²
$(squaredegree(MetricArcsecond))

julia> squaredegree(MetricGradian) # gon²
$(squaredegree(MetricGradian))
```
""" squaredegree

# time

@doc """
$(unitext(:second,"time(𝟏,U,Metric)"))

Unit of `time` defined by `hyperfine` transition frequency of Cs-133 atom (s).
```Julia
julia> second(Metric) # s
$(second(Metric))

julia> second(MPH) # h
$(second(MPH))

julia> second(IAU) # D
$(second(IAU))
```
""" second

@doc """
$(unitext(:minute,"𝟐^2*𝟑*𝟓*second(U)"))

Unit of `time` defined by 60 `second` intervals (s).
```Julia
julia> minute(Metric) # s
$(minute(Metric))

julia> minute(MPH) # h
$(minute(MPH))

julia> minute(IAU) # D
$(minute(IAU))
```
""" minute

@doc """
$(unitext(:hour,"𝟐^2*𝟑*𝟓*minute(U)"))

Unit of `time` defined by 60 `minute` intervals (s).
```Julia
julia> hour(Metric) # s
$(hour(Metric))

julia> hour(MPH) # h
$(hour(MPH))

julia> hour(IAU) # D
$(hour(IAU))
```
""" hour, HOUR

@doc """
$(unitext(:day,"𝟐^3*𝟑*hour(U)"))

Unit of `time` defined by 24 `hour` intervals (s).
```Julia
julia> day(Metric) # s
$(day(Metric))

julia> day(MPH) # h
$(day(MPH))

julia> day(IAU) # D
$(day(IAU))
```
""" day, DAY

@doc """
$(unitext(:year,"aⱼ*day(U)"))

Unit of `time` defined by Julian calendar year interval (s).
```Julia
julia> year(Metric) # s
$(year(Metric))

julia> year(MPH) # h
$(year(MPH))

julia> year(IAU) # D
$(year(IAU))
```
""" year, aⱼ

@doc """
$(unitext(:radarmile,"𝟐*nauticalmile(U)/lightspeed(U)"))

Unit of `time` delay from a two-way `nauticalmile` radar return (s).
```Julia
julia> radarmile(Metric)
$(radarmile(Metric))
```
""" radarmile

# length

@doc """
$(unitext(:meter,"length(𝟏,U,Metric)"))

Metric unit of `length` (m or ft).
```Julia
julia> meter(CGS) # cm
$(meter(CGS))

julia> meter(English) # ft
$(meter(English))

julia> meter(Meridian) # em
$(meter(Meridian))
```
""" meter

@doc """
$(unitext(:earthmeter,"greatcircle(U)/𝟐^9/𝟓^7"))

Meridian unit of `length` as originally defined in France (m or ft).
```Julia
julia> earthmeter(CGS) # cm
$(earthmeter(CGS))

julia> earthmeter(English) # ft
$(earthmeter(English))

julia> earthmeter(Meridian) # em
$(earthmeter(Meridian))
```
""" earthmeter

@doc """
$(unitext(:angstrom,"hecto*pico*meter(U)"))

Metric unit of `length` (m or ft).
```Julia
julia> angstrom(CGS) # cm
$(angstrom(CGS))

julia> angstrom(English) # ft
$(angstrom(English))

julia> angstrom(IPS) # in
$(angstrom(IPS))
```
""" angstrom

@doc """
$(unitext(:foot,"length(𝟏,U,English)"))

English unit of `length` (m or ft).
```Julia
julia> foot(Metric) # m
$(foot(Metric))

julia> foot(Survey) # ftUS
$(foot(Survey))

julia> foot(IPS) # in
$(foot(IPS))
```
""" foot

@doc """
$(unitext(:surveyfoot,"length(𝟏,U,Survey)"))

Survey unit of `length` (m or ft).
```Julia
julia> surveyfoot(Metric) # m
$(surveyfoot(Metric))

julia> surveyfoot(English) # ft
$(surveyfoot(English))

julia> surveyfoot(IPS) # in
$(surveyfoot(IPS))
```
""" surveyfoot

@doc """
$(unitext(:inch,"length(𝟏,U,IPS)"))

English unit of `length` (m or ft).
```Julia
julia> inch(Metric) # m
$(inch(Metric))

julia> inch(English) # ft
$(inch(English))

julia> inch(IPS) # in
$(inch(IPS))
```
""" inch

#=@doc """
$(unitext(:rackunit,"length($(ft*𝟕/𝟐^4/𝟑),U,English)"))

Height unit of a 19-inch rack frame (m or ft).
```Julia
julia> rackunit(Metric) # m
$(rackunit(Metric))

julia> rackunit(CGS) # cm
$(rackunit(CGS))

julia> rackunit(English) # ft
$(rackunit(English))
```
""" rackunit=#

@doc """
$(unitext(:yard,"𝟑*foot(U)"))

English unit of `length` (m or ft).
```Julia
julia> yard(Metric) # m
$(yard(Metric))

julia> yard(English) # ft
$(yard(English))

julia> yard(IPS) # in
$(yard(IPS))
```
""" yard

@doc """
$(unitext(:statutemile,"length(𝟐^5*𝟑*𝟓*𝟏𝟏,U,Survey)"))

Statute `Survey` mile (m or ft).
```Julia
julia> statutemile(Metric) # m
$(statutemile(Metric))

julia> statutemile(English) # ft
$(statutemile(English))

julia> statutemile(Survey) # ftUS
$(statutemile(Survey))
```
""" statutemile

@doc """
$(unitext(:astronomicalunit,"length(𝟏,U,IAU)"))

Standard astronomical unit from the International Astronomical Union (m or ft).
```Julia
julia> astronomicalunit(Metric) # m
$(astronomicalunit(Metric))

julia> astronomicalunit(English) # ft
$(astronomicalunit(English))

julia> astronomicalunit(Metric)/lightspeed(Metric) # s
$(astronomicalunit(Metric)/lightspeed(Metric))
```
""" astronomicalunit, au

@doc """
$(unitext(:lunardistance,"length(𝟏,U,IAUE)"))

Standard distance between the Earth and the Moon (m or ft).
```Julia
julia> lunardistance(Metric) # m
$(lunardistance(Metric))

julia> lunardistance(Nautical) # nm
$(lunardistance(Nautical))

julia> lunardistance(Metric)/lightspeed(Metric) # s
$(lunardistance(Metric)/lightspeed(Metric))
```
""" lunardistance, LD

@doc """
$(unitext(:jupiterdistance,"length(𝟏,U,IAUJ)"))

Standard distance between the Sun and the planet Jupiter (m or ft).
```Julia
julia> jupiterdistance(Metric) # m
$(jupiterdistance(Metric))

julia> jupiterdistance(IAU) # au
$(jupiterdistance(IAU))

julia> jupiterdistance(Metric)/lightspeed(Metric) # s
$(jupiterdistance(Metric)/lightspeed(Metric))
```
""" jupiterdistance, JD

@doc """
$(unitext(:mile,"length(𝟏,U,MPH)"))

Statute `English` mile (m or ft).
```Julia
julia> mile(Metric) # m
$(mile(Metric))

julia> mile(English) # ft
$(mile(English))

julia> mile(Nautical) # nm
$(mile(Survey))
```
""" mile

@doc """
$(unitext(:admiraltymile,"length(𝟐^6*𝟓*𝟏𝟗,U,English)"))

Historic nautical mile as defined by the Clarke authalic radius (m or ft).
```Julia
julia> admiraltymile(Metric) # m
$(admiraltymile(Metric))

julia> admiraltymile(English) # ft
$(admiraltymile(English))

julia> admiraltymile(Nautical) # nm
$(admiraltymile(Nautical))
```
""" admiraltymile

@doc """
$(unitext(:meridianmile,"length(𝟐^4*𝟓^5/𝟑^3,U,Metric)"))

Historic nautical mile as defined by naive meridian assumption (m or ft).
```Julia
julia> meridianmile(Metric) # m
$(meridianmile(Metric))

julia> meridianmile(English) # ft
$(meridianmile(English))

julia> meridianmile(Nautical) # nm
$(meridianmile(Nautical))
```
""" meridianmile

@doc """
$(unitext(:nauticalmile,"greatcircle(U)/𝟐^5/𝟑^3/𝟓^2"))

Standard `nauticalmile` as defined by `earthradius` (m or ft).
```Julia
julia> nauticalmile(Metric) # m
$(nauticalmile(Metric))

julia> nauticalmile(Meridian) # em
$(nauticalmile(Meridian))

julia> nauticalmile(English) # ft
$(nauticalmile(English))
```
""" nauticalmile, nm

@doc """
$(unitext(:lightyear,"year(U)*lightspeed(U)"))

Unit of `length` defined by distance traveled by light in 1 `year` unit.
```Julia
julia> lightyear(Metric) # m
$(lightyear(Metric))

julia> lightyear(English) # ft
$(lightyear(English))

julia> lightyear(IAU) # au
$(lightyear(IAU))
```
""" lightyear, ly

@doc """
$(unitext(:parsec,"astronomicalunit(U)*𝟐^2*𝟑^4*𝟓^3/τ"))

Unit of `length` defined at which 1 `astronomicalunit` subtends an angle of 1 arcsecond.
```Julia
julia> parsec(Metric) # m
$(parsec(Metric))

julia> parsec(English) # ft
$(parsec(English))

julia> parsec(IAU) # au
$(parsec(IAU))
```
""" parsec, pc

# area

@doc """
$(unitext(:barn,"area((𝟐*𝟓)^-28,U,Metric)"))

Unit of `area` defined by `100` square femto-meters (m² or ft²).
```Julia
julia> barn(Metric) # m²
$(barn(Metric))

julia> barn(CGS) # cm²
$(barn(CGS))

julia> barn(English) # ft²
$(barn(English))
```
""" barn

@doc """
$(unitext(:hectare,"area(hecto^2,U,Metric)"))

Metric unit of land `area` defined by `100` square meters (m² or ft²).
```Julia
julia> hectare(Metric) # m²
$(hectare(Metric))

julia> hectare(English) # ft²
$(hectare(English))

julia> hectare(Survey) # ftUS²
$(hectare(Survey))
```
""" hectare

@doc """
$(unitext(:acre,"area(𝟐^4*𝟓^4,U,Metric)"))

English unit of land `area` (m² or ft²).
```Julia
julia> acre(Metric) # m²
$(acre(Metric))

julia> acre(English) # ft²
$(acre(English))

julia> acre(Survey) # ftUS²
$(acre(Survey))
```
""" acre

@doc """
$(unitext(:surveyacre,"area(𝟐^3*𝟑^2*𝟓*𝟏𝟏^2,U,Survey)"))

Survey unit of land `area` (m² or ft²).
```Julia
julia> surveyacre(Metric) # m²
$(surveyacre(Metric))

julia> surveyacre(English) # ft²
$(surveyacre(English))

julia> surveyacre(Survey) # ftUS²
$(surveyacre(Survey))
```
""" surveyacre

# volume

@doc """
$(unitext(:gallon,"volume(𝟕*𝟏𝟏/𝟐^2,U,English)"))

Unit of `volume` derived from the US liquid `gallon` in cubic inches (m³ or ft³).
```Julia
julia> gallon(Metric) # m³
$(gallon(Metric))

julia> gallon(CGS) # cm³
$(gallon(CGS))

julia> gallon(IPS) # in³
$(gallon(IPS))
```
""" gallon, gal

@doc """
$(unitext(:liter,"volume(𝟏𝟎^-3,U,Metric)"))

Unit of `volume` derived from 1 cubic decimeter (m³ or ft³).
```Julia
julia> liter(Metric) # m³
$(liter(Metric))

julia> liter(CGS) # cm³
$(liter(CGS))

julia> liter(IPS) # in³
$(liter(IPS))
```
""" liter

@doc """
$(unitext(:quart,"gallon(U)/𝟐^2"))

English unit of `volume` (m³ or ft³).
```Julia
julia> quart(Metric) # m³
$(quart(Metric))

julia> quart(CGS) # cm³
$(quart(CGS))

julia> quart(IPS) # in³
$(quart(IPS))
```
""" quart

@doc """
$(unitext(:pint,"quart(U)/𝟐"))

English unit of `volume` (m³ or ft³).
```Julia
julia> pint(Metric) # m³
$(pint(Metric))

julia> pint(CGS) # cm³
$(pint(CGS))

julia> pint(IPS) # in³
$(pint(IPS))
```
""" pint

@doc """
$(unitext(:cup,"pint(U)/𝟐"))

English unit of `volume` (m³ or ft³).
```Julia
julia> cup(Metric) # m³
$(cup(Metric))

julia> cup(CGS) # cm³
$(cup(CGS))

julia> cup(IPS) # in³
$(cup(IPS))
```
""" cup

@doc """
$(unitext(:fluidounce,"cup(U)/𝟐^3"))

English unit of `volume` (m³ or ft³).
```Julia
julia> fluidounce(Metric) # m³
$(fluidounce(Metric))

julia> fluidounce(CGS) # cm³
$(fluidounce(CGS))

julia> fluidounce(IPS) # in³
$(fluidounce(IPS))
```
""" fluidounce

@doc """
$(unitext(:teaspoon,"𝟓*milli*liter(U)"))

Measuring `teaspoon` unit of `volume` (m³ or ft³).
```Julia
julia> teaspoon(Metric) # m³
$(teaspoon(Metric))

julia> teaspoon(CGS) # cm³
$(teaspoon(CGS))

julia> teaspoon(IPS) # in³
$(teaspoon(IPS))
```
""" teaspoon

@doc """
$(unitext(:tablespoon,"𝟑*teaspoon(U)"))

Measuring `tablespoon` unit of `volume` (m³ or ft³).
```Julia
julia> tablespoon(Metric) # m³
$(tablespoon(Metric))

julia> tablespoon(CGS) # cm³
$(tablespoon(CGS))

julia> tablespoon(IPS) # in³
$(tablespoon(IPS))
```
""" tablespoon

# speed

@doc """
$(unitext(:bubnoff,"meter(U)/year(U)"))

Reference unit of erosion `speed` (m⋅s⁻¹ or ft⋅s⁻¹).
```Julia
julia> bubnoff(CGS) # cm⋅s⁻¹
$(bubnoff(CGS))

julia> bubnoff(English) # ft⋅s⁻¹
$(bubnoff(English))
```
""" bubnoff

@doc """
$(unitext(:ips,"inch(U)/second(U)"))

Inch per second unit of `speed` (m⋅s⁻¹ or ft⋅s⁻¹).
```Julia
julia> ips(CGS) # cm⋅s⁻¹
$(ips(CGS))

julia> ips(English) # ft⋅s⁻¹
$(ips(English))
```
""" ips

@doc """
$(unitext(:fps,"feet(U)/second(U)"))

Feet per second unit of `speed` (m⋅s⁻¹ or ft⋅s⁻¹).
```Julia
julia> fps(Metric) # m⋅s⁻¹
$(fps(Metric))

julia> fps(KKH) # km⋅h⁻¹
$(fps(KKH))

julia> fps(MPH) # mi⋅h⁻¹
$(fps(MPH))
```
""" fps

@doc """
$(unitext(:fpm,"feet(U)/minute(U)"))

Feet per minute unit of `speed` (m⋅s⁻¹ or ft⋅s⁻¹).
```Julia
julia> fpm(CGS) # cm⋅s⁻¹
$(fpm(CGS))

julia> fpm(IPS) # in⋅s⁻¹
$(fpm(IPS))

julia> fpm(English) # ft⋅s⁻¹
$(fpm(English))
```
""" fpm

@doc """
$(unitext(:ms,"meter(U)/second(U)"))

Meters per second unit of `speed` (m⋅s⁻¹ or ft⋅s⁻¹).
```Julia
julia> ms(KKH) # km⋅h⁻¹
$(ms(KKH))

julia> ms(MPH) # mi⋅h⁻¹
$(ms(MPH))

julia> ms(Nautical) # nm⋅h⁻¹
$(ms(Nautical))
```
""" ms

@doc """
$(unitext(:kmh,"kilo(U)*meter(U)/hour(U)"))

Kilometers per hour unit of `speed` (m⋅s⁻¹ or ft⋅s⁻¹).
```Julia
julia> kmh(Metric) # m⋅s⁻¹
$(kmh(Metric))

julia> kmh(MPH) # mi⋅h⁻¹
$(kmh(MPH))

julia> kmh(Nautical) # nm⋅h⁻¹
$(kmh(Nautical))
```
""" kmh

@doc """
$(unitext(:mph,"mile(U)/hour(U)"))

Miles per hour unit of `speed` (m⋅s⁻¹ or ft⋅s⁻¹).
```Julia
julia> mph(Metric) # m⋅s⁻¹
$(mph(Metric))

julia> mph(KKH) # km⋅h⁻¹
$(mph(KKH))

julia> mph(Nautical) # nm⋅h⁻¹
$(mph(Nautical))
```
""" mph

@doc """
$(unitext(:knot,"nauticalmile(U)/hour(U)"))

Nautical miles per hour unit of `speed` (m⋅s⁻¹ or ft⋅s⁻¹).
```Julia
julia> knot(Metric) # m⋅s⁻¹
$(knot(Metric))

julia> knot(KKH) # km⋅h⁻¹
$(knot(KKH))

julia> knot(MPH) # mi⋅h⁻¹
$(knot(MPH))
```
""" knot

@doc """
$(unitext(:mps,"mile(U)/second(U)"))

Miles per second unit of `speed` (m⋅s⁻¹ or ft⋅s⁻¹).
```Julia
julia> mps(KKH) # km⋅h⁻¹
$(mps(KKH))

julia> mps(MPH) # mi⋅h⁻¹
$(mps(MPH))

julia> mps(Nautical) # nm⋅h⁻¹
$(mps(Nautical))
```
""" mps

# mass

@doc """
$(unitext(:grain,"milli(U)*pound(U)/𝟕"))

Ideal `grain` seed of cereal, unit of `mass` (kg or lb).
```Julia
julia> grain(Metric) # kg
$(grain(Metric))

julia> grain(CGS) # g
$(grain(CGS))

julia> grain(English) # lb
$(grain(English))

julia> grain(QCD) # mₚ
$(grain(QCD))
```
""" grain

@doc """
$(unitext(:gram,"mass(𝟏,U,Gauss)"))

Metric `gram` unit of `mass` (kg or lb).
```Julia
julia> gram(Metric) # kg
$(gram(Metric))

julia> gram(CGS) # g
$(gram(CGS))

julia> gram(English) # lb
$(gram(English))

julia> gram(British) # slug
$(gram(British))

julia> gram(GravitationalMetric) # hyl
$(gram(GravitationalMetric))
```
""" gram

@doc """
$(unitext(:earthgram,"mass(milli,U,Meridian)"))

Meridian `gram` unit of `mass` based on `earthmeter` (kg or lb).
```Julia
julia> earthgram(Meridian) # keg
$(earthgram(Meridian))

julia> earthgram(CGS) # g
$(earthgram(CGS))

julia> earthgram(English) # lb
$(earthgram(English))

julia> earthgram(British) # slug
$(earthgram(British))

julia> earthgram(GravitationalMetric) # hyl
$(earthgram(GravitationalMetric))
```
""" earthgram

@doc """
$(unitext(:kilogram,"mass(𝟏,U,Metric)"))

Metric `kilogram` unit of `mass` (kg or lb).
```Julia
julia> kilogram(Metric) # kg
$(kilogram(Metric))

julia> kilogram(CGS) # g
$(kilogram(CGS))

julia> kilogram(English) # lb
$(kilogram(English))

julia> kilogram(British) # slug
$(kilogram(British))

julia> kilogram(GravitationalMetric) # hyl
$(kilogram(GravitationalMetric))
```
""" kilogram

@doc """
$(unitext(:tonne,"mass(𝟏,U,MTS)"))

Metric `tonne` unit of `mass` (kg or lb).
```Julia
julia> tonne(Metric) # kg
$(tonne(Metric))

julia> tonne(MTS) # t
$(tonne(MTS))

julia> tonne(English) # lb
$(tonne(English))

julia> tonne(British) # slug
$(tonne(British))

julia> tonne(GravitationalMetric) # hyl
$(tonne(GravitationalMetric))
```
""" tonne

@doc """
$(unitext(:ton,"mass(𝟐*kilo,U,English)"))

English `ton` unit of `mass` (kg or lb).
```Julia
julia> ton(Metric) # kg
$(ton(Metric))

julia> ton(MTS) # t
$(ton(MTS))

julia> ton(English) # lb
$(ton(English))

julia> ton(British) # slug
$(ton(British))

julia> ton(GravitationalMetric) # hyl
$(ton(GravitationalMetric))
```
""" ton

@doc """
$(unitext(:pound,"mass(𝟏,U,English)"))

English `pound` unit of `mass` (kg or lb).
```Julia
julia> pound(Metric) # kg
$(pound(Metric))

julia> pound(CGS) # g
$(pound(CGS))

julia> pound(English) # lb
$(pound(English))

julia> pound(British) # slug
$(pound(British))

julia> pound(GravitationalMetric) # hyl
$(pound(GravitationalMetric))
```
""" pound

@doc """
$(unitext(:ounce,"pound(U)/𝟐^4"))

English `ounce` unit of `mass` (kg or lb).
```Julia
julia> ounce(Metric) # kg
$(ounce(Metric))

julia> ounce(CGS) # g
$(ounce(CGS))

julia> ounce(English) # lb
$(ounce(English))

julia> ounce(British) # slug
$(ounce(British))

julia> ounce(GravitationalMetric) # hyl
$(ounce(GravitationalMetric))
```
""" ounce

@doc """
$(unitext(:slug,"mass(𝟏,U,British)"))

British gravitational `slug` unit of `mass` (kg or lb).
```Julia
julia> slug(Metric) # kg
$(slug(Metric))

julia> slug(CGS) # g
$(slug(CGS))

julia> slug(English) # lb
$(slug(English))

julia> slug(British) # slug
$(slug(British))

julia> slug(GravitationalMetric) # hyl
$(slug(GravitationalMetric))
```
""" slug

@doc """
$(unitext(:slinch,"mass(𝟏,U,IPS)"))

British gravitational `slinch` unit of `mass` (kg or lb).
```Julia
julia> slinch(Metric) # kg
$(slinch(Metric))

julia> slinch(CGS) # g
$(slinch(CGS))

julia> slinch(English) # lb
$(slinch(English))

julia> slinch(British) # slug
$(slinch(British))

julia> slinch(GravitationalMetric) # hyl
$(slinch(GravitationalMetric))
```
""" slinch

@doc """
$(unitext(:hyl,"mass(𝟏,U,GravitationalMetric)"))

Gravitational Metric `hyl` unit of `mass` (kg or lb).
```Julia
julia> hyl(Metric) # kg
$(hyl(Metric))

julia> hyl(CGS) # g
$(hyl(CGS))

julia> hyl(English) # lb
$(hyl(English))

julia> hyl(British) # slug
$(hyl(British))

julia> hyl(GravitationalMetric) # hyl
$(hyl(GravitationalMetric))
```
""" hyl

# force

@doc """
$(unitext(:dyne,"force(𝟏,U,Gauss)"))

Historical `dyne` unit of `force` (N or lb).
```Julia
julia> dyne(Metric) # N
$(dyne(Metric))

julia> dyne(CGS) # dyn
$(dyne(CGS))

julia> dyne(English) # lb
$(dyne(English))

julia> dyne(FPS) # pdl
$(dyne(FPS))

julia> dyne(MetricEngineering) # kp
$(dyne(MetricEngineering))
```
""" dyne

@doc """
$(unitext(:newton,"force(𝟏,U,Metric)"))

Metric `newton` unit of `force` (N or lb).
```Julia
julia> newton(Metric) # N
$(newton(Metric))

julia> newton(CGS) # dyn
$(newton(CGS))

julia> newton(English) # lb
$(newton(English))

julia> newton(FPS) # pdl
$(newton(FPS))

julia> newton(MetricEngineering) # kp
$(newton(MetricEngineering))
```
""" newton

@doc """
$(unitext(:poundal,"force(𝟏,U,FPS)"))

Absolute English `poundal` unit of `force` (N or lb).
```Julia
julia> poundal(Metric) # N
$(poundal(Metric))

julia> poundal(CGS) # dyn
$(poundal(CGS))

julia> poundal(English) # lb
$(poundal(English))

julia> poundal(FPS) # pdl
$(poundal(FPS))

julia> poundal(MetricEngineering) # kp
$(poundal(MetricEngineering))
```
""" poundal

@doc """
    kilopond(U::UnitSystem) = force(𝟏,U,MetricEngineering)

Gravitational `kilopond` unit of `force` used in engineering systems (N or lb).
```Julia
julia> kilopond(Metric) # N
$(kilopond(Metric))

julia> kilopond(CGS) # dyn
$(kilopond(CGS))

julia> kilopond(English) # lb
$(kilopond(English))

julia> kilopond(FPS) # pdl
$(kilopond(FPS))

julia> kilopond(MetricEngineering) # kp
$(kilopond(MetricEngineering))
```
""" kilopond

@doc """
$(unitext(:poundforce,"force(𝟏,U,English)"))

English `poundforce` unit of `force` used in engineering systems (N or lb).
```Julia
julia> poundforce(Metric) # N
$(poundforce(Metric))

julia> poundforce(CGS) # dyn
$(poundforce(CGS))

julia> poundforce(English) # lb
$(poundforce(English))

julia> poundforce(FPS) # pdl
$(poundforce(FPS))

julia> poundforce(MetricEngineering) # kp
$(poundforce(MetricEngineering))
```
""" poundforce

# pressure

@doc """
$(unitext(:pascal,"pressure(𝟏,U,Metric)"))

Metric unit of `pressure` (Pa or lb⋅ft⁻²).
```Julia
julia> pascal(Metric) # Pa
$(pascal(Metric))

julia> pascal(English) # lb⋅ft⁻²
$(pascal(English))

julia> pascal(IPS) # lb⋅in⁻²
$(pascal(IPS))
```
""" pascal

@doc """
$(unitext(:bar,"pressure(hecto*kilo,U,Metric)"))

Reference unit of `pressure` (Pa or lb⋅ft⁻²).
```Julia
julia> bar(Metric) # Pa
$(bar(Metric))

julia> bar(English) # lb⋅ft⁻²
$(bar(English))

julia> bar(IPS) # lb⋅in⁻²
$(bar(IPS))
```
""" bar

@doc """
$(unitext(:barye,"pressure(𝟏,U,Gauss)"))

Historical unit of `pressure` (Pa or lb⋅ft⁻²).
```Julia
julia> barye(Metric) # Pa
$(barye(Metric))

julia> barye(English) # lb⋅ft⁻²
$(barye(English))

julia> barye(IPS) # lb⋅in⁻²
$(barye(IPS))
```
""" barye

@doc """
$(unitext(:psi,"pressure(𝟏,U,IPS)"))

English unit of `pressure` (Pa or lb⋅ft⁻²).
```Julia
julia> psi(Metric) # Pa
$(psi(Metric))

julia> psi(English) # lb⋅ft⁻²
$(psi(English))

julia> psi(IPS) # lb⋅in⁻²
$(psi(IPS))
```
""" psi

@doc """
$(unitext(:technicalatmosphere,"kilopond(U)/(centi*meter(U))^2"))

Gravitational Metric unit of `pressure` (Pa or lb⋅ft⁻²).
```Julia
julia> technicalatmosphere(Metric) # Pa
$(technicalatmosphere(Metric))

julia> technicalatmosphere(English) # lb⋅ft⁻²
$(technicalatmosphere(English))

julia> technicalatmosphere(IPS) # lb⋅in⁻²
$(technicalatmosphere(IPS))
```
""" technicalatmosphere

@doc """
$(unitext(:atmosphere,"pressure($atm,U)"))

Standard `pressure` reference level of one atmosphere `atm` (Pa or lb⋅ft⁻²).
```Julia
julia> atmosphere(Metric) # Pa
$(atmosphere(Metric))

julia> atmosphere(English) # lb⋅ft⁻²
$(atmosphere(English))

julia> atmosphere(IPS) # lb⋅in⁻²
$(atmosphere(IPS))
```
""" atmosphere, atm

@doc """
$(unitext(:inchmercury,"pressure(inHg,U,Metric)"))

Unit of `pressure` exerted by 1 inch of mercury at standard atmospheric conditions.
```Julia
juila> inchmercury(Metric) # Pa
$(inchmercury(Metric))

julia> inchmercury(English) # lb⋅ft⁻²
$(inchmercury(English))

julia> inchmercury(IPS) # lb⋅in⁻²
$(inchmercury(IPS))
```
""" inchmercury, inHg

@doc """
$(unitext(:torr,"pressure(atm/𝟐^3/𝟓/𝟏𝟗,U,Metric)"))

Unit of `pressure` exerted by 1 mm of mercury at standard atmospheric conditions.
```Julia
juila> torr(Metric) # Pa
$(torr(Metric))

julia> torr(English) # lb⋅ft⁻²
$(torr(English))

julia> torr(IPS) # lb⋅in⁻²
$(torr(IPS))
```
""" torr

# energy

@doc """
$(unitext(:erg,"energy(𝟏,U,Gauss)"))

Historical unit of `energy` (J or lb⋅ft).
```Julia
julia> erg(Metric) # J
$(erg(Metric))

julia> erg(CGS) # erg
$(erg(CGS))

julia> erg(British) # lb⋅ft
$(erg(British))
```
""" erg

@doc """
$(unitext(:joule,"energy(𝟏,U,Metric)"))

Metric unit of `energy` (J or lb⋅ft).
```Julia
julia> joule(Metric) # J
$(joule(Metric))

julia> joule(CGS) # erg
$(joule(CGS))

julia> joule(British) # lb⋅ft
$(joule(British))
```
""" joule

@doc """
$(unitext(:tontnt,"giga*calorie(U)"))

Ton TNT equivalent reference unit of `energy` (J or lb⋅ft).
```Julia
julia> tontnt(Metric) # J
$(tontnt(Metric))

julia> tontnt(CGS) # erg
$(tontnt(CGS))

julia> tontnt(British) # lb⋅ft
$(tontnt(British))
```
""" tontnt

@doc """
$(unitext(:gasgallon,"𝟐*𝟑*𝟏𝟗*kilo*thermalunit(U)"))

Gasoline gallon equivalent reference unit of `energy` (J or lb⋅ft).
```Julia
julia> gasgallon(Metric) # J
$(gasgallon(Metric))

julia> gasgallon(CGS) # erg
$(gasgallon(CGS))

julia> gasgallon(British) # lb⋅ft
$(gasgallon(British))
```
""" gasgallon

@doc """
$(unitext(:footpound,"poundforce(U)*foot(U)"))

English unit of `energy` in gravitational and engineering systems (J or lb⋅ft).
```Julia
julia> footpound(Metric) # J
$(footpound(Metric))

julia> footpound(CGS) # erg
$(footpound(CGS))

julia> footpound(British) # lb⋅ft
$(footpound(British))
```
""" footpound

@doc """
$(unitext(:kilocalorie,"energy(𝟐^5*𝟓^4*𝟑^2/𝟒𝟑,U,International)"))

Heat energy required to raise 1 kg of water by 1 Kelvin (`kcal`) in `International` units.
```Julia
julia> kilocalorie(International) # J
$(kilocalorie(International))

julia> kilocalorie(Metric) # J
$(kilocalorie(Metric))

julia> kilocalorie(English) # ft⋅lb
$(kilocalorie(English))
```
""" kilocalorie, kcal

@doc """
$(unitext(:calorie,"kilocalorie(U)/𝟐^3/𝟓^3"))

Heat energy required to raise 1 g of water by 1 Kelvin (`cal`) in `International` units.
```Julia
julia> calorie(International) # J
$(calorie(International))

julia> calorie(Metric) # J
$(calorie(Metric))

julia> calorie(English) # ft⋅lb
$(calorie(English))
```
""" calorie, cal

@doc """
$(unitext(:meancalorie,"energy(𝟐^2*𝟓*𝟑^2/𝟒𝟑,U,InternationalMean)"))

Heat energy required to raise 1 g of water by 1 Kelvin (`cal`) in `InternationalMean` units.
```Julia
julia> meancalorie(InternationalMean) # J
$(meancalorie(InternationalMean))

julia> meancalorie(Metric) # J
$(meancalorie(Metric))

julia> meancalorie(English) # ft⋅lb
$(meancalorie(English))
```
""" meancalorie

@doc """
$(unitext(:earthcalorie,"calorie(U)*(sqrt(g₀/GME)/τ)^3*𝟐^27*𝟓^21"))

Heat energy required to raise 1 `earthgram` of water by 1 `kelvin` in `Meridian` units.
```Julia
julia> earthcalorie(Meridian) # J
$(earthcalorie(Meridian))

julia> earthcalorie(Metric) # J
$(earthcalorie(Metric))

julia> earthcalorie(British) # ft⋅lb
$(earthcalorie(British))
```
""" earthcalorie

@doc """
$(unitext(:thermalunit,"kilocalorie(U)*𝟑^2/𝟓/lb"))

Heat energy required to raise 1 lb of water by 1 Rankine (`BTU`) in `International` units.
```Julia
julia> thermalunit(British) # J
$(thermalunit(British))

julia> thermalunit(International) # J
$(thermalunit(International))

julia> thermalunit(Metric) # ft⋅lb
$(thermalunit(Metric))
```
""" thermalunit, BTU, BTUJ, BTUftlb

@doc """
$(unitext(:electronvolt,"elementarycharge(U)*volt(U)"))

Unit of `energy` gained by a rest electron accelerated by 1 `volt` in vacuum (J or lb⋅ft).
```Julia
julia> electronvolt(SI2019) # J
$(electronvolt(SI2019))

julia> electronvolt(SI2019)/lightspeed(SI2019) # kg⋅m⋅s⁻¹
$(electronvolt(SI2019)/lightspeed(SI2019))

julia> electronvolt(SI2019)/lightspeed(SI2019)^2 # kg
$(electronvolt(SI2019)/lightspeed(SI2019)^2)

julia> electronvolt(SI2019)/planck(SI2019)/lightspeed(SI2019) # m⁻¹
$(electronvolt(SI2019)/planck(SI2019)/lightspeed(SI2019))

julia> electronvolt(SI2019)/boltzmann(SI2019) # K
$(electronvolt(SI2019)/boltzmann(SI2019))
```
""" electronvolt, eV

# power

@doc """
$(unitext(:watt,"power(𝟏,U,Metric)"))

Metric `watt` unit of `power` (W or lb⋅ft⋅s⁻¹).
```Julia
julia> watt(Metric) # W
$(watt(Metric))

julia> watt(English) # lb⋅ft⋅s⁻¹
$(watt(English))

julia> watt(MetricEngineering) # kgf⋅m⋅s⁻¹
$(watt(MetricEngineering))
```
""" watt

@doc """
$(unitext(:tonsrefrigeration,"frequency(𝟐*𝟓/𝟑,U,Metric)*thermalunit(U)"))

Unit of `power` derived from melting of 1 short ton of ice in 24 hours.
```Julia
julia> tonsrefrigeration(British) # lb⋅ft⋅s⁻¹
$(tonsrefrigeration(British))

julia> tonsrefrigeration(Metric) # W
$(tonsrefrigeration(Metric))

julia> tonsrefrigeration(MetricEngineering) # kgf⋅m⋅s⁻¹
$(tonsrefrigeration(MetricEngineering))
```
""" tonsrefrigeration

@doc """
$(unitext(:boilerhorsepower,"frequency(1339/𝟐^4/𝟑^2,U,Metric)*thermalunit(U)"))

Unit of `power` derived from evaporating 34.5 lb of boiling water in 1 hour.
```Julia
julia> boilerhorsepower(British) # lb⋅ft⋅s⁻¹
$(boilerhorsepower(British))

julia> boilerhorsepower(Metric) # W
$(boilerhorsepower(Metric))

julia> boilerhorsepower(MetricEngineering) # kgf⋅m⋅s⁻¹
$(boilerhorsepower(MetricEngineering))
```
""" boilerhorsepower

@doc """
$(unitext(:horsepower,"power(𝟐*𝟓^2*𝟏𝟏,U,British)"))

Unit of `power` derived from raising 550 lb by 1 ft in 1  in 1 s.
```Julia
julia> horsepower(British) # lb⋅ft⋅s⁻¹
$(horsepower(British))

julia> horsepower(Metric) # W
$(horsepower(Metric))

julia> horsepower(MetricEngineering) # kgf⋅m⋅s⁻¹
$(horsepower(MetricEngineering))
```
""" horsepower, HP

@doc """
$(unitext(:horsepowerwatt,"power(𝟐^4*𝟑^3/𝟓*τ,U,British)"))

Unit of `power` derived from Watt's exact original horse power estimate.
```Julia
julia> horsepowerwatt(British) # lb⋅ft⋅s⁻¹
$(horsepowerwatt(British))

julia> horsepowerwatt(Metric) # W
$(horsepowerwatt(Metric))

julia> horsepowerwatt(MetricEngineering) # kgf⋅m⋅s⁻¹
$(horsepowerwatt(MetricEngineering))
```
""" horsepowerwatt

@doc """
$(unitext(:horsepowermetric,"power(𝟑*𝟓^2,U,GravitationalMetric)"))

Unit of `power` derived from raising 75 kp by 1 m in 1  in 1 s.
```Julia
julia> horsepowermetric(British) # lb⋅ft⋅s⁻¹
$(horsepowermetric(British))

julia> horsepowermetric(Metric) # W
$(horsepowermetric(Metric))

julia> horsepowermetric(MetricEngineering) # kgf⋅m⋅s⁻¹
$(horsepowermetric(MetricEngineering))
```
""" horsepowermetric

@doc """
$(unitext(:electricalhorsepower,"power(746,U,Metric)"))

Unit of `power` for electrical motors in the United States.
```Julia
julia> electricalhorsepower(British) # lb⋅ft⋅s⁻¹
$(electricalhorsepower(British))

julia> electricalhorsepower(Metric) # W
$(electricalhorsepower(Metric))

julia> electricalhorsepower(MetricEngineering) # kgf⋅m⋅s⁻¹
$(electricalhorsepower(MetricEngineering))
```
""" electricalhorsepower

# electromagnetic

@doc """
$(unitext(:coulomb,"charge(𝟏,U,Metric)"))

Metric unit of `charge` (C).
```Julia
julia> coulomb(Metric) # C
$(coulomb(Metric))

julia> coulomb(EMU) # abC
$(coulomb(EMU))

julia> coulomb(ESU) # statC
$(coulomb(ESU))
```
""" coulomb

@doc """
$(unitext(:ampere,"current(𝟏,U,Metric)"))

Metric unit of `current` (C⋅s⁻¹).
```Julia
julia> ampere(Metric) # C⋅s⁻¹
$(ampere(Metric))

julia> ampere(EMU) # abC⋅s⁻¹
$(ampere(EMU))

julia> ampere(ESU) # statC⋅s⁻¹
$(ampere(ESU))
```
""" ampere

@doc """
$(unitext(:volt,"electricpotential(𝟏,U,Metric)"))

Metric unit of `electricpotential` (V).
```Julia
julia> volt(Metric) # V
$(volt(Metric))

julia> volt(EMU) # abV
$(volt(EMU))

julia> volt(ESU) # statV
$(volt(ESU))
```
""" volt

@doc """
$(unitext(:henry,"inductance(𝟏,U,Metric)"))

Metric unit of `inductance` (H).
```Julia
julia> henry(Metric) # H
$(henry(Metric))

julia> henry(EMU) # abH
$(henry(EMU))

julia> henry(ESU) # statH
$(henry(ESU))
```
""" henry

@doc """
$(unitext(:ohm,"resistance(𝟏,U,Metric)"))

Metric unit of `resistance` (Ω).
```Julia
julia> ohm(Metric) # Ω
$(ohm(Metric))

julia> ohm(EMU) # abΩ
$(ohm(EMU))

julia> ohm(ESU) # statΩ
$(ohm(ESU))
```
""" ohm

@doc """
$(unitext(:siemens,"conductance(𝟏,U,Metric)"))

Metric unit of `conductance` (S).
```Julia
julia> siemens(Metric) # S
$(siemens(Metric))

julia> siemens(EMU) # abS
$(siemens(EMU))

julia> siemens(ESU) # statS
$(siemens(ESU))
```
""" siemens

@doc """
$(unitext(:farad,"capacitance(𝟏,U,Metric)"))

Metric unit of `capacitance` (F).
```Julia
julia> farad(Metric) # F
$(farad(Metric))

julia> farad(EMU) # abF
$(farad(EMU))

julia> farad(ESU) # statF
$(farad(ESU))
```
""" farad

@doc """
$(unitext(:weber,"magneticflux(𝟏,U,Metric)"))

Metric unit of `magneticflux` (Wb).
```Julia
julia> weber(Metric) # Wb
$(weber(Metric))

julia> weber(EMU) # Mx
$(weber(EMU))

julia> weber(ESU) # statWb
$(weber(ESU))
```
""" weber

@doc """
$(unitext(:tesla,"magneticfluxdensity(𝟏,U,Metric)"))

Metric unit of `magneticfluxdensity` (T).
```Julia
julia> tesla(Metric) # T
$(tesla(Metric))

julia> tesla(EMU) # G
$(tesla(EMU))

julia> tesla(ESU) # statT
$(tesla(ESU))
```
""" tesla

@doc """
$(unitext(:abcoulomb,"charge(𝟏,U,EMU)"))

Electromagnetic unit of `charge` (C).
```Julia
julia> abcoulomb(Metric) # C
$(abcoulomb(Metric))

julia> abcoulomb(EMU) # abC
$(abcoulomb(EMU))

julia> abcoulomb(ESU) # statC
$(abcoulomb(ESU))
```
""" abcoulomb

@doc """
$(unitext(:abampere,"current(𝟏,U,EMU)"))

Electromagnetic unit of `current` (C⋅s⁻¹).
```Julia
julia> abampere(Metric) # C⋅s⁻¹
$(abampere(Metric))

julia> abampere(EMU) # abC⋅s⁻¹
$(abampere(EMU))

julia> abampere(ESU) # statC⋅s⁻¹
$(abampere(ESU))
```
""" abampere

@doc """
$(unitext(:abvolt,"electricpotential(𝟏,U,EMU)"))

Electromagnetic unit of `electricpotential` (V).
```Julia
julia> abvolt(Metric) # V
$(abvolt(Metric))

julia> abvolt(EMU) # abV
$(abvolt(EMU))

julia> abvolt(ESU) # statV
$(abvolt(ESU))
```
""" abvolt

@doc """
$(unitext(:abhenry,"inductance(𝟏,U,EMU)"))

Electromagnetic unit of `inductance` (H).
```Julia
julia> abhenry(Metric) # H
$(abhenry(Metric))

julia> abhenry(EMU) # abH
$(abhenry(EMU))

julia> abhenry(ESU) # statH
$(abhenry(ESU))
```
""" abhenry

@doc """
$(unitext(:abohm,"resistance(𝟏,U,EMU)"))

Electromagnetic unit of `resistance` (Ω).
```Julia
julia> abohm(Metric) # Ω
$(abohm(Metric))

julia> abohm(EMU) # abΩ
$(abohm(EMU))

julia> abohm(ESU) # statΩ
$(abohm(ESU))
```
""" abohm

@doc """
$(unitext(:abmho,"conductance(𝟏,U,EMU)"))

Electromagnetic unit of `conductance` (S).
```Julia
julia> abmho(Metric) # S
$(abmho(Metric))

julia> abmho(EMU) # abS
$(abmho(EMU))

julia> abmho(ESU) # statS
$(abmho(ESU))
```
""" abmho

@doc """
$(unitext(:abfarad,"capacitance(𝟏,U,EMU)"))

Electromagnetic unit of `capacitance` (F).
```Julia
julia> abfarad(Metric) # F
$(abfarad(Metric))

julia> abfarad(EMU) # abF
$(abfarad(EMU))

julia> abfarad(ESU) # statF
$(abfarad(ESU))
```
""" abfarad

@doc """
$(unitext(:maxwell,"magneticflux(𝟏,U,EMU)"))

Electromagnetic unit of `magneticflux` (Wb).
```Julia
julia> maxwell(Metric) # Wb
$(maxwell(Metric))

julia> maxwell(EMU) # Mx
$(maxwell(EMU))

julia> maxwell(ESU) # statWb
$(maxwell(ESU))
```
""" maxwell

@doc """
$(unitext(:gauss,"magneticfluxdensity(𝟏,U,EMU)"))

Electromagnetic unit of `magneticfluxdensity` (T).
```Julia
julia> gauss(Metric) # T
$(gauss(Metric))

julia> gauss(EMU) # G
$(gauss(EMU))

julia> gauss(ESU) # statT
$(gauss(ESU))
```
""" gauss

@doc """
$(unitext(:oersted,"magneticfield(𝟏,U,EMU)"))

Electromagnetic unit of `magneticfield` (Oe).
```Julia
julia> oersted(Metric) # A⋅m⁻¹
$(oersted(Metric))

julia> oersted(EMU) # Oe
$(oersted(EMU))

julia> oersted(ESU) # statA⋅cm⁻¹
$(oersted(ESU))
```
""" oersted

@doc """
$(unitext(:gilbert,"abampere(U)/𝟐/turn(U)"))

Electromagnetic unit of magnetization (Gb).
```Julia
julia> gilbert(Metric) # A⋅rad⁻¹
$(gilbert(Metric))

julia> gilbert(EMU) # Gb
$(gilbert(EMU))

julia> gilbert(ESU) # statA⋅rad⁻¹
$(gilbert(ESU))
```
""" gilbert

@doc """
$(unitext(:statcoulomb,"charge(𝟏,U,ESU)"))

Electrostatic unit of `charge` (C).
```Julia
julia> statcoulomb(Metric) # C
$(statcoulomb(Metric))

julia> statcoulomb(EMU) # abC
$(statcoulomb(EMU))

julia> statcoulomb(ESU) # statC
$(statcoulomb(ESU))
```
""" statcoulomb

@doc """
$(unitext(:statampere,"current(𝟏,U,ESU)"))

Electrostatic unit of `current` (C⋅s⁻¹).
```Julia
julia> statampere(Metric) # C⋅s⁻¹
$(statampere(Metric))

julia> statampere(EMU) # abC⋅s⁻¹
$(statampere(EMU))

julia> statampere(ESU) # statC⋅s⁻¹
$(statampere(ESU))
```
""" statampere

@doc """
$(unitext(:statvolt,"electricpotential(𝟏,U,ESU)"))

Electrostatic unit of `electricpotential` (V).
```Julia
julia> statvolt(Metric) # V
$(statvolt(Metric))

julia> statvolt(EMU) # abV
$(statvolt(EMU))

julia> statvolt(ESU) # statV
$(statvolt(ESU))
```
""" statvolt

@doc """
$(unitext(:stathenry,"inductance(𝟏,U,ESU)"))

Electrostatic unit of `inductance` (H).
```Julia
julia> stathenry(Metric) # H
$(stathenry(Metric))

julia> stathenry(EMU) # abH
$(stathenry(EMU))

julia> stathenry(ESU) # statH
$(stathenry(ESU))
```
""" stathenry

@doc """
$(unitext(:statohm,"resistance(𝟏,U,ESU)"))

Electrostatic unit of `resistance` (Ω).
```Julia
julia> statohm(Metric) # Ω
$(statohm(Metric))

julia> statohm(EMU) # abΩ
$(statohm(EMU))

julia> statohm(ESU) # statΩ
$(statohm(ESU))
```
""" statohm

@doc """
$(unitext(:statmho,"conductance(𝟏,U,ESU)"))

Electrostatic unit of `conductance` (S).
```Julia
julia> statmho(Metric) # S
$(statmho(Metric))

julia> statmho(EMU) # abS
$(statmho(EMU))

julia> statmho(ESU) # statS
$(statmho(ESU))
```
""" statmho

@doc """
$(unitext(:statfarad,"capacitance(𝟏,U,ESU)"))

Electrostatic unit of `capacitance` (F).
```Julia
julia> statfarad(Metric) # F
$(statfarad(Metric))

julia> statfarad(EMU) # abF
$(statfarad(EMU))

julia> statfarad(ESU) # statF
$(statfarad(ESU))
```
""" statfarad

@doc """
$(unitext(:statweber,"magneticflux(𝟏,U,ESU)"))

Electrostatic unit of `magneticflux` (Wb).
```Julia
julia> statweber(Metric) # Wb
$(statweber(Metric))

julia> statweber(EMU) # Mx
$(statweber(EMU))

julia> statweber(ESU) # statWb
$(statweber(ESU))
```
""" statweber

@doc """
$(unitext(:stattesla,"magneticfluxdensity(𝟏,U,ESU)"))

Electrostatic unit of `magneticfluxdensity` (T).
```Julia
julia> stattesla(Metric) # T
$(stattesla(Metric))

julia> stattesla(EMU) # G
$(stattesla(EMU))

julia> stattesla(ESU) # statT
$(stattesla(ESU))
```
""" stattesla

@doc """
$(unitext(:earthcoulomb,"charge(𝟏,U,Meridian)"))

Meridian unit of `charge` (C).
```Julia
julia> earthcoulomb(Metric) # C
$(earthcoulomb(Metric))

julia> earthcoulomb(EMU) # abC
$(earthcoulomb(EMU))

julia> earthcoulomb(ESU) # statC
$(earthcoulomb(ESU))
```
""" earthcoulomb

# temperature

@doc """
$(unitext(:kelvin,"temperature(𝟏,U,Metric)"))

Metric unit of `temperature` (K or °R).
```Julia
julia> kelvin(Metric) # K
$(kelvin(Metric))

julia> kelvin(SI2019) # K
$(kelvin(SI2019))

julia> kelvin(British) # °R
$(kelvin(British))
```
""" kelvin

@doc """
$(unitext(:rankine,"temperature(𝟏,U,English)"))

English unit of `temperature` (K or °R).
```Julia
julia> rankine(Metric) # K
$(rankine(Metric))

julia> rankine(SI2019) # K
$(rankine(SI2019))

julia> rankine(British) # °R
$(rankine(British))
```
""" rankine

@doc """
$(unitext(:celsius,"temperature(T₀,U,Metric)"))

Metric unit of `temperature` (K or °R).
```Julia
julia> celsius(Metric) # K
$(celsius(Metric))

julia> celsius(SI2019) # K
$(celsius(SI2019))

julia> celsius(British) # °R
$(celsius(British))
```
""" celsius, T₀

@doc """
$(unitext(:fahrenheit,"temperature(Constant(459.67),U,English)"))

English unit of `temperature` (K or °R).
```Julia
julia> fahrenheit(Metric) # K
$(fahrenheit(Metric))

julia> fahrenheit(SI2019) # K
$(fahrenheit(SI2019))

julia> fahrenheit(British) # °R
$(fahrenheit(British))
```
""" fahrenheit

#=@doc """
$(unitext(:delisle,"temperature(𝟐/𝟑,U,Metric)"))

Historical unit of `temperature` (K or °R).
```Julia
julia> delisle(Metric) # K
$(delisle(Metric))

julia> delisle(SI2019) # K
$(delisle(SI2019))

julia> delisle(British) # °R
$(delisle(British))
```
""" delisle

@doc """
$(unitext(:reaumur,"temperature(𝟓/𝟐^2,U,Metric)"))

Historical unit of `temperature` (K or °R).
```Julia
julia> reaumur(Metric) # K
$(reaumur(Metric))

julia> reaumur(SI2019) # K
$(reaumur(SI2019))

julia> reaumur(British) # °R
$(reaumur(British))
```
""" reaumur

@doc """
$(unitext(:freezing,"temperature(T₀-milli,U)"))

Standard `temperature` reference at `freezing` point of water (K or °R).
```Julia
julia> freezing(Metric) # K
$(freezing(Metric))

julia> freezing(SI2019) # K
$(freezing(SI2019))

julia> freezing(English) # °R
$(freezing(English))
```
""" freezing=#

@doc """
$(unitext(:boiling,"temperature(T₀+Constant(99.9839),U)"))

Standard `temperature` reference at `boiling` point of water (K or °R).
```Julia
julia> boiling(Metric) # K
$(boiling(Metric))

julia> boiling(SI2019) # K
$(boiling(SI2019))

julia> boiling(English) # °R
$(boiling(English))
```
""" boiling

@doc """
$(unitext(:sealevel,"temperature(T₀+𝟑*𝟓,U)"))

Standard `temperature` reference at `sealevel` (K or °R).
```Julia
julia> sealevel(Metric) # K
$(sealevel(Metric))

julia> sealevel(SI2019) # K
$(sealevel(SI2019))

julia> sealevel(English) # °R
$(sealevel(English))
```
""" sealevel

# mole

@doc """
$(unitext(:mole,"molaramount(𝟏,U,Metric)"))

Molecular `molaramount` unit (mol or lb-mol).
```Julia
julia> mole(Metric) # mol
$(mole(Metric))

julia> mole(English) # lb-mol
$(mole(English))

julia> mole(British) # slug-mol
$(mole(British))
```
""" mole

@doc """
$(unitext(:earthmole,"molaramount(𝟏,U,Meridian)"))

Molecular `molaramount` unit (mol or lb-mol).
```Julia
julia> earthmole(Metric) # mol
$(earthmole(Metric))

julia> earthmole(English) # lb-mol
$(earthmole(English))

julia> earthmole(British) # slug-mol
$(earthmole(British))
```
""" earthmole

@doc """
$(unitext(:poundmole,"molaramount(𝟏,U,English)"))

Molecular `molaramount` unit (mol or lb-mol).
```Julia
julia> poundmole(Metric) # mol
$(poundmole(Metric))

julia> poundmole(English) # lb-mol
$(poundmole(English))

julia> poundmole(British) # slug-mol
$(poundmole(British))
```
""" poundmole

@doc """
$(unitext(:slinchmole,"molaramount(𝟏,U,IPS)"))

Molecular `molaramount` unit (mol or lb-mol).
```Julia
julia> slinchmole(Metric) # mol
$(slinchmole(Metric))

julia> slinchmole(English) # lb-mol
$(slinchmole(English))

julia> slinchmole(British) # slug-mol
$(slinchmole(British))
```
""" slinchmole

@doc """
$(unitext(:slugmole,"molaramount(𝟏,U,British)"))

Molecular `molaramount` unit (mol or lb-mol).
```Julia
julia> slugmole(Metric) # mol
$(slugmole(Metric))

julia> slugmole(English) # lb-mol
$(slugmole(English))

julia> slugmole(British) # slug-mol
$(slugmole(British))
```
""" slugmole

# photometric

@doc """
$(unitext(:lumen,"luminousflux(𝟏,U,Metric)"))

Common unit of `luminousflux` (lm).
```Julia
julia> lumen(Metric) # lm
$(lumen(Metric))

julia> lumen(CGS) # lm
$(lumen(CGS))

julia> lumen(English) # lm
$(lumen(English))
```
""" lumen

@doc """
$(unitext(:candela,"luminousintensity(𝟏,U,Metric)"))

Common unit of `luminousintensity` (cd).
```Julia
julia> candela(MetricEngineering) # lm⋅rad⁻²
$(candela(MetricEngineering))

julia> candela(MetricDegree) # lm⋅deg⁻²
$(candela(MetricDegree))

julia> candela(MetricGradian) # lm⋅gon⁻²
$(candela(MetricGradian))

julia> candela(CGS) # cd
$(candela(CGS))

julia> candela(English) # cd
$(candela(English))
```
""" candela

@doc """
$(unitext(:lux,"illuminance(𝟏,U,Metric)"))

Metric unit of `illuminance` (lx).
```Julia
julia> lux(Metric) # lx
$(lux(Metric))

julia> lux(CGS) # ph
$(lux(CGS))

julia> lux(English) # fc
$(lux(English))
```
""" lux

@doc """
$(unitext(:footcandle,"illuminance(𝟏,U,English)"))

English unit of `illuminance` (lx).
```Julia
julia> footcandle(Metric) # lx
$(footcandle(Metric))

julia> footcandle(CGS) # ph
$(footcandle(CGS))

julia> footcandle(English) # fc
$(footcandle(English))
```
""" footcandle

@doc """
$(unitext(:phot,"illuminance(𝟏,U,Gauss)"))

Historic unit of `illuminance` (lx).
```Julia
julia> phot(Metric) # lx
$(phot(Metric))

julia> phot(CGS) # ph
$(phot(CGS))

julia> phot(English) # fc
$(phot(English))
```
""" phot

@doc """
$(unitext(:nit,"luminance(𝟏,U,Metric)"))

Metric unit of `luminance` (lx⋅rad⁻²).
```Julia
julia> nit(MetricEngineering) # nt
$(nit(MetricEngineering))

julia> nit(MetricDegree) # lm⋅m⁻²deg⁻²
$(nit(MetricDegree))

julia> nit(MetricGradian) # lm⋅m⁻²gon⁻²
$(nit(MetricGradian))

julia> nit(CGS) # sb
$(nit(CGS))

julia> nit(English) # fc
$(nit(English))
```
""" nit

@doc """
$(unitext(:apostilb,"luminance(𝟐/turn(U),U,Metric)"))

Metric unit of `luminance` (lx⋅rad⁻²).
```Julia
julia> apostilb(MetricEngineering) # nt
$(apostilb(MetricEngineering))

julia> apostilb(MetricDegree) # lm⋅m⁻²deg⁻²
$(apostilb(MetricDegree))

julia> apostilb(MetricGradian) # lm⋅m⁻²gon⁻²
$(apostilb(MetricGradian))

julia> apostilb(CGS) # sb
$(apostilb(CGS))

julia> apostilb(English) # fc
$(apostilb(English))
```
""" apostilb

@doc """
$(unitext(:stilb,"luminance(𝟏,U,Gauss)"))

Historic unit of `luminance` (lx⋅rad⁻²).
```Julia
julia> stilb(MetricEngineering) # nt
$(stilb(MetricEngineering))

julia> stilb(MetricDegree) # lm⋅m⁻²deg⁻²
$(stilb(MetricDegree))

julia> stilb(MetricGradian) # lm⋅m⁻²gon⁻²
$(stilb(MetricGradian))

julia> stilb(CGS) # sb
$(stilb(CGS))

julia> stilb(English) # fc
$(stilb(English))
```
""" stilb

@doc """
$(unitext(:lambert,"luminance(𝟐/turn(U),U,Gauss)"))

Historic unit of `luminance` (nt).
```Julia
julia> lambert(MetricEngineering) # nt
$(lambert(MetricEngineering))

julia> lambert(MetricDegree) # lm⋅m⁻²deg⁻²
$(lambert(MetricDegree))

julia> lambert(MetricGradian) # lm⋅m⁻²gon⁻²
$(lambert(MetricGradian))

julia> lambert(CGS) # sb
$(lambert(CGS))

julia> lambert(English) # fc
$(lambert(English))
```
""" lambert

@doc """
$(unitext(:footlambert,"luminance(𝟐/turn(U),U,English)"))

English unit of `luminance` (nt).
```Julia
julia> footlambert(MetricEngineering) # nt
$(footlambert(MetricEngineering))

julia> footlambert(MetricDegree) # lm⋅m⁻²deg⁻²
$(footlambert(MetricDegree))

julia> footlambert(MetricGradian) # lm⋅m⁻²gon⁻²
$(footlambert(MetricGradian))

julia> footlambert(CGS) # sb
$(footlambert(CGS))

julia> footlambert(English) # fc
$(footlambert(English))
```
""" footlambert

@doc """
$(unitext(:bril,"centi*nano*lambert(U)"))

Reference unit of `luminance` (nt).
```Julia
julia> bril(MetricEngineering) # nt
$(bril(MetricEngineering))

julia> bril(MetricDegree) # lm⋅m⁻²deg⁻²
$(bril(MetricDegree))

julia> bril(MetricGradian) # lm⋅m⁻²gon⁻²
$(bril(MetricGradian))

julia> bril(CGS) # sb
$(bril(CGS))

julia> bril(English) # fc
$(bril(English))
```
""" bril

# special

@doc """
$(unitext(:hertz,"𝟏/second(U)"))

Metric unit of `frequency` (s⁻¹).
```Julia
julia> hertz(MetricEngineering) # rad⋅s⁻¹
$(hertz(MetricEngineering))

julia> hertz(IAU) # D⁻¹
$(hertz(IAU))
```
""" hertz

@doc """
$(unitext(:apm,"𝟏/minute(U)"))

Actions per minute `apm` unit of `frequency` (s⁻¹).
```Julia
julia> apm(Metric) # s⁻¹
$(apm(Metric))

julia> apm(MPH) # h⁻¹
$(apm(MPH))

julia> apm(IAU) # D⁻¹
$(apm(IAU))
```
""" apm

@doc """
$(unitext(:rpm,"turn(U)/minute(U)"))

Revolutions per minute `rpm` unit of `angularfrequency` (rad⋅s⁻¹).
```Julia
julia> rpm(MetricEngineering) # rad⋅s⁻¹
$(rpm(MetricEngineering))

julia> rpm(MetricGradian) # gon⋅s⁻¹
$(rpm(MetricGradian))

julia> rpm(MetricDegree) # deg⋅s⁻¹
$(rpm(MetricDegree))

julia> rpm(MetricArcminute) # amin⋅s⁻¹
$(rpm(MetricArcminute))

julia> rpm(MetricArcsecond) # asec⋅s⁻¹
$(rpm(MetricArcsecond))

julia> rpm(MPH) # rad⋅h⁻¹
$(rpm(MPH))

julia> rpm(IAU) # rad⋅D⁻¹
$(rpm(IAU))
```
""" rpm

@doc """
$(unitext(:kayser,"wavenumber(𝟏,U,Gauss)"))

Metric unit of `wavenumber` or curvature (m⁻¹ or ft⁻¹).
```Julia
julia> kayser(Metric) # m⁻¹
$(kayser(Metric))

julia> kayser(CGS) # cm⁻¹
$(kayser(CGS))

julia> kayser(English) # ft⁻¹
$(kayser(English))
```
""" kayser

@doc """
$(unitext(:diopter,"wavenumber(𝟏,U,Metric)"))

Metric unit of `wavenumber` or curvature (m⁻¹ or ft⁻¹).
```Julia
julia> diopter(Metric) # m⁻¹
$(diopter(Metric))

julia> diopter(CGS) # cm⁻¹
$(diopter(CGS))

julia> diopter(English) # ft⁻¹
$(diopter(English))
```
""" diopter

@doc """
$(unitext(:gforce,"specificforce(𝟏,U,English)"))

Standard gravity `specificforce` `g₀` at geodetic reference latitude (m⋅s⁻² or ft⋅s⁻²).
```Julia
julia> gforce(CGS) # gal
$(gforce(CGS))

julia> gforce(British) # ft⋅s⁻²
$(gforce(British))

julia> gforce(English) # lbf⋅lbm⁻¹
$(gforce(English))
```
""" gforce, g₀, g0, lbm

@doc """
$(unitext(:galileo,"specificforce(𝟏,U,Gauss)"))

Metric unit of `specificforce` used in gravimetry (m⋅s⁻² or ft⋅s⁻²).
```Julia
julia> galileo(Metric) # m⋅s⁻²
$(galileo(Metric))

julia> galileo(CGS) # gal
$(galileo(CGS))

julia> galileo(English) # lbf⋅lbm⁻¹
$(galileo(English))
```
""" galileo

@doc """
$(unitext(:eotvos,"specificforce(nano,U,Gauss)/length(𝟏,U,Gauss)"))

Metric unit of `specificforce` per `length` used in gravimetry (s⁻² or gal⋅cm⁻¹).
```Julia
julia> eotvos(Metric) # s⁻²
$(eotvos(Metric))

julia> eotvos(CGS) # gal⋅cm⁻¹
$(eotvos(CGS))

julia> eotvos(English) # lbf⋅lbm⁻¹ft⁻¹
$(eotvos(English))
```
""" eotvos

@doc """
$(unitext(:darcy,"area(milli/atm,U,Gauss)"))

Metric unit of permeability (m² or ft²).
```Julia
julia> darcy(Metric) # m²
$(darcy(Metric))

julia> darcy(CGS) # cm²
$(darcy(CGS))

julia> darcy(English) # ft²
$(darcy(English))
```
""" darcy

@doc """
$(unitext(:poise,"viscosity(𝟏,U,Gauss)"))

Metric unit of `viscosity` (kg⋅m⁻¹⋅s⁻¹ or lb⋅s⋅ft⁻²).
```Julia
julia> poise(Metric) # kg⋅m⁻¹⋅s⁻¹
$(poise(Metric))

julia> poise(CGS) # g⋅cm⁻¹⋅s⁻¹
$(poise(CGS))

julia> poise(English) # lb⋅s⋅ft⁻²
$(poise(English))
```
""" poise

@doc """
$(unitext(:reyn,"viscosity(𝟏,U,IPS)"))

IPS unit of `viscosity` named after Reynolds (kg⋅m⁻¹⋅s⁻¹ or lb⋅s⋅ft⁻²).
```Julia
julia> reyn(Metric) # kg⋅m⁻¹⋅s⁻¹
$(reyn(Metric))

julia> reyn(CGS) # g⋅cm⁻¹⋅s⁻¹
$(reyn(CGS))

julia> reyn(English) # lb⋅s⋅ft⁻²
$(reyn(English))
```
""" reyn

@doc """
$(unitext(:stokes,"diffusivity(𝟏,U,Gauss)"))

Metric unit of `diffusivity` (m²⋅s⁻¹ or ft²⋅s⁻¹).
```Julia
julia> stokes(Metric) # m²⋅s⁻¹
$(stokes(Metric))

julia> stokes(CGS) # cm²⋅s⁻¹
$(stokes(CGS))

julia> stokes(English) # ft²⋅s⁻¹
$(stokes(English))
```
""" stokes

@doc """
$(unitext(:rayl,"specificimpedance(𝟏,U,Metric)"))

Metric unit of `specificimpedance` (kg⋅m⁻²⋅s⁻¹ or lb⋅s⋅ft⁻³).
```Julia
julia> rayl(Metric) # kg⋅m⁻²⋅s⁻¹
$(rayl(Metric))

julia> rayl(CGS) # g⋅cm⁻²⋅s⁻¹
$(rayl(CGS))

julia> rayl(English) # lb⋅s⋅ft⁻³
$(rayl(English))
```
""" rayl

@doc """
$(unitext(:katal,"catalysis(𝟏,U,Metric)"))

Metric unit of `catalysis` (mol⋅s⁻¹ or lb-mol⋅s⁻¹).
```Julia
julia> katal(Metric) # mol⋅s⁻¹
$(katal(Metric))

julia> katal(English) # lb-mol⋅s⁻¹
$(katal(English))

julia> katal(British) # slug-mol⋅s⁻¹
$(katal(British))
```
""" katal

@doc """
$(unitext(:mpge,"mile(U)/gasgallon(U)"))

Equivalent `mile` per `gasgallon` reference unit of `length` per `energy` (N⁻¹ or lb⁻¹).
```Julia
julia> mpge(Metric) # N⁻¹
$(mpge(Metric))

julia> mpge(CGS) # dyn⁻¹
$(mpge(CGS))

julia> mpge(English) # lb⁻¹
$(mpge(English))
```
""" mpge

@doc """
$(unitext(:langley,"calorie(U)/(centi*meter(U))^2"))

Solar radiation unit (kg⋅s⁻² or lb⋅ft⁻¹).
```Julia
julia> langley(Metric) # kg⋅s⁻²
$(langley(Metric))

julia> langley(CGS) # g⋅s⁻²
$(langley(CGS))

julia> langley(English) # lb⋅ft⁻¹
$(langley(English))
```
""" langley

@doc """
$(unitext(:jansky,"fluence(deci^-26,U,Metric)"))

Reference unit of spectral irradiance (kg⋅s⁻² or lb⋅ft⁻¹).
```Julia
julia> jansky(Metric) # kg⋅s⁻²
$(jansky(Metric))

julia> jansky(CGS) # g⋅s⁻²
$(jansky(CGS))

julia> jansky(English) # lb⋅ft⁻¹
$(jansky(English))
```
""" jansky

@doc """
$(unitext(:solarflux,"hecto^2*jansky(U)"))

Reference unit of spectral irradiance (kg⋅s⁻² or lb⋅ft⁻¹).
```Julia
julia> solarflux(Metric) # kg⋅s⁻²
$(solarflux(Metric))

julia> solarflux(CGS) # g⋅s⁻²
$(solarflux(CGS))

julia> solarflux(English) # lb⋅ft⁻¹
$(solarflux(English))
```
""" solarflux

@doc """
$(unitext(:curie,"frequency(𝟏,U,Metric)"))

Reference unit of radioactivity (Bq or s⁻¹).
```Julia
julia> curie(Metric) # Bq
$(curie(Metric))

julia> curie(IAU) # D⁻¹
$(curie(IAU))
```
""" curie

@doc """
$(unitext(:sievert,"energy(𝟏,U,Metric)/mass(U,Metric)"))

Metric unit of radioactivity (Sv or m²⋅s⁻²).
```Julia
julia> sievert(Metric) # Sv
$(sievert(Metric))
```
""" sievert

@doc """
    rem(U::UnitSystem) = centi*sievert(U)

Obsolete unit of radioactivity (Sv or m²⋅s⁻²).
```Julia
julia> rem(Metric) # Sv
$(rem(Metric))
```
""" rem

@doc """
$(unitext(:roentgen,"chargedensity(𝟏,U,ESU)/density(Constant(1.293),U,Metric)"))

Legacy unit of ionisation `exposure` (R or C⋅kg⁻¹).
```Julia
julia> roentgen(Metric) # R
$(roentgen(Metric))
```
""" roentgen
