using UnitSystems, Test

@test molarmass(Natural) == molarmass(CGS) == 1000molarmass(Metric)
@test molarmass(CGS2019) == 1000molarmass(SI2019)
