
#   This file is part of UnitSystems.jl
#   It is licensed under the MIT license
#   UnitSystems Copyright (C) 2022 Michael Reed
#       _           _                         _
#      | |         | |                       | |
#   ___| |__   __ _| | ___ __ __ ___   ____ _| | __ _
#  / __| '_ \ / _` | |/ / '__/ _` \ \ / / _` | |/ _` |
# | (__| | | | (_| |   <| | | (_| |\ V / (_| | | (_| |
#  \___|_| |_|\__,_|_|\_\_|  \__,_| \_/ \__,_|_|\__,_|
#
#   https://github.com/chakravala
#   https://crucialflow.com

export Constant, constant

# constant

struct Constant{D} <: Real
    @pure Constant{D}() where D = new{D}()
end

@pure param(::Constant{D}) where D = D
Base.float(x::Constant) = float(constant(x))
Base.convert(::Type{Float64},c::Constant) = float(c)
logdb(x::Constant{D}) where D = Constant{logdb(D)}()
Base.:+(a::Constant,b::Constant) = Constant(constant(a)+constant(b))
Base.:-(a::Constant,b::Constant) = Constant(constant(a)-constant(b))
Base.:+(a::Constant{D},::Constant{D}) where D = 𝟐*a
Base.:-(a::Constant{D},::Constant{D}) where D = Constant(0)
Base.:+(a::Number,b::Constant) = a+constant(b)
Base.:+(a::Constant,b::Number) = constant(a)+b
Base.:-(a::Number,b::Constant) = a-constant(b)
Base.:-(a::Constant,b::Number) = constant(a)-b
Base.:*(a::Real,b::Constant) = a*constant(b)
Base.:*(a::Constant,b::Real) = constant(a)*b
Base.:*(a::Constant{A},b::Constant{B}) where {A,B} = Constant{A*B}()
Base.:/(a::Constant{A},b::Constant{B}) where {A,B} = Constant{A/B}()
Base.:/(a::Number,b::Constant) = a*inv(b)
Base.:/(a::Constant,b::Number) = a*inv(b)
Base.inv(a::Constant{D}) where D = Constant{inv(D)}()
Base.sqrt(a::Constant{D}) where D = Constant{sqrt(D)}()
Base.cbrt(a::Constant{D}) where D = Constant{cbrt(D)}()
Base.log(x::Constant{D}) where D = Constant{log(D)}()
Base.log2(x::Constant{D}) where D = Constant{log2(D)}()
Base.log10(x::Constant{D}) where D = Constant{log10(D)}()
Base.log(b::Number,x::Constant{D}) where D = Constant{log(b,D)}()
Base.exp(x::Constant{D}) where D = Constant{exp(D)}()
Base.exp2(x::Constant{D}) where D = Constant{exp2(D)}()
Base.exp10(x::Constant{D}) where D = Constant{exp10(D)}()
Base.:^(a::Number,b::Constant{D}) where D = Constant{a^D}()
Base.:^(a::Constant{D},b::Number) where D = Constant{D^b}()
Base.:^(a::Constant{D},b::Integer) where D = Constant{D^b}()
Base.:^(a::Constant{D},b::Rational{Int}) where D = Constant{D^b}()

Base.isone(x::Constant{D}) where D = isone(D)
Base.iszero(x::Constant{D}) where D = iszero(D)

Base.:(==)(a::Real,b::Constant) = a == constant(b)
Base.:(==)(a::Constant,b::Real) = constant(a) == b
Base.:(==)(a::Constant,b::Constant) = constant(a) == constant(b)
Base.isapprox(a::Real,b::Constant) = isapprox(a,constant(b))
Base.isapprox(a::Constant,b::Real) = isapprox(constant(a),b)
Base.isapprox(a::Constant,b::Constant) = isapprox(constant(a),constant(b))
