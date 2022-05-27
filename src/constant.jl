
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

struct Constant{N} <: Real
    @pure Constant{N}() where N = new{N}()
end

@pure isconstant(::Constant) = true
@pure Constant(N::Constant) = N
@pure Constant(N::Irrational) = Constant(float(N))

@pure param(::Constant{N}) where N = N
Base.abs(::Constant{N}) where N = Constant{abs(N)}()
Base.float(c::Constant) = float(constant(c))
Base.convert(::Type{Float64},c::Constant) = float(c)
logdb(::Constant{N}) where N = Constant{logdb(N)}()
Base.:<(a::Real,b::Constant) = a<constant(b)
Base.:<(a::Constant,b::Real) = constant(a)<b
Base.:+(a::Constant,b::Constant) = Constant(constant(a)+constant(b))
Base.:-(a::Constant,b::Constant) = Constant(constant(a)-constant(b))
Base.:+(a::Constant{N},::Constant{N}) where N = 𝟐*a
Base.:-(a::Constant{N},::Constant{N}) where N = Constant(0)
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
Base.inv(a::Constant{N}) where N = Constant{inv(N)}()
Base.sqrt(a::Constant{N}) where N = Constant{sqrt(N)}()
Base.cbrt(a::Constant{N}) where N = Constant{cbrt(N)}()
Base.log(x::Constant{N}) where N = Constant{log(N)}()
Base.log2(x::Constant{N}) where N = Constant{log2(N)}()
Base.log10(x::Constant{N}) where N = Constant{log10(N)}()
Base.log(b::Number,x::Constant{N}) where N = Constant{log(b,N)}()
Base.exp(x::Constant{N}) where N = Constant{exp(N)}()
Base.exp2(x::Constant{N}) where N = Constant{exp2(N)}()
Base.exp10(x::Constant{N}) where N = Constant{exp10(N)}()
Base.:^(a::Number,b::Constant{N}) where N = Constant{a^N}()
Base.:^(a::Constant{N},b::Number) where N = Constant{N^b}()
Base.:^(a::Constant{N},b::Integer) where N = Constant{N^b}()
Base.:^(a::Constant{N},b::Rational{Int}) where N = Constant{N^b}()

Base.:^(a::Constant,b::Irrational) = a^Constant(b)
Base.:^(a::Irrational,b::Constant) = Constant(a)^b
Base.:*(a::Irrational,b::Constant) = Constant(a)*b
Base.:*(a::Constant,b::Irrational) = a*Constant(b)
Base.:/(a::Irrational,b::Constant) = Constant(a)/b
Base.:/(a::Constant,b::Irrational) = a/Constant(b)
Base.:+(a::Irrational,b::Constant) = Constant(a)+b
Base.:+(a::Constant,b::Irrational) = a+Constant(b)
Base.:-(a::Irrational,b::Constant) = Constant(a)-b
Base.:-(a::Constant,b::Irrational) = a-Constant(b)

Base.isone(x::Constant{N}) where N = isone(N)
Base.iszero(x::Constant{N}) where N = iszero(N)

Base.:(==)(a::Real,b::Constant) = a == constant(b)
Base.:(==)(a::Constant,b::Real) = constant(a) == b
Base.:(==)(a::Constant,b::Constant) = constant(a) == constant(b)
Base.isapprox(a::Real,b::Constant) = isapprox(a,constant(b))
Base.isapprox(a::Constant,b::Real) = isapprox(constant(a),b)
Base.isapprox(a::Constant,b::Constant) = isapprox(constant(a),constant(b))
