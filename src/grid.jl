abstract AbstractGrid{D}

type Uniform{D,S} <: AbstractGrid{D} end
Uniform{D,T<:Integer}(x::NTuple{D,T}) = Uniform{D,Val{x}}
Uniform(x::Int) = Uniform{1,Val{tuple(x)}}

type Irregular{D,N} <: AbstractGrid{D} end
Irregular(D::Integer,N::Integer) = Irregular{D,N}
Irregular(N::Integer) = Irregular{1,N}
