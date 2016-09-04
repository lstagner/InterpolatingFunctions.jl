abstract AbstractGrid{D}

abstract AbstractBoundaryCondition

type Uniform{D,S} <: AbstractGrid{D} end
Uniform{D,T<:Integer}(x::NTuple{D,T}) = Uniform{D,Val{x}}
Uniform(x::Int) = Uniform{1,Val{tuple(x)}}

type Irregular{D,N} <: AbstractGrid{D} end
Irregular(D::Integer,N::Integer) = Irregular{D,N}
Irregular(N::Integer) = Irregular{1,N}

type Value{V} <: AbstractBoundaryCondition end
Value(V) = Value{Val{V}}
type Nearest <: AbstractBoundaryCondition end
type Reflect <: AbstractBoundaryCondition end
type None <: AbstractBoundaryCondition end
type Error <: AbstractBoundaryCondition end

type BoundaryConditions{D,C} end
BoundaryConditions(x::Tuple) = BoundaryConditions{length(x),Val{x}}
BoundaryConditions(x) = BoundaryConditions{1,Val{tuple(x)}}

abstract AbstractInterpolation
typealias AI AbstractInterpolation
typealias AG AbstractGrid
typealias ABC AbstractBoundaryCondition

function interpolate{N,T,I<:AI}(f::AbstractArray{T}, axes::NTuple{N,AbstractArray}, ::Type{I})
    G,BC = default_settings(I,axes)
    interpolate(f,axes,I,G,BC)
end

@generated function interpolate{N,T,I<:AI,G<:AG,BC<:ABC}(f::AbstractArray{T}, axes::NTuple{N,AbstractVector},::Type{I},::Type{G},::Type{BC})
    bc_expr = apply_bc(I, BC, :interp, :x)
    quote
        interp = preprocess($I, axes, f)
        function InterpolatingFunction(x)
            $(bc_expr)

            inds = get_indices($I, $G, interp, x)
            val = get_value($I, interp, x, inds)
            return val
        end
    end
end
