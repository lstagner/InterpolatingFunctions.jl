immutable Polynomial <: AbstractInterpolation
end
Polynomial() = Polynomial

function default_settings{D,T}(::Type{Polynomial},axes::NTuple{D,T})
    N = length(axes[1])
    return Irregular{1,N}, None
end

type PolynomialInterpolation{T,N}
  coeff::Vector{T}
  x::Vector{T}
  y::Vector{T}
end

function preprocess(::Type{Polynomial}, axes, f)
    x = convert(Array, axes[1])
    y = convert(Array, f)
    nx = length(x)
    coeff = ones(nx)

    for i=1:nx
        for j = 1:nx
            i == j && continue
            coeff[i] = coeff[i]*(x[i] - x[j])
        end
        coeff[i] = (y[i]/coeff[i])
    end
    x,y,coeff = promote(x,y,coeff)
    return PolynomialInterpolation{eltype(x),nx}(coeff,x,y)
end

function apply_bc(::Type{Polynomial}, ::Type{None}, I::Symbol, x::Symbol)
    :()
end

@inline function get_indices{D,N}(::Type{Polynomial}, G::Type{Irregular{D,N}}, I, x)
    return ()
end

@inline function get_value{T,N}(::Type{Polynomial}, I::PolynomialInterpolation{T,N}, x, inds)
    coeff = I.coeff
    xx = I.x
    y = 0.0
    @inbounds for i=1:N
        tmp = 1.0
        for j = 1:N
            j == i && continue
            tmp = tmp * (x - xx[j])
        end
        y = y + tmp*coeff[i]
    end
    return y
end
