immutable LinearSpline <: AbstractInterpolation
    LinearSpline() = LinearSpline
end

immutable LinearSplineInterpolation{T,N}
    x::Vector{T}
    y::Vector{T}
    m::Vector{T}
end

function default_settings{D,T}(::Type{LinearSpline},axes::NTuple{D,T})
    N = length(axes[1])
    return Irregular{D,N}, Nearest
end

function LinearSpline(x::AbstractVector,y::AbstractVector)
  issorted(x) || throw(ArgumentError("x points must be in ascending order"))
  nx = length(x)
  x = convert(Array,x)
  y = convert(Array,y)
  m = zeros(nx)

  for i=1:nx-1
      m[i] = (y[i+1]-y[i])/(x[i+1] - x[i])
  end
  x, y, m = promote(x,y,m)
  return LinearSplineInterpolation{eltype(x),nx}(x,y,m)
end

@inline function preprocess(::Type{LinearSpline}, axes, f)
    return LinearSpline(axes...,f)
end

function apply_bc(::Type{LinearSpline},BC::Type{Nearest}, I::Symbol, x::Symbol)
    apply_nearest(I,1,x,1)
end

function apply_bc(::Type{LinearSpline},BC::Type{Reflect}, I::Symbol, x::Symbol)
    apply_reflect(I,1,x,1)
end

function apply_bc{V}(::Type{LinearSpline},BC::Type{Value{V}}, I::Symbol, x::Symbol)
    apply_value(BC, I, 1, x, 1)
end

function apply_bc(::Type{LinearSpline},BC::Type{None}, I::Symbol, x::Symbol)
    apply_none()
end

function apply_bc(::Type{LinearSpline},BC::Type{Error}, I::Symbol, x::Symbol)
    apply_error(I, 1, x, 1)
end

@inline function get_indices{S}(::Type{LinearSpline}, G::Type{Uniform{1,S}}, I, x)
    dx = I.x[2]-I.x[1]
    i1 = clamp(Int(floor((x[1] - I.x[1])/dx)) + 1, 1, S.parameters[1][1]-1)
    return i1
end

@inline function get_indices{N}(::Type{LinearSpline}, G::Type{Irregular{1,N}}, I, x)
    xr = searchsortedlast(I.x,x[1],Base.Order.Forward)
    i1 = clamp(xr, 1, N-1)
    return i1
end

@inline function get_value(::Type{LinearSpline}, I, x, ix)
    return I.y[ix] + I.m[ix]*(x[1] - I.x[ix])
end
