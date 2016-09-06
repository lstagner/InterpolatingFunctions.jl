type Monotonic end
type NonMonotonic end

immutable CubicSpline{M} <: AbstractInterpolation end
CubicSpline() = CubicSpline{NonMonotonic}

immutable CubicSplineInterpolation{T,N,M}
    x::Vector{T}
    y::Vector{T}
    m::Vector{T}
end

function default_settings{D,T,M}(::Type{CubicSpline{M}},axes::NTuple{D,T})
    N = length(axes[1])
    return Irregular{D,N}, Nearest
end

function CubicSpline(x::AbstractVector,y::AbstractArray,M::Type{NonMonotonic})
  issorted(x) || throw(ArgumentError("x points must be in ascending order"))
  nx = length(x)
  x = convert(Array,x)
  y = convert(Array,y)
  m = zeros(nx)
  m[1] = (y[2] - y[1])/(x[2]-x[1])
  for i=2:nx-1
    m[i] = 0.5*((y[i+1]-y[i])/(x[i+1]-x[i]) + (y[i]-y[i-1])/(x[i]-x[i-1]))
  end
  m[nx] = (y[nx]-y[nx-1])/(x[nx]-x[nx-1])
  x,y,m = promote(x,y,m)
  return CubicSplineInterpolation{eltype(x),length(x),M}(x,y,m)
end

function CubicSpline(x::AbstractVector,y::AbstractVector,M::Type{Monotonic})
  issorted(x) || throw(ArgumentError("x points must be in ascending order"))
  nx = length(x)
  m = zeros(nx)
  x = convert(Array,x)
  y = convert(Array,y)

  m[1] = (y[2] - y[1])/(x[2]-x[1])
  for i=2:nx-1
    hi = x[i+1]-x[i]
    hi_1 = x[i]-x[i-1]
    di = (y[i+1]-y[i])/hi
    di_1 = (y[i]-y[i-1])/hi_1
    m[i] = sign(di) != sign(di_1) ? 0.0 : 3*(hi_1+hi)*(((2*hi+hi_1)/di_1)+((hi+2*hi_1)/di))^(-1.0)
  end
  m[nx] = (y[nx]-y[nx-1])/(x[nx]-x[nx-1])

  x,y,m = promote(x,y,m)
  return CubicSplineInterpolation{eltype(x),length(x),M}(x,y,m)
end

@inline function preprocess{M}(::Type{CubicSpline{M}}, axes, f)
    return CubicSpline(axes...,f, M)
end

function apply_bc{M}(::Type{CubicSpline{M}},BC::Type{Nearest}, I::Symbol, x::Symbol)
    apply_nearest(I,1,x,1)
end

function apply_bc{M}(::Type{CubicSpline{M}},BC::Type{Reflect}, I::Symbol, x::Symbol)
    apply_reflect(I,1,x,1)
end

function apply_bc{M,V}(::Type{CubicSpline{M}},BC::Type{Value{V}}, I::Symbol, x::Symbol)
    apply_value(BC,I,1,x,1)
end

function apply_bc{M}(::Type{CubicSpline{M}},BC::Type{None}, I::Symbol, x::Symbol)
    apply_none()
end

function apply_bc{M}(::Type{CubicSpline{M}},BC::Type{Error}, I::Symbol, x::Symbol)
    apply_error(I,1,x,1)
end

@inline function get_indices{M,S}(::Type{CubicSpline{M}}, G::Type{Uniform{1,S}}, I, x)
    dx = I.x[2]-I.x[1]
    i1 = clamp(Int(floor((x - I.x[1])/dx)) + 1, 1, S.parameters[1][1]-1)
    return i1
end

@inline function get_indices{M,N}(::Type{CubicSpline{M}}, G::Type{Irregular{1,N}}, I, x)
    xr = searchsortedlast(I.x,x[1],Base.Order.Forward)
    i1 = clamp(xr, 1, N-1)
    return i1
end

@inline function get_value{M}(::Type{CubicSpline{M}}, I, x, ind)
    i1 = ind
    i2 = i1 + 1
    dx = (I.x[i2] - I.x[i1])
    t = (x[1] - I.x[i1])/dx
    h00 = 2*t^3 - 3*t^2 + 1
    h10 = t^3 - 2*t^2 + t
    h01 = -2*t^3 + 3*t^2
    h11 = t^3 - t^2
    return h00*I.y[i1] + h10*dx*I.m[i1] + h01*I.y[i2] + h11*dx*I.m[i2]
end
