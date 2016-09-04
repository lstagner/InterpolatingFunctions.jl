immutable LinearSpline{D} <: AbstractInterpolation end
LinearSpline() = LinearSpline{1}
LinearSpline(D) = LinearSpline{D}

immutable LinearSplineInterpolation{T,D,N}
    x::Vector{T}
    y::Vector{T}
    m::Vector{T}
end

function default_settings{D,T}(::Type{LinearSpline{D}},axes::NTuple{D,T})
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
  return LinearSplineInterpolation{eltype(x),1,nx}(x,y,m)
end

@inline function preprocess{D}(::Type{LinearSpline{D}}, axes, f)
    return LinearSpline(axes...,f)
end

function apply_bc{D}(::Type{LinearSpline{D}},BC::Type{Nearest}, I::Symbol, x::Symbol)
    :($x = clamp($x,$I.x[1],$I.x[end]))
end

function apply_bc{D}(::Type{LinearSpline{D}},BC::Type{Reflect}, I::Symbol, x::Symbol)
    quote
        if $x < $I.x[1]
            delX = $I.x[end] - $I.x[1]
            n = floor(($I.x[1] - $x)/delX)
            $x = $I.x[end] - (($I.x[1] - $x) - n*delX)
        elseif $x > $I.x[end]
            delX = $I.x[end] - $I.x[1]
            n = floor(($x - $I.x[end])/delX)
            x = $x[1] + ($x - $I.x[end] - n*delX)
        end
    end
end

function apply_bc{D,V}(::Type{LinearSpline{D}},BC::Type{Value{V}}, I::Symbol, x::Symbol)
    quote
        if $x < $I.x[1] || $x > $I.x[end]
            return $V.parameters[1]
        end
    end
end

function apply_bc{D}(::Type{LinearSpline{D}},BC::Type{None}, I::Symbol, x::Symbol)
    :()
end

function apply_bc{D}(::Type{LinearSpline{D}},BC::Type{Error}, I::Symbol, x::Symbol)
    quote
        if $x < $I.x[1] || $x > $I.x[end]
            throw(ArgumentError("$x is out of range"))
        end
    end
end

@inline function get_indices{D,S}(::Type{LinearSpline{D}}, G::Type{Uniform{1,S}}, I, x)
    dx = I.x[2]-I.x[1]
    i1 = clamp(Int(floor((x - I.x[1])/dx)) + 1, 1, S.parameters[1][1]-1)
    return i1, i1+1
end

@inline function get_indices{D,N}(::Type{LinearSpline{D}}, G::Type{Irregular{1,N}}, I, x)
    xr = searchsortedlast(I.x,x)
    i1 = clamp(xr, 1, N-1)
    return i1, i1+1
end

@inline function get_value{D}(::Type{LinearSpline{D}}, I, x, inds)
    i1, i2 = inds
    return I.y[i1] + I.m[i1]*(x - I.x[i1])
end
