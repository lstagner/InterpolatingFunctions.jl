immutable BilinearSpline <: AbstractInterpolation
    BilinearSpline() = BilinearSpline
end

immutable BilinearSplineInterpolation{T}
    x::Vector{T}
    y::Vector{T}
    z::Array{T,2}
end

function default_settings{D,T}(::Type{BilinearSpline},axes::NTuple{D,T})
    nx = length(axes[1])
    ny = length(axes[2])
    return Irregular{2,Val{(nx,ny)}}, Nearest
end

function BilinearSpline(x::AbstractVector, y::AbstractVector, z::AbstractArray)
  issorted(x) || throw(ArgumentError("x points must be in ascending order"))
  issorted(y) || throw(ArgumentError("y points must be in ascending order"))
  x = convert(Array{Float64,1},x)
  y = convert(Array{Float64,1},y)
  z = convert(Array{Float64,2},z)

  return BilinearSplineInterpolation{eltype(x)}(x,y,z)
end

@inline function preprocess(::Type{BilinearSpline}, axes, f)
    return BilinearSpline(axes...,f)
end

function apply_bc(::Type{BilinearSpline},BC::Type{Nearest}, I::Symbol, x::Symbol)
    xbc = apply_nearest(I,1, x, 1)
    ybc = apply_nearest(I,2, x, 2)
    quote
        $(xbc)
        $(ybc)
    end
end

function apply_bc(::Type{BilinearSpline},BC::Type{Reflect}, I::Symbol, x::Symbol)
    xbc = apply_reflect(I,1,x,1)
    ybc = apply_reflect(I,2,x,2)
    quote
        $(xbc)
        $(ybc)
    end
end

function apply_bc{V}(::Type{BilinearSpline},BC::Type{Value{V}}, I::Symbol, x::Symbol)
    xbc = apply_value(BC,I,1,x,1)
    ybc = apply_value(BC,I,2,x,2)
    quote
        $(xbc)
        $(ybc)
    end
end

function apply_bc(::Type{BilinearSpline},BC::Type{None}, I::Symbol, x::Symbol)
    apply_none()
end

function apply_bc(::Type{BilinearSpline},BC::Type{Error}, I::Symbol, x::Symbol)
    xbc = apply_error(I,1,x,1)
    ybc = apply_error(I,2,x,2)
    quote
        $(xbc)
        $(ybc)
    end
end

@inline function get_indices{S}(::Type{BilinearSpline}, G::Type{Uniform{2,S}}, I, x)
    dx = I.x[2]-I.x[1]
    ix = clamp(Int(floor((x[1] - I.x[1])/dx)) + 1, 1, S.parameters[1][1]-1)
    dy = I.y[2]-I.y[1]
    iy = clamp(Int(floor((x[2] - I.y[1])/dy)) + 1, 1, S.parameters[1][2]-1)
    return ix, iy
end

@inline function get_indices{S}(::Type{BilinearSpline}, G::Type{Irregular{2,S}}, I, x)
    xr = searchsortedlast(I.x,x[1],Base.Order.Forward)
    ix = clamp(xr, 1, S.parameters[1][1]-1)
    yr = searchsortedlast(I.y,x[2],Base.Order.Forward)
    iy = clamp(yr, 1, S.parameters[1][2]-1)
    return ix, iy
end

@inline function get_value(::Type{BilinearSpline}, I, x, inds)
    ix, iy = inds
    xp = x[1]
    yp = x[2]
    x1, x2 = I.x[ix], I.x[ix+1]
    y1, y2 = I.y[iy], I.y[iy+1]
    dxdy = (x2-x1)*(y2-y1)

    b11 = (x2 - xp)*(y2 - yp)
    b21 = (xp - x1)*(y2 - yp)
    b12 = (x2 - xp)*(yp - y1)
    b22 = (xp - x1)*(yp - y1)
    v = b11*I.z[ix,iy] + b12*I.z[ix,iy+1] + b21*I.z[ix+1,iy] + b22*I.z[ix+1,iy+1]
    return v/dxdy
end
