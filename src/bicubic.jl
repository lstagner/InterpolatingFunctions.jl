immutable BicubicSpline <: AbstractInterpolation
    BicubicSpline() = BicubicSpline
end

immutable BicubicSplineInterpolation{T}
    x::Vector{T}
    y::Vector{T}
    z::Array{T,2}
    A::Array{T,4}
end

function default_settings{D,T}(::Type{BicubicSpline},axes::NTuple{D,T})
    nx = length(axes[1])
    ny = length(axes[2])
    return Irregular{2,Val{(nx,ny)}}, Nearest
end

function deriv(x,y,i)
    nx = length(x)
    if i == 1
        dy = y[i+1] - y[i]
        return dy
    end
    if i == nx
        dy = y[i] - y[i-1]
        return dy
    end
    dy1 = y[i+1] - y[i]
    dy2 = y[i] - y[i-1]
    return 0.5*(dy1 + dy2)
end

function BicubicSpline(x::AbstractVector,y::AbstractVector, z::AbstractArray)
    issorted(x) || throw(ArgumentError("x points must be in ascending order"))
    issorted(y) || throw(ArgumentError("y points must be in ascending order"))
    x = convert(Array{Float64,1},x)
    y = convert(Array{Float64,1},y)
    z = convert(Array{Float64,2},z)
    nx = length(x)
    ny = length(y)
    nx == size(z,1) || throw(ArgumentError("length(x) != size(z,1)"))
    ny == size(z,2) || throw(ArgumentError("length(y) != size(z,2)"))

    px = similar(z)
    py = similar(z)
    pxy = similar(z)

    for i=1:nx
        for j = 1:ny
            px[i,j] = deriv(x,z[:,j],i)
            py[i,j] = deriv(y,z[i,:],j)
        end
        for j=1:ny
            pxy[i,j] = deriv(y,px[i,:],j)
        end
    end
    C = [1 0 0 0; 0 0 1 0; -3 3 -2 -1; 2 -2 1 1]
    A = zeros(4,4,nx,ny)
    for i=1:nx-1
        for j=1:ny-1
            F = [z[i,j]    z[i,j+1]    py[i,j]    py[i,j+1];
                 z[i+1,j]  z[i+1,j+1]  py[i+1,j]  py[i+1,j+1];
                 px[i,j]   px[i,j+1]   pxy[i,j]   pxy[i,j+1];
                 px[i+1,j] px[i+1,j+1] pxy[i+1,j] pxy[i+1,j+1]]
            A[:,:,i,j] = C*(F*C')
        end
    end

    x,y,z, A = promote(x,y,z,A)
    return BicubicSplineInterpolation{eltype(x)}(x,y,z,A)
end

@inline function preprocess(::Type{BicubicSpline}, axes, f)
    return BicubicSpline(axes...,f)
end

function apply_bc(::Type{BicubicSpline},BC::Type{Nearest}, I::Symbol, x::Symbol)
    xbc = apply_nearest(I, 1, x, 1)
    ybc = apply_nearest(I, 2, x, 2)
    quote
        $(xbc)
        $(ybc)
    end
end

function apply_bc(::Type{BicubicSpline},BC::Type{Reflect}, I::Symbol, x::Symbol)
    xbc = apply_reflect(I,1,x,1)
    ybc = apply_reflect(I,2,x,2)
    quote
        $(xbc)
        $(ybc)
    end
end

function apply_bc{V}(::Type{BicubicSpline},BC::Type{Value{V}}, I::Symbol, x::Symbol)
    xbc = apply_value(BC,I,1,x,1)
    ybc = apply_value(BC,I,2,x,2)
    quote
        $(xbc)
        $(ybc)
    end
end

function apply_bc(::Type{BicubicSpline},BC::Type{None}, I::Symbol, x::Symbol)
    apply_none()
end

function apply_bc(::Type{BicubicSpline},BC::Type{Error}, I::Symbol, x::Symbol)
    xbc = apply_error(I,1,x,1)
    ybc = apply_error(I,2,x,2)
    quote
        $(xbc)
        $(ybc)
    end
end

@inline function get_indices{S}(::Type{BicubicSpline}, G::Type{Uniform{2,S}}, I, x)
    dx = I.x[2]-I.x[1]
    ix = clamp(Int(floor((x - I.x[1])/dx)) + 1, 1, S.parameters[1][1]-1)
    dy = I.y[2]-I.y[1]
    iy = clamp(Int(floor((y - I.y[1])/dy)) + 1, 1, S.parameters[1][2]-1)
    return ix, iy
end

@inline function get_indices{S}(::Type{BicubicSpline}, G::Type{Irregular{2,S}}, I, x)
    xr = searchsortedlast(I.x,x[1],Base.Order.Forward)
    ix = clamp(xr, 1, S.parameters[1][1]-1)
    yr = searchsortedlast(I.y,x[2],Base.Order.Forward)
    iy = clamp(yr, 1, S.parameters[1][2]-1)
    return ix, iy
end

@inline function get_value(::Type{BicubicSpline}, I, x, inds)
    ix, iy = inds
    A = I.A[:,:,ix,iy]
    dx = I.x[ix+1] - I.x[ix]
    xx = (x[1] - I.x[ix])/dx
    dy = I.y[iy+1] - I.y[iy]
    yy = (x[2] - I.y[iy])/dy
    v = (1.0) *(A[1,1] + A[1,2]*yy + A[1,3]*(yy^2) + A[1,4]*(yy^3)) +
        (xx)  *(A[2,1] + A[2,2]*yy + A[2,3]*(yy^2) + A[2,4]*(yy^3)) +
        (xx^2)*(A[3,1] + A[3,2]*yy + A[3,3]*(yy^2) + A[3,4]*(yy^3)) +
        (xx^3)*(A[4,1] + A[4,2]*yy + A[4,3]*(yy^2) + A[4,4]*(yy^3))
    return v
end
