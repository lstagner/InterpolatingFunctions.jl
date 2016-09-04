immutable PolyharmonicSpline{K,S} <: AbstractInterpolation
end

PolyharmonicSpline() = PolyharmonicSpline{2,Val{0}}
PolyharmonicSpline(K) = PolyharmonicSpline{K,Val{0}}
PolyharmonicSpline(K,S) = PolyharmonicSpline{K,Val{S}}

typealias ThinPlateSpline{S} PolyharmonicSpline{2,S}
ThinPlateSpline() = ThinPlateSpline{Val{0}}
ThinPlateSpline(S) = ThinPlateSpline{Val{S}}

function default_settings{D,T,K,S}(::Type{PolyharmonicSpline{K,S}},axes::NTuple{D,T})
    N = length(axes[1])
    return Irregular{D,N}, None
end

type PolyharmonicSplineInterpolation{T,D,N}
  coeff::Vector{T}
  values::Vector{T}
  centers::Array{T,2}
end

@inline function polyharmonicK(r,K)
    if iseven(K)
        p = r < 1.0 ? (r^(K-1))*log(r^r) : (r^K)*log(r)
    else
        p = r^K
    end
    return p
end

function preprocess{K,S}(phs::Type{PolyharmonicSpline{K,S}}, axes, f)
    values = convert(Array, f)
    centers = hcat(axes...)'
    n, m = size(centers)
    m != length(values) && throw(DimensionMismatch())

    M = zeros(m,m)
    N = zeros(m,n+1)

    for i=1:m
        N[i,1] = 1
        N[i,2:end] = centers[:,i]
        for j=i:m
            M[i,j] = polyharmonicK(norm(centers[:,i] .- centers[:,j]),K)
            M[j,i] = M[i,j]
        end
    end
    s = S.parameters[1]
    M = M .+ s.*eye(m)
    L = [[M N]; [N' zeros(n+1,n+1)]]

    w = pinv(L)*[values; zeros(n+1)]

    w, centers, values = promote(w, centers, values)
    return PolyharmonicSplineInterpolation{eltype(w),n,length(w)}(w,values,centers)
end

function apply_bc{K,S}(::Type{PolyharmonicSpline{K,S}}, ::Type{None}, I::Symbol, x::Symbol)
    :()
end

@inline function get_indices{K,S,D,N}(::Type{PolyharmonicSpline{K,S}}, G::Type{Irregular{D,N}}, I, x)
    return ()
end

@inline function get_value{K,S,T,D,N}(::Type{PolyharmonicSpline{K,S}}, I::PolyharmonicSplineInterpolation{T,D,N}, x, inds)
    v = zero(Float64)
    l = N - (D+1)
    centers = I.centers
    coeff = I.coeff
    @inbounds for i=1:l
        v = v + coeff[i]*polyharmonicK(norm(x .- centers[:,i]),K)
    end
    v = v + coeff[l+1]
    @inbounds for i=2:D+1
        v = v + coeff[l+i]*x[i-1]
    end
    return v
end
