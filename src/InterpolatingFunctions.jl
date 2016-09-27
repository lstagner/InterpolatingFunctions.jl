module InterpolatingFunctions

abstract AbstractInterpolation

include("grid.jl")
include("boundary.jl")
include("linear.jl")
include("bilinear.jl")
include("cubic.jl")
include("bicubic.jl")
include("polyharmonic.jl")
include("polynomial.jl")

typealias AI AbstractInterpolation
typealias AG AbstractGrid
typealias ABC AbstractBoundaryCondition

function interpolate_gen(I,G,BC)
    bc_expr = apply_bc(I, BC, :interp, :x)
    ifunc = quote
        interp = preprocess($I, axes, f)
        function InterpolatingFunction(args...)
            if length(args) == 1
                if typeof(args[1]) <: AbstractArray
                    x = args[1]
                else
                    x = [args[1]]
                end
            else
                x = collect(args)
            end
            $(bc_expr)

            inds = get_indices($I, $G, interp, x)
            val = get_value($I, interp, x, inds)
            return val
        end
    end
    return ifunc
end

function interpolate{N,T,I<:AI}(f::AbstractArray{T}, axes::NTuple{N,AbstractArray}, ::Type{I})
    G,BC = default_settings(I,axes)
    interpolate(f,axes,I,G,BC)
end

@generated function interpolate{N,T,I<:AI,G<:AG,BC<:ABC}(f::AbstractArray{T}, axes::NTuple{N,AbstractVector},::Type{I},::Type{G},::Type{BC})
    interpolate_gen(I,G,BC)
end

export interpolate, InterpolatingFunction
export Uniform, Irregular
export Nearest, Value, Reflect, Error, None
export CubicSpline, Monotonic, NonMonotonic
export PolyharmonicSpline, polyharmonic
export Polynomial, LinearSpline
export BicubicSpline, BilinearSpline

end # module
