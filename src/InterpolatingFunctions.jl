module InterpolatingFunctions

include("interpolation.jl")
include("cubic_spline.jl")
include("polyharmonic_spline.jl")
include("linear.jl")
include("polynomial.jl")

export interpolate
export Uniform, Irregular
export Nearest, Value, Reflect, Error, None
export CubicSpline, Monotonic, NonMonotonic
export PolyharmonicSpline, polyharmonic
export Polynomial, LinearSpline

end # module
