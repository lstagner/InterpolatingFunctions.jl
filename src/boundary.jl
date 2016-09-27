abstract AbstractBoundaryCondition

type Value{V} <: AbstractBoundaryCondition end
Value(V) = Value{Val{V}}
function apply_value{V}(BC::Type{Value{V}},I::Symbol, fi::Int, x::Symbol, xi::Int)
    quote
        xmin, xmax = extrema(getfield($I,$fi))
        if $x[$xi] < xmin || $x[$xi] > xmax
            return $V.parameters[1]
        end
    end
end

type Nearest <: AbstractBoundaryCondition end
function apply_nearest(I::Symbol, fi::Int, x::Symbol, xi::Int)
    :($x[$xi] = clamp($x[$xi], extrema(getfield($I,$fi))...))
end

type Reflect <: AbstractBoundaryCondition end
function apply_reflect(I::Symbol, fi::Int, x::Symbol, xi::Int)
    quote
        xmin, xmax = extrema(getfield($I,$fi))
        if $x[$xi] < xmin
            delX = xmax - xmin
            n = floor((xmin - $x[$xi])/delX)
            $x[$xi] = xmax - ((xmin - $x[$xi]) - n*delX)
        elseif $x[$xi] > xmax
            delX = xmax - xmin
            n = floor(($x[$xi] - xmax)/delX)
            $x[$xi] = xmin + ($x[$xi] - xmax - n*delX)
        end
    end
end

type None <: AbstractBoundaryCondition end
function apply_none()
    :()
end

type Error <: AbstractBoundaryCondition end
function apply_error(I::Symbol,fi::Int,x::Symbol,xi::Int)
    quote
        xmin, xmax = extrema(getfield($I,$fi))
        if $x[$xi] < xmin || $x[$xi] > xmax
            return throw(ArgumentError("Out of range"))
        end
    end
end
