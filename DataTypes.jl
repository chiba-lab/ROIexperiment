abstract type AbstractDataType end


abstract type AbstractArrayData <: AbstractDataType end

function Base.values(dt::AbstractArrayData)
end

function axes(dt::AbstractArrayData)
end

function axes(dt::AbstractArrayData, ax)
    return axes(dt::AbstractArrayData)[ax]
end

function shape(dt::AbstractArrayData) 
    return shape(values(dt))
end

function size(dt::AbstractArrayData) 
    return size(values(dt))
end

# function (da::T)(axs...)  where T<:AbstractArrayData
#     return T(da, axs)
# end

# function (da::T)(newdata)  where T<:AbstractArrayData
#     return T(da, newdata)
# end

# function (da::T)(newdata, axs...)  where T<:AbstractArrayData
#     return T(da, newdata, axs)
# end



Dict([1,2,3].=> [4,5,6])













