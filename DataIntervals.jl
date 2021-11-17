abstract type AbstractDataInterval <: AbstractDataSetObject end

function start_time(i::AbstractDataInterval) end

function end_time(i::AbstractDataInterval) end

function duration(i::AbstractDataInterval)
  return end_time(i) - start_time(i)
end



function Base.in(i₁::AbstractDataInterval, i₂::AbstractDataInterval) 
    if i₁.file != i₂.file
        return false
    elseif (i₁.start_time ≥ i₂.start_time) && (i₁.end_time ≤ i₂.end_time)
        return true
    else
        return false
    end
end

# function base.in(i₁::AbstractDataInterval, i₂::AbstractDataInterval) 
#     if i₁.file != i₂.file
#         return false
#     elseif (i₁.start_time ≤ i₂.start_time) || (i₁.end_time ≥ i₂.end_time)
#         return true
#     else
#         return false
#     end
# end

function Base.:∈(i::AbstractDataInterval, d::AbstractDataSetObject) 
    if d.file != i.file
        return false
    elseif (i.start_time ≥ start_time(d)) && (end_time(d) ≤ i.end_time)
        return true
    else
        return false
    end
end


function Base.:∋(i::AbstractDataSetObject, d::AbstractDataInterval) 
  if d.file != i.file
      return false
  elseif (start_time(d) ≥ start_time(i)) && (end_time(i) ≤ end_time(d))
      return true
  else
      return false
  end
end

# function base.in(d::AbstractDataSetObject, i::AbstractDataInterval) 
#     if d.file != i.file
#         return false
#     elseif (i.start_time ≤ start_time(d)) && (end_time(d) ≥ i.end_time)
#         return true
#     else
#         return false
#     end
# end

# function getinterval(i::AbstractDataInterval, d::AbstractExpDataType, pre=0, post=0) end





