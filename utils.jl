using Dates

function findnearest(A::AbstractArray,t) 
    findmin(abs(A-t))
end

function nearest(A::AbstractArray,t) 
    A[findnearest(A,t)]
end


function minstr2seconds(s)
    return sum(parse.(Float64, split(s, ":")) .* [60,1])
end

function fstr2date(s)
    sarr = parse.(Int, split(s, "_")[2:4])
    return Date(sarr[3] + 2000, sarr[1], sarr[2])
end

function genid(f,s,e)
    f*"_"*string(round(s,digits=1))*"_"*string(round(e,digits=1))
end




