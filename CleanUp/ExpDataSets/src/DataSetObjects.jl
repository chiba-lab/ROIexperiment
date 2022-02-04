abstract type AbstractDataSetObject end 

function id(o::AbstractDataSetObject) end

function data(o::AbstractDataSetObject) end 

function source(o::AbstractDataSetObject) end

function meta(x) 
    return  nothing 
end

function meta(o::AbstractDataSetObject) 
    return meta(o), meta.(source(o))
end



















