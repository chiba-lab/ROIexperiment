# (×)(x, y) = Base.Iterators.product(x, y)

# [1:1:10] × [1:1:10]

# #logical Operations
# (∨)(x, y) -> x || y
# (∧)(x, y) -> x && y
# #others implemented in base


# struct mystruct{T,N}
#     a::T
#     b::N
# end


# #algebraic operations
# import Base.Iterators.product

# (×)(x, y) = Base.Iterators.product((x...), (y...))
# (⊗)(f::Function, 𝐗) = map(x -> (f(x), x), 𝐗)
# import Base.(/)(X, R) =
#     T = [1:1:10] × [2:1:5] |> collect
# t = ((X, Y) -> [(x, y) for x in X, y in Y])([(:a => 1), (:a => 2)], [(:b => 3), (:b => 4)])

# nt = NamedTuple.(t)
# using TypedTables

# Table(nt)



# t = Table([(friend1 = (name = "Alice", age = 25), friend2 = (name = "Bob", age = 42))], [(friend1 = (name = "Charlie", age = 37), friend2 = (name = "Bob", age = 42))])




# # @generated function myprod(A::Array{T,N}, B::{U, V}) where {T,N,U,V} 
# #     quote
# #         s = zero()
# #         @nloops $N i A begin
# #             s += @nref $N A i
# #         end
# #         s
# #     end
# # enMeta.show_sexpr(:(mystruct(1 + 2, 2 + 3)))





# #Sets Operations
# ⊂(x, y) = (x) ∩ (y) == (x)
# ⊃(x, y) = (x) ∪ (y) == (y)




# # (×)(x, y) = Base.Iterators.product(x, y)

# # [1:1:10] × [1:1:10]


# using AbstractTrees
# import AbstractTrees: children

# function AbstractTrees.children(tree::T) where {T<:Observation}
#     [ismissing(getfield(tree, Symbol(f))) ? "Missing" : getfield(tree, Symbol(f)) for f in fieldnames(T)]
# end

# function AbstractTrees.printnode(io::IO, node::BehavioralEvent)
#     print(io, "Event")
# end

# function AbstractTrees.printnode(io::IO, node::Session)
#     print(io, "Session")
# end

# function AbstractTrees.printnode(io::IO, node::Behavior)
#     print(io, node.name)
# end
# function AbstractTrees.printnode(io::IO, node::Rat)
#     print(io, node.name)
# end
# function AbstractTrees.printnode(io::IO, node::Date)
#     print(io, node)
# end

# import DataStructures: SortedMultiDict

# struct Foo2
#     properties::Dict{Any,Any}
# end

# (f::Foo2)(x) = getproperty(f, x)



# Foo2() = Foo2(Dict{Any,Any}())

# Base.getproperty(x::Foo2, property::Symbol) = getfield(x, :properties)[property]
# Base.setproperty!(x::Foo2, property::Symbol, value) = getfield(x, :properties)[property] = value
# Base.propertynames(x::Foo2) = keys(getfield(x, :properties))

# using DomainSet
# import DomainSets: Domain
# using InverseFunctions


# abstract type Map end

# function dom(m::Map)::Domain end
# function codom(m::Map)::Domain end
# function fun(m::Map)::Function end
# function inverse(m::Map)::Map end


# function image(m::Map, x)
#     @assert y ⊆ dom(map) "Undefined: input outside domain"
#     fun(m).(x)
# end

# function preimage(m::Map, y)
#     @assert y ⊆ codom(map) "Undefined: output outside codomain"
#     fun(m).(y)
# end





# (1, 2) ⊆ [1, 2, 3]
# map(x -> x * 2, s)

# using Dates
# struct Dates3 <: Domain{Date}
#     Domain
#     Dates3() = Domain{Date}()
# end

# const DateDom = Domain(Date)

