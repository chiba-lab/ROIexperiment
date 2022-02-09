⊂(x, y) = (x) ∩ (y) == (x)
⊃(x, y) = (x) ∪ (y) == (y)




# (×)(x, y) = Base.Iterators.product(x, y)

# [1:1:10] × [1:1:10]


using AbstractTrees
import AbstractTrees: children

function AbstractTrees.children(tree::T) where {T<:Observation}
    [ismissing(getfield(tree, Symbol(f))) ? "Missing" : getfield(tree, Symbol(f)) for f in fieldnames(T)]
end

function AbstractTrees.printnode(io::IO, node::BehavioralEvent)
    print(io, "Event")
end

function AbstractTrees.printnode(io::IO, node::Session)
    print(io, "Session")
end

function AbstractTrees.printnode(io::IO, node::Behavior)
    print(io, node.name)
end
function AbstractTrees.printnode(io::IO, node::Rat)
    print(io, node.name)
end
function AbstractTrees.printnode(io::IO, node::Date)
    print(io, node)
end



