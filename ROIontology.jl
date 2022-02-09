using TypedTables, Unitful, Dates, Parameters, FourierAnalysis, Glob, CSV, StructArrays, TypedTables, JLD2, IntervalSets
import Unitful: Hz, s as 𝓈


##==============================================================================
#Experiment Types
#
#Brain Regions
struct Region
    name
end

const AMG = Region("Amygdala")
const CA2 = Region("Ca2")
const MOB = Region("MOB")
const INS = Region("Insula")

#Agents
abstract type Agent end
struct Rat <: Agent
    name
end
struct Robot <: Agent
    name
end
struct Object <: Agent
    name
end


#Conditions 
struct TrialCondition{T}
    name
end

const EE(T) = TrialCondition{T}("Empty Empty")
const OF(T) = TrialCondition{T}("Open Field")
const FR(T) = TrialCondition{T}("Free Roam")
const HBT(T) = TrialCondition{T}("Habituation")
const ITR(T) = TrialCondition{T}("Interaction")

#Behaviors 
struct Behavior
    name
end

const Groom = Behavior("Groom")
const Rear = Behavior("Rear")
const Baseline = Behavior("Baseline")
const Immobility = Behavior("Immobility")
const Sniff = Behavior("Sniff")
const Top = Behavior("Top")
const Approach = Behavior("Approach")
const Retreat = Behavior("Retreat")



#
##==============================================================================
# Observation Epochs
abstract type Observation end

abstract type ObservationInterval <: Observation end
@with_kw struct Session <: ObservationInterval
    filename
    date
    start_time
    end_time
end

@with_kw struct Trial <: ObservationInterval
    session
    start_time
    end_time
    condition
end



##=============================================================================
#  Data Types
abstract type ObservationData <: Observation end

@with_kw struct BehavioralEvent <: ObservationData
    session
    start_time
    end_time
    behavior
    actor
    receiver
end

@with_kw struct LFPRecording <: ObservationData
    session
    start_time
    end_time
    rat
    region
    fs
    lfp
end
get_interval_lfp(x::LFPRecording, s, e) = x.lfp[round(Int, (s-x.start_time)*1010.1 : (e-x.start_time)*1010.1)]

##==============================================================================
#  Time Handlers 
root_interval(t::T) where {T<:Session} = t
root_interval(t::T) where {T<:Observation} = t.session
date(t::T) where {T<:Session} = t.date
date(t::T) where {T<:Observation} = date(root_interval(t))
start_time(t::T) where {T<:Observation} = t.start_time
end_time(t::T) where {T<:Observation} = t.end_time
time_interval(t::T) where {T<:Observation} = start_time(t) .. end_time(t)
duration(t::T) where {T<:Observation} = end_time(t) - start_time(t)
⊂(x, y) = (x) ∩ (y) == (x)
in_interval(i::T, j::N) where {N,T} = root_interval(i) === root_interval(j) ? time_interval(i) ⊂ time_interval(j) : false
##=============================================================================
# Data Handlers
# abstract type DataShape end
# struct ArrayData <: DataShape end
# struct ScalarData <: DataShape end
# datashape(::T) where T = ScalarData() #fall-back, by default everything is scalar
# datashape(::LFPRecording) = ArrayData()

# data(x::T) where T = data(datashape(T), x)

# data(::ScalarData, x) = value(x)
# data(::ArrayData, x) = (value(x), dims(x))

# value(t::LFPRecording) = t.lfp
# dims(t::LFPRecording) = (:time, start_time(t):1/t.fs:end_time(t))

##==============================================================================
# Links
##==============================================================================



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



