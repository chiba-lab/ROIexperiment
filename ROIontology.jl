using TypedTables, Unitful, Dates, Parameters, Glob, CSV, StructArrays, TypedTables, JLD2, IntervalSets, AutoHashEquals
import FourierAnalysis as FA
import Unitful: Hz, s as 𝓈



##===========================================================================================================================

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
@auto_hash_equals struct Rat <: Agent
    name
end
@auto_hash_equals struct Robot <: Agent
    name
end
@auto_hash_equals struct Object <: Agent
    name
end

const EG7 = Rat("EG7")
const RRSD18 = Rat("RRSD18")
const RRSD28 = Rat("RRSD28")
const RRSD17 = Rat("RRSD17")


#Conditions 
struct TrialCondition{T}
    name
end

function agenttype(t::TrialCondition{T}) where T
    return T
end

const EE(T) = TrialCondition{T}("Empty Empty")
const OF(T) = TrialCondition{T}("Open Field")
const FR(T) = TrialCondition{T}("Free Roam")
const HBT(T) = TrialCondition{T}("Habituation")
const ITR(T) = TrialCondition{T}("Interaction")

#Behaviors 
@auto_hash_equals struct Behavior
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

const FBands = Dict(())

#
##==============================================================================
# Observation Epochs
abstract type Observation end

abstract type ObservationInterval <: Observation end
@auto_hash_equals struct Session <: ObservationInterval
    filename
    date
    start_time
    end_time
end

@auto_hash_equals struct Trial <: ObservationInterval
    session
    start_time
    end_time
    condition
end



##=============================================================================
#  Data Types
abstract type ObservationData <: Observation end

struct BehavioralEvent <: ObservationData
    session
    start_time
    end_time
    behavior
    actor
    receiver
end

struct LFPRecording <: ObservationData
    session
    start_time
    end_time
    rat
    region
    fs
    lfp
end

get_interval(x::LFPRecording, s, e) =  x.lfp[round(Int, (s - x.start_time) * 1010.1):round(Int, (e - x.start_time) * 1010.1)]
struct DataWindow 
    onset
    offset
    data::ObservationData
end

function get_data(w::DataWindow, pre, post) 
    get_interval(w.data, max(w.onset - pre, start_time(w.data)), min(w.offset + post, end_time(w.data)))
end



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



