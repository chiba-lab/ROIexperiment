
using TypedTables, Unitful, AxisKeys, Dates, Parameters, FourierAnalysis, Glob, CSV, StructArrays, TypedTables, JLD2
import Unitful: Hz, s as 𝓈

##==============================================================================
#Experiment Types


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

##==============================================================================
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

##==============================================================================
#  Time Handlers 
root_interval(t::T) where {T<:Session} = t
root_interval(t::T) where {T<:Observation} = t.session
date(t::T) where {T<:Session} = t.date
date(t::T) where {T<:Observation} = date(root_interval(t))
start_time(t::T) where {T<:Observation} = t.start_time
end_time(t::T) where {T<:Observation} = t.end_time
duration(t::T) where {T<:Observation} = end_time(t) - start_time(t)
in_interval(::T, j::T) where {T<:Observation} = root_interval(i) === root_interval(j) ? (start_time(j) <= start_time(i)) && (end_time(i) >= end_time(j)) : false

##=============================================================================
# Data Handlers
abstract type DataShape end
struct ArrayData <: DataShape end
struct ScalarData <: DataShape end
datashape(t::T) where {T} = ScalarData()

data(t::DataShape) = data(datashape(t), t)
data(::ScalarData, t) = value(t)
data(::ArrayData, t) = KeyedArray(value(t), dims(t)...)

datashape(t::LFPRecording) = ArrayData()
value(t::LFPRecording) = t.lfp
dims(t::LFPRecording) = (time = start_time(t):1/t.fs:end_time(t))

##==============================================================================
# Data Filters
##==============================================================================

#







