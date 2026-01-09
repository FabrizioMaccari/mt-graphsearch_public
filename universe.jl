# needed external packages 
using JSMDInterfaces.Errors
using CairoMakie
using ColorSchemes
using GLMakie
using CSV
using JSON3 
using IterTools
using LinearAlgebra
using MixedTypesContainers
using SciMLBase: init, reinit!, solve!
using StaticArrays
using NonlinearSolve
using Graphs, SimpleWeightedGraphs  #, StaticGraphs
using Tables
using Printf
using DataFrames

include(joinpath("src", "ModelObjects.jl"))
include(joinpath("src", "Conversions.jl"))
include(joinpath("src", "Transfers.jl"))
include(joinpath("src", "MoonsSetup.jl"))
include(joinpath("src", "DesignSpaceDiscr.jl"))
include(joinpath("src", "HoppingGraph.jl"))
include(joinpath("src", "Wrappers.jl"))
include(joinpath("src", "Plots.jl"))


println("Universe defined")
