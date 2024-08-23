using Plots, SciMLSensitivity, ModelingToolkit, DifferentialEquations
using LinearAlgebra, ToeplitzMatrices, Plots.PlotMeasures
using ForwardDiff, SparseArrays, Roots, QuadGK, Trapz, Random, Lux, ComponentArrays
using ReverseDiff, Zygote, DiffEqFlux
using BSON: @save, @load
using Parameters: @unpack

# Include external files
include("./equations/uti_cheb.jl")
include("VIdata.jl")
Npa = 30 #pseudo dimension
Na = Ns = Nc = Np = 5   #dummy index in thickness dimension
include("para_battery.jl")
include("./equations/uti_pf.jl")
include("./equations/uti_BMINN.jl")

struct BatteryParams
    Npa::Int
    Na::Int
    Ns::Int
    Nc::Int
    Np::Int
    MM::SparseMatrixCSC
    VV::SparseMatrixCSC
    theta::Vector{Float64}
    K::Float64
    β::Float64
    tol::Float64
    MV::SparseMatrixCSC
end

# Initialize solver and adjoint sensitivity algorithm
solver = Rodas4()
sensealg = InterpolatingAdjoint(autojacvec=ZygoteVJP(true))

# Initialize battery parameters
params = init_battery_params(sparse(MM), sparse(VV))

# Load neural network parameters
begin
    @load "dl_training.bson" learned
    myNN, st = setup_neural_network(n=4)
    p1 = learned
    γ = ComponentArray{Float32}()
    γ = ComponentArray(γ;p1)
    theta1 = params.theta[1]
    theta2 = params.theta[2]
    theta3 = params.theta[3]
    theta4 = params.theta[4]
    γ = ComponentArray(γ;theta1)
    γ = ComponentArray(γ;theta2)
    γ = ComponentArray(γ;theta3)
    γ = ComponentArray(γ;theta4)
    γ = ComponentArray(γ;params.K)
    γ = ComponentArray(γ;params.β) 
end

# Run the simulation and plot voltage curves
solve_BMINN(γ, solver, params)

