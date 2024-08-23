using Plots, SciMLSensitivity
using ModelingToolkit, DifferentialEquations
using LinearAlgebra, ToeplitzMatrices
using Plots.PlotMeasures
using Parameters: @unpack
using ForwardDiff
using SparseArrays, Roots, QuadGK, Trapz, OptimizationOptimJL, OptimizationOptimisers
using Random, Lux, ReverseDiff, ComponentArrays, Optimization, DataInterpolations
using OptimizationFlux, ReverseDiff, Zygote, DiffEqFlux, LaTeXStrings, SpecialFunctions
using BSON: @save, @load

include("uti_cheb.jl")
include("VIdata.jl")

# Physics-based model specifications
begin
    Npa = 30 #pseudo dimension
    Na = Ns = Nc = Np = 5   #dummy index in thickness dimension
    include("para_battery.jl")
    include("uti_pf.jl")  
    solver = Rodas4()
    sensealg = InterpolatingAdjoint(autojacvec=ZygoteVJP(true)) 
    theta1 = 60000
    theta2 = 27.5
    theta3 = 1e4
    theta4 = 50
    K = 0.05
    β = 2
    tol = 1e-4
    MV = sparse(MM*VV)
    Cavg(u) = (sum(MV * Array(u), dims = 1)./(4/3*pi))[:]
end

# Neural network specifications
begin
    n = 4
    myNN = Lux.Chain(Lux.Dense(1,n),Lux.Dense(n,2*n,tanh),Lux.Dense(2*n,n,tanh),Lux.Dense(n,1))
    p1, st = Lux.setup(Random.default_rng(), myNN)
    p = ComponentArray{Float32}()
    p = ComponentArray(p;p1)
    c = Float32.(collect(range(0.0, stop=1.0, length=200)))
    cc = reshape(c,1,:)
    sc(p) = myNN(cc, p, st)[1][:] 
    @load "./Data/dl_training60.bson" learned 
    p1 = learned
    γ = ComponentArray{Float32}()
    γ = ComponentArray(γ;p1)
    γ = ComponentArray(γ;theta1)
    γ = ComponentArray(γ;theta2)
    γ = ComponentArray(γ;theta3)
    γ = ComponentArray(γ;theta4)
    γ = ComponentArray(γ;K)
    γ = ComponentArray(γ;β)
end

# Plot training results
function plotalpha(θθ)
    x_de = solve(remake(prob_CC_de; p = θθ), solver,
        callback=Vcb_de, reltol = tol, abstol = tol);
    fde = const1*f_muN_de(x_de,θθ)
    I0_de = @. 3*θθ.theta2*sqrt(x_de[end,:]*(1-x_de[end,:]))
    V_de = -fde .- 2*const1*asinh.(Iapp./I0_de*θθ.theta3) .- Iapp*θθ.theta4

    x_in = solve(remake(prob_CC_in; p = θθ), solver, 
        callback=Vcb_in, reltol = tol, abstol = tol);
    fin = const1*f_muN_in(x_in,θθ)
    I0_in = @. 3*θθ.theta2*sqrt(x_in[end,:]*(1-x_in[end,:]))
    V_in = -fin .+ 2*const1*asinh.(Iapp./I0_in*θθ.theta3) .+ Iapp*θθ.theta4
    plt = plot(x_de.t, V_de, label = "C/40 simulated de-intercalation", 
        linecolor=:blue, linestyle=:dash, dpi=150, ylims=(0.0, 0.6))
    Plots.plot!(plt, x_in.t, V_in, label = "C/40 simulated intercalation", 
        linecolor=:red, linestyle=:dash)
    Plots.plot!(plt, x_de.t, muh60_de.(x_de.t*3600), 
        label = "C/40 exp. de-intercalation", linecolor=:blue)
    Plots.plot!(plt, x_in.t, muh60_in.(x_in.t*3600), 
        label = "C/40 exp. intercalation", linecolor=:red)
    return plt
end

# Plot testing results
function plot30(θθ,plt)
    x_de = solve(remake(prob_CC_de; p = θθ), solver,
        callback=Vcb_de, reltol = tol, abstol = tol);
    fde = const1*f_muN_de(x_de,θθ)
    I0_de = @. 3*θθ.theta2*sqrt(x_de[end,:]*(1-x_de[end,:]))
    V_de = -fde .- 2*const1*asinh.(Iapp./I0_de*θθ.theta3) .- Iapp*θθ.theta4

    x_in = solve(remake(prob_CC_in; p = θθ), solver, 
        callback=Vcb_in, reltol = tol, abstol = tol);
    fin = const1*f_muN_in(x_in,θθ)
    I0_in = @. 3*θθ.theta2*sqrt(x_in[end,:]*(1-x_in[end,:]))
    V_in = -fin .+ 2*const1*asinh.(Iapp./I0_in*θθ.theta3) .+ Iapp*θθ.theta4
    plot!(plt, x_de.t, V_de, label = "C/20 simulated de-intercalation", 
        linecolor=:blue, alpha = 0.5, linestyle=:dash, dpi=150, ylims=(0.0, 0.6))
    Plots.plot!(plt, x_in.t, V_in, label = "C/20 simulated intercalation", 
        linecolor=:red, alpha = 0.5, linestyle=:dash)
    Plots.plot!(plt, x_de.t, muh30_de.(x_de.t*3600), alpha = 0.5, 
        label = "C/20 exp. de-intercalation", linecolor=:blue)
    Plots.plot!(plt, x_in.t, muh30_in.(x_in.t*3600), alpha = 0.5, 
        label = "C/20 exp. intercalation", linecolor=:red)
    return plt
end

# BMINN setup and initialization
begin 
    q0 = 0.08
    q1 = 0.87
    condition(u,t,integrator) = sum(MV * u) /(4/3*pi) <= q0
    affect!(integrator) = terminate!(integrator)
    Vcb_de = DiscreteCallback(condition,affect!)
    const1 = kB*T/eV
    condition1(u,t,integrator) = sum(MV * u) /(4/3*pi) >= q1
    affect1!(integrator) = terminate!(integrator)
    Vcb_in = DiscreteCallback(condition1,affect1!)
    tspan = (0, 1e6)
    tol = 5e-4
    y0_de = ones(Npa+1) * q1
    y0_in = ones(Npa+1) * q0
    prob_CC_de = ODEProblem(ODEFunction(get_phis_de!), y0_de, tspan, α)
    prob_CC_in = ODEProblem(ODEFunction(get_phis_in!), y0_in, tspan, α)
end

# Solve
begin
    solver = Rosenbrock23(autodiff=false)
    θθ = γ
    α = γ
    Iapp = 0.0;
    probSS = SteadyStateProblem(ODEFunction(get_phis_de!), y0_de, α)
    solss_de = solve(probSS, DynamicSS(solver), dt = 0.1)
    prob_CC_de = ODEProblem(ODEFunction(get_phis_de!), solss_de.u, tspan, α)

    Iapp = -5.0e-3 /39 
    @time x_de = solve(remake(prob_CC_de; p = α), solver,saveat=0.1,sensealg=sensealg, 
        callback=Vcb_de, reltol = tol, abstol = tol);
    fde = const1*f_muN_de(x_de,α)
    I0_de = @. 3*θθ.theta2*sqrt(x_de[end,:]*(1-x_de[end,:]))
    V_de = -fde .- 2*const1*asinh.(Iapp./I0_de*θθ.theta3) .- Iapp*θθ.theta4
    plot(x_de.t, V_de, ylims=(0.0,0.6))
    plot!(x_de.t, muh60_de.(x_de.t*3600))

    Iapp = 0.0;
    probSS = SteadyStateProblem(ODEFunction(get_phis_in!), y0_in, α)
    solss_in = solve(probSS, DynamicSS(solver), dt = 0.1)
    prob_CC_in = ODEProblem(ODEFunction(get_phis_in!), solss_in.u, tspan, α)

    Iapp = -5.0e-3 /39 
    @time x_in = solve(remake(prob_CC_in; p = α), solver,saveat=0.1,sensealg=sensealg, 
        callback=Vcb_in, reltol = tol, abstol = tol);
    fin = const1*f_muN_in(x_in,α)
    I0_in = @. 3*θθ.theta2*sqrt(x_in[end,:]*(1-x_in[end,:]))
    V_in = -fin .+ 2*const1*asinh.(Iapp./I0_in*θθ.theta3) .+ Iapp*θθ.theta4
    plot!(x_in.t, V_in)
    plot!(x_in.t, muh60_in.(x_in.t*3600))

    tspan_de = (0.0, floor(x_de.t[end]; digits=1))
    tspan_in = (0.0, floor(x_in.t[end]; digits=1))

    x_de = solve(remake(prob_CC_de; p = θθ, tspan = tspan_de), solver, saveat=0.1,
        sensealg=sensealg, reltol = tol, abstol = tol);
    x_in = solve(remake(prob_CC_in; p = θθ, tspan = tspan_in), solver, sensealg=sensealg,
        saveat=0.1, reltol = tol, abstol = tol);
    tde = x_de.t
    tin = x_in.t
    Mde = muh60_de.(tde*3600)
    Min = muh60_in.(tin*3600)
end

function loss(θθ)
    x_de = solve(remake(prob_CC_de; p = θθ, tspan = tspan_de), solver, saveat=0.1,
        sensealg=sensealg, reltol = tol, abstol = tol);
    fde = const1*f_muN_de(x_de,θθ)
    I0_de = @. 3*θθ.theta2*sqrt(x_de[end,:]*(1-x_de[end,:]))
    V_de = -fde .- 2*const1*asinh.(Iapp./I0_de*θθ.theta3) .- Iapp*θθ.theta4

    x_in = solve(remake(prob_CC_in; p = θθ, tspan = tspan_in), solver, sensealg=sensealg,
        saveat=0.1, reltol = tol, abstol = tol);
    fin = const1*f_muN_in(x_in,θθ)
    I0_in = @. 3*θθ.theta2*sqrt(x_in[end,:]*(1-x_in[end,:]))
    V_in = -fin .+ 2*const1*asinh.(Iapp./I0_in*θθ.theta3) .+ Iapp*θθ.theta4
    return sum(abs2, Mde .- V_de) .+ sum(abs2, Min .- V_in) 
end

# Compute the gradients and plot the errors
begin
    @time grad_reverse = ReverseDiff.gradient(loss, α); #ReverseDiff or Zygote
    @time grad_forward = ForwardDiff.gradient(loss, α);
    plot(collect(grad_reverse[1]) .- collect(grad_forward))
    plotalpha(γ)
    loss(γ)
end

# Plot the voltage curves
begin
    Iapp = -5.0e-3 /39
    @time loss(α)
    plt = plotalpha(α)
    
    Iapp = -5.0e-3 /39 *2
    plot30(α,plt)
    display(plt)    
end


function sgld!(∇L, θᵢ, t) 
    ϵ = a1*(b1 + t)^-γ1
    η = ϵ.*randn(size(θᵢ))
    Δθᵢ = .5ϵ*∇L .+ η
    θᵢ .-= Δθᵢ
end