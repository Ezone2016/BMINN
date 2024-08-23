# Initialize battery parameters
function init_battery_params(M, V; Npa=30, Na=5, Ns=5, Nc=5, Np=5)
    theta = [60000.0, 27.5, 1e4, 50.0]  # theta1, theta2, theta3, theta4
    K = 0.05
    β = 2.0
    tol = 1e-4
    MV = sparse(M * V) # Mass matrix
    return BatteryParams(Npa, Na, Ns, Nc, Np, M, V, theta, K, β, tol, MV)
end

# Function to compute SOC
@inline function Cavg(u::AbstractArray, params::BatteryParams)
    @views sum(params.MV * u, dims=1) ./ (4/3 * pi)
end

# Define functions to setup neural networks
function setup_neural_network(;n=4)
    myNN = Lux.Chain(Lux.Dense(1, n), Lux.Dense(n, 2n, tanh), Lux.Dense(2n, n, tanh), Lux.Dense(n, 1))
    rng = Random.default_rng()
    p1, st = Lux.setup(rng, myNN)
    return myNN, st
end

# Function for solving ODE problems (De-intercalation & Intercalation)
function solve_ode(prob, γ, solver, callback, tol)
    solve(remake(prob; p = γ), solver, callback=callback, reltol=tol, abstol=tol)
end

# Function for computing the voltage curves
function compute_voltage_curve(x_de, x_in, Iapp, θθ)
    fde = kB * T / eV * f_muN_de(x_de, θθ)
    I0_de = @. 3 * θθ.theta2 * sqrt(x_de[end, :] * (1 - x_de[end, :]))
    V_de = -fde .- 2 * kB * T / eV * asinh.(Iapp ./ I0_de * θθ.theta3) .- Iapp * θθ.theta4

    fin = kB * T / eV * f_muN_in(x_in,θθ)
    I0_in = @. 3 * θθ.theta2 * sqrt(x_in[end, :] * (1 - x_in[end,:]))
    V_in = -fin .+ 2 * kB * T / eV * asinh.(Iapp ./ I0_in * θθ.theta3) .+ Iapp * θθ.theta4
    return V_de, V_in
end

# Plotting function with modular structure
function plot_voltage_40(x_de, V_de, x_in, V_in, params, exp_data_de, exp_data_in)
    plt = plot(x_de.t, V_de, label="C/40 exp. de-intercalation", linecolor=:blue, linestyle=:dash, dpi=150, ylims=(0.0, 0.6))
    plot!(plt, x_in.t, V_in, label="C/40 exp. intercalation", linecolor=:red, linestyle=:dash)
    plot!(plt, x_de.t, exp_data_de, label="C/40 simulated de-intercalation", linecolor=:blue)
    plot!(plt, x_in.t, exp_data_in, label="C/40 simulated intercalation", linecolor=:red)
    return plt
end

function plot_voltage_20(plt, x_de, V_de, x_in, V_in, params, exp_data_de, exp_data_in)
    plot!(plt, x_de.t, V_de, label="C/20 exp. de-intercalation", linecolor=:blue, linestyle=:dot, dpi=150, ylims=(0.0, 0.6))
    plot!(plt, x_in.t, V_in, label="C/20 exp. intercalation", linecolor=:red, linestyle=:dot)
    plot!(plt, x_de.t, exp_data_de, label="C/20 simulated de-intercalation", linecolor=:blue, lw=4, alpha = 0.5)
    plot!(plt, x_in.t, exp_data_in, label="C/20 simulated intercalation", linecolor=:red, lw=4, alpha = 0.5)
    return plt
end

# Discrete callback setup to avoid excessive termination checks
function setup_callbacks(q0, q1, params)
    condition_de(u, t, integrator) = sum(params.MV * u) / (4/3 * pi) <= q0
    condition_in(u, t, integrator) = sum(params.MV * u) / (4/3 * pi) >= q1
    affect!(integrator) = terminate!(integrator)
    Vcb_de = DiscreteCallback(condition_de, affect!)
    Vcb_in = DiscreteCallback(condition_in, affect!)
    return Vcb_de, Vcb_in
end

# Define main plotting and solving function
function Initialize_BMINN(params, solver, γ)
    # Setup battery params, callbacks, and initial conditions
    tspan = (0.0, 1e6)
    y0_de = ones(params.Npa + 1) * 0.87
    y0_in = ones(params.Npa + 1) * 0.08
    prob_CC_de = ODEProblem(ODEFunction(get_phis_de!), y0_de, tspan, γ)
    prob_CC_in = ODEProblem(ODEFunction(get_phis_in!), y0_in, tspan, γ)

    # Solve steady state problems
    global Iapp = 0.0;
    probSS_de = SteadyStateProblem(ODEFunction(get_phis_de!), y0_de, γ)
    probSS_in = SteadyStateProblem(ODEFunction(get_phis_in!), y0_in, γ)
    solss_de = solve(probSS_de, DynamicSS(solver), dt=0.1)
    solss_in = solve(probSS_in, DynamicSS(solver), dt=0.1)

    # Update initial condition with steady state
    prob_CC_de = remake(prob_CC_de; u0 = solss_de.u)
    prob_CC_in = remake(prob_CC_in; u0 = solss_in.u)
    return prob_CC_de, prob_CC_in
end

function solve_BMINN(γ, solver, params)
   # Set applied current and solve ODE problems
   prob_CC_de, prob_CC_in = Initialize_BMINN(params, solver, γ)
   Vcb_de, Vcb_in = setup_callbacks(0.08, 0.87, params)

   # Compute voltage curves and plot the results
   global Iapp = -5.0e-3 /39 
   x_de = solve_ode(prob_CC_de, γ, solver, Vcb_de, params.tol)
   x_in = solve_ode(prob_CC_in, γ, solver, Vcb_in, params.tol)
   V_de, V_in = compute_voltage_curve(x_de, x_in, Iapp, γ)
   plt = plot_voltage_40(x_de, V_de, x_in, V_in, params, muh60_de.(x_de.t * 3600), muh60_in.(x_in.t * 3600))

   global Iapp = -5.0e-3 /39 * 2
   x_de = solve_ode(prob_CC_de, γ, solver, Vcb_de, params.tol)
   x_in = solve_ode(prob_CC_in, γ, solver, Vcb_in, params.tol)
   V_de, V_in = compute_voltage_curve(x_de, x_in, Iapp, γ)
   plot_voltage_20(plt, x_de, V_de, x_in, V_in, params, muh30_de.(x_de.t * 3600), muh30_in.(x_in.t * 3600))
   display(plt)
end