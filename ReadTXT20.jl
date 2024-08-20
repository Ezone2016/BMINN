using DelimitedFiles, Dates, Plots
using Plots.PlotMeasures
using DataInterpolations

# Read the data from the file
data = readdlm("paper4/exp/cell1_20.txt")[2:end-1,1:3]

function time_to_seconds(time_strings::Vector{String})
    seconds_array = Int[]
    for time_string in time_strings
        parts = split(time_string, ":")
        hours = parse(Int, parts[1])
        minutes = parse(Int, parts[2])
        seconds = parse(Int, parts[3])
        total_seconds = hours * 3600 + minutes * 60 + seconds
        push!(seconds_array, total_seconds)
    end
    return seconds_array
end

function find_changes(arr)
    changes = []
    for i in 2:length(arr)
        if arr[i] != arr[i-1]
            push!(changes, (i, arr[i-1], arr[i]))
        end
    end
    return changes
end

t_cyc = Float64.(time_to_seconds(string.(data[:,1])))
plot(t_cyc)

V_cyc = Float64.(data[:,2])
I_cyc = Float64.(data[:,3])
I_ind = [i for i in 2:length(I_cyc) if abs(I_cyc[i] - I_cyc[i-1]) > 0.0001]
I_ind./3600

plot(t_cyc/3600, V_cyc, xlabel="Time [hour]", ylabel="Voltage [V]", label="Voltage", 
    color=:blue, linewidth=2, box=:on, legend=false,
    y_guidefontcolor=:blue,
    y_foreground_color_axis=:blue,
    y_foreground_color_text=:blue,dpi=200)
plot!(twinx(), t_cyc/3600, I_cyc/5e-3, ylabel="Current [C-rate]", label="Current", 
        color=:red, linewidth=2, y_guidefontcolor=:red,
        y_foreground_color_axis=:red,
        y_foreground_color_text=:red,
        right_margin=0.3cm,dpi=200,
        legend=false)

# savefig("IV.png")

ind_40_charge = I_ind[4]:I_ind[5]-1
plot(t_cyc[ind_40_charge]/3600, V_cyc[ind_40_charge], xlabel="Time [hour]", ylabel="Voltage [V]", label="Voltage", 
        color=:blue, linewidth=2, box=:on,
        size=(300, 600))

ind_40_discharge = I_ind[6]:length(t_cyc)
plot!(t_cyc[ind_40_discharge]/3600, V_cyc[ind_40_discharge], xlabel="Time [hour]", ylabel="Voltage [V]", label="Voltage", 
        color=:blue, linewidth=2, box=:on,
        size=(300, 600))

# Calibrate stoichometry at 1/40 C 
# Figure 3A DOI:10.1149/1945-7111/ac7e77 for intercalation
# reduction peaks 0.203V and 0.074V corresponding to x = 1/12 and 0.5 respectively

in40V = V_cyc[ind_40_charge]
de40V = V_cyc[ind_40_discharge]

in40Q = t_cyc[ind_40_charge]/3600 * 5 / 20 #[mAh]
in40Q = in40Q .- in40Q[1]
de40Q = t_cyc[ind_40_discharge]/3600 * 5 / 20 #[mAh]

in30t = (t_cyc[ind_40_charge] .- t_cyc[ind_40_charge[1]])/1.07
de30t = t_cyc[ind_40_discharge] .- t_cyc[ind_40_discharge[1]]

using Plots; gr()
plot(in40Q, in40V, ylims = (0,0.4), label = "charge", ylabel="Voltage [V]")
plot!(in40Q[end] .- (de40Q .- de40Q[1]) , de40V, label = "discharge",xlabel="Capacity [mAh]")
plot!(in40Q[end] .- (de40Q .- de40Q[1])*1.07 , de40V, label = "discharge at\n 100% Faraday efficiency", dpi = 150)

savefig("./Exp/20.png")

# savefig("./paper4/Exp/30.png")

muh30_de = LinearInterpolation(de40V,
    de30t; extrapolate = true)
muh30_in = LinearInterpolation(in40V,
    in30t; extrapolate = true)

# using Statistics
# dQ = in40Q[500] .- in40Q[499]
# dQdV = dQ ./ (in40V[2:end] .- in40V[1:end-1])
# window_size = 1000
# moving_avg = [mean(in40V[i:min(i+window_size-1, end)]) for i in 1:length(in40V)-window_size+1]
# plot(in40Q[Int64.(0.5e5:1.9e5)], moving_avg[Int64.(0.5e5:1.9e5)], lw=3)

x50 = in40Q[findall(x -> abs(x .- 0.08) < 0.001, 
    in40V)]
x50 = x50[1]
x00 = in40Q[findall(x -> abs(x .- 0.25) < 0.001, 
    in40V)]
x00 = x00[1]

LiCx40charge = (in40Q .- x00)*(0.5-1/12)/(x50 .- x00) .+ 1/12
scale_charge = (0.5 - 1/12) / (x50 .- x00)

LiCx40discharge = LiCx40charge[end] .- scale_charge[1] * (de40Q .- de40Q[1]) 

plot(LiCx40charge, in40V, label = "intercalation", ylims = (0,0.4))
plot!(LiCx40discharge, de40V, label = "de-intercalation")

using DataInterpolations
muh1 = LinearInterpolation(reverse(de40V)[2:100:end-1],reverse(LiCx40discharge)[2:100:end-1], extrapolate=true)
plot(cavg2, -muh1.(cavg2),label = "de-intercalation", ylims = (-0.4,0.0), dpi=200)

muh = LinearInterpolation(in40V[1:100:end],LiCx40charge[1:100:end], extrapolate=true)
plot!(cavg1, -muh.(cavg1),label = "intercalation")
