using DelimitedFiles, Dates
using Plots.PlotMeasures
using DataInterpolations

# Read the data from the file
data = readdlm("./exp/cell1_20.txt")[2:end-1,1:3]
# data = readdlm("../exp/cell1_20.txt")[2:end-1,1:3]
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
V_cyc = Float64.(data[:,2])
I_cyc = Float64.(data[:,3])
I_ind = [i for i in 2:length(I_cyc) if abs(I_cyc[i] - I_cyc[i-1]) > 0.0001]

ind_40_charge = 1:I_ind[1]
ind_40_discharge = I_ind[2]:I_ind[3]

# Calibrate stoichometry at 1/40 C 
# Figure 3A DOI:10.1149/1945-7111/ac7e77 for intercalation
# reduction peaks 0.203V and 0.074V corresponding to x = 1/12 and 0.5 respectively

in40V = V_cyc[ind_40_charge]
de40V = V_cyc[ind_40_discharge]

in40Q = t_cyc[ind_40_charge]/3600 * 5 / 40 #[mAh]
de40Q = t_cyc[ind_40_discharge]/3600 * 5 / 40 #[mAh]

plot(in40Q, in40V, ylims = (0,0.4), label = "charge", title = "1/60 C",ylabel="Voltage [V]")
plot!(in40Q[end] .- (de40Q .- de40Q[1]) , de40V, label = "discharge",xlabel="Q [mAh]")
plot!(in40Q[end] .- (de40Q .- de40Q[1])*1.08 , de40V, label = "adjusted discharge", dpi = 150)

in30t = (t_cyc[ind_40_charge] .- t_cyc[ind_40_charge[1]])/1.08
de30t = t_cyc[ind_40_discharge] .- t_cyc[ind_40_discharge[1]]

x50 = in40Q[findall(x -> abs(x .- 0.08) < 0.001, 
    in40V)]
x50 = x50[1]
x00 = in40Q[findall(x -> abs(x .- 0.22) < 0.001, 
    in40V)]
x00 = x00[1]

LiCx40charge = (in40Q .- x00)*(0.5-1/12)/(x50 .- x00) .+ 1/12
LiCx40charge = LiCx40charge/1.08 
scale_charge = (0.5 - 1/12) / (x50 .- x00)

LiCx40discharge = LiCx40charge[end] .- scale_charge[1] * (de40Q .- de40Q[1]) 


V_ind = [i for i in 1:length(in40V) if abs(in40V[i]) > 0.4631]
in30t = in30t[V_ind[end]:end]
in30t = in30t .- in30t[1]
in40V = in40V[V_ind[end]:end]
LiCx40charge = LiCx40charge[V_ind[end]:end]

V_ind = [i for i in 1:length(de40V) if abs(de40V[i]) > 0.4631]
de30t = de30t[1:V_ind[1]]
de30t = de30t .- de30t[1]
de40V = de40V[1:V_ind[1]]
LiCx40discharge = LiCx40discharge[1:V_ind[1]]

plot(in30t, in40V, ylims = (0,0.4631))
plot!(de30t, reverse(de40V), ylims = (0,0.4))

muh60_de = LinearInterpolation(de40V[2:100:end-1],
    de30t[2:100:end-1])
muh60_in = LinearInterpolation(in40V[2:100:end-1],
    in30t[2:100:end-1])

muh60_de = LinearInterpolation(de40V[2:100:end-1],
    de30t[2:100:end-1]; extrapolate = true)
muh60_in = LinearInterpolation(in40V[2:100:end-1],
    in30t[2:100:end-1]; extrapolate = true)
    
plot(LiCx40charge .+ 0.035, in40V, label = "intercalation", ylims = (0,0.4631),
    xticks = ([0.08, 0.87], ["0.08", "0.87"]), )
plot!(LiCx40discharge .+ 0.035, de40V, label = "de-intercalation")
# plot!(c, -get_mu(c,θθ.p1))
muh60_deQ = LinearInterpolation(reverse(de40V)[2:100:end-1],reverse(LiCx40discharge)[2:100:end-1] .+ 0.035, extrapolate=true)
muh60_inQ = LinearInterpolation(in40V[1:100:end],LiCx40charge[1:100:end] .+ 0.035, extrapolate=true)

data = readdlm("./exp/cell1_20.txt")[2:end-1,1:3]
# data = readdlm("../exp/cell1_20.txt")[2:end-1,1:3]

t_cyc = Float64.(time_to_seconds(string.(data[:,1])))
V_cyc = Float64.(data[:,2])
I_cyc = Float64.(data[:,3])
I_ind = [i for i in 2:length(I_cyc) if abs(I_cyc[i] - I_cyc[i-1]) > 0.0001]

ind_40_charge = I_ind[4]:I_ind[5]-1
ind_40_discharge = I_ind[6]:length(t_cyc)

# Calibrate stoichometry at 1/40 C 
# Figure 3A DOI:10.1149/1945-7111/ac7e77 for intercalation
# reduction peaks 0.203V and 0.074V corresponding to x = 1/12 and 0.5 respectively

in40V = V_cyc[ind_40_charge]
de40V = V_cyc[ind_40_discharge]

in40Q = t_cyc[ind_40_charge]/3600 * 5 / 40 #[mAh]
in40Q = in40Q .- in40Q[1]
de40Q = t_cyc[ind_40_discharge]/3600 * 5 / 40 #[mAh]

in30t = (t_cyc[ind_40_charge] .- t_cyc[ind_40_charge[1]])/1.07
de30t = t_cyc[ind_40_discharge] .- t_cyc[ind_40_discharge[1]]

V_ind = [i for i in 1:length(in40V) if abs(in40V[i]) > 0.4631]
in30t = in30t[V_ind[end]:end]
in30t = in30t .- in30t[1]
in40V = in40V[V_ind[end]:end]

V_ind = [i for i in 1:length(de40V) if abs(de40V[i]) > 0.4631]
de30t = de30t[1:V_ind[1]]
de30t = de30t .- de30t[1]
de40V = de40V[1:V_ind[1]]

muh30_de = LinearInterpolation(de40V[2:100:end-1],
    de30t[2:100:end-1]; extrapolate = true)
muh30_in = LinearInterpolation(in40V[2:100:end-1],
    in30t[2:100:end-1]; extrapolate = true)

ocvfit = [0.08385 -0.38692;
0.08556 -0.37398;
0.08740 -0.36354;
0.08924 -0.35316;
0.09108 -0.34279;
0.09292 -0.33237;
0.09475 -0.32187;
0.09659 -0.31140;
0.09843 -0.30092;
0.10027 -0.29048;
0.10211 -0.28011;
0.10395 -0.26971;
0.10579 -0.25927;
0.10763 -0.24877;
0.10947 -0.23831;
0.11131 -0.22800;
0.11315 -0.21803;
0.11533 -0.20652;
0.12029 -0.19411;
0.13124 -0.19275;
0.14442 -0.20301;
0.14840 -0.21391;
0.15738 -0.22633;
0.17383 -0.21983;
0.21245 -0.16971;
0.22012 -0.16022;
0.22777 -0.14967;
0.23787 -0.14091;
0.27020 -0.12584;
0.29053 -0.12334;
0.31089 -0.12280;
0.35168 -0.13006;
0.37207 -0.13385;
0.39246 -0.13662;
0.41284 -0.13917;
0.45357 -0.14087;
0.47391 -0.13840;
0.50040 -0.12082;
0.51179 -0.09513;
0.52381 -0.08528;
0.54219 -0.08033;
0.56254 -0.07972;
0.66444 -0.09136;
0.68483 -0.09466;
0.70522 -0.09734;
0.72559 -0.09893;
0.74597 -0.10065;
0.76634 -0.10189;
0.80708 -0.10425;
0.82744 -0.10459;
0.84778 -0.10226;
0.86071 -0.09316;
0.86324 -0.08242;
0.86637 -0.06917;
0.86966 -0.05930;
0.87149 -0.04719;
0.87348 -0.03560;
0.87419 -0.02564;
0.87702 -0.01700]

# ocv_fit = LinearInterpolation(ocvfit[:,2],ocvfit[:,1])
ocv_fit = LinearInterpolation(ocvfit[:,2],ocvfit[:,1]; extrapolate = true)
