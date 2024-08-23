begin
    # Constants
    N = 150
    R = 8.314  # Gas constant [J.K-1.mol-1]
    F = 96487  # Faraday constant [C.mol-1]
    kB = 1.38064852e-23  # [m2 kg s-2 K-1]
    NA = 6.0221409e23  # mol-1
    T = 273 + 25  # K
    V_theta = 3.42  # standard potential [V] voltage plateau
    eV = 1.60218e-19
    constant = 0.025693  # R*T/F 25.693mV
    Q = 7e-3 * 3.6e3
    q0 = 0.07
    q1 = 0.86
    avg_pd(mid, f) = trapz(mid, mid .* f)
    c0(c,q0) = @. sqrt(1/2/pi/(v2))*exp(-1/(2*(v2))*(c-q0)^2)
    c = Float32.(collect(range(0.0, stop=1.0, length=N)))
    c_mid = (c[1:end-1] + c[2:end]) / 2
    cc = reshape(c_mid,1,:)
    sc(p) = myNN(cc, p, st)[1][:]
    sc2const = (1 .- c_mid) .* c_mid 
    sc2(p) = sc2const .* myNN2(cc, p, st2)[1][:]
    logc0 = 10*exp.(-c_mid*100) .+ 1;
    logc1 = 10*exp.(-reverse(c_mid)*100) .+ 1;
    dc = c[2] - c[1]

    plot_font = "Computer Modern"
    default(
        fontfamily=plot_font,
        linewidth=3, 
        framestyle=:box, 
        grid=false
    )
end

get_mu(yy,pp) = myNN(reshape(yy,1,:), pp, st)[1][:]
f_muN_de = (us1,pp) -> myNN(reshape(Array(us1)[end,:],1,:), pp.p1, st)[1][:] .- 
    K *eV/(pp.Rp^2*pp.cmax*kB*T) *
    (2/rr[end]*pp.β .+ 2*(us1[end-1,:] .-us1[end,:] .-dr*pp.β)/dr^2)

f_muN_in = (us1,pp) -> myNN(reshape(Array(us1)[end,:],1,:), pp.p1, st)[1][:] .- 
    K *eV/(pp.Rp^2*pp.cmax*kB*T) *
    (-2/rr[end]*pp.β .+ 2*(us1[end-1,:] .-us1[end,:] .+dr*pp.β)/dr^2)
    

D2Cs = sparse(D2C[2:end-1,:]) ./dr^2
D1Cs = sparse(D1C[2:end-1,:]) ./dr
rrs = collect(rr[2:Npa])
D1Cs = 1 ./rrs .* D1Cs
D12C = D2Cs .+ D1Cs

model_r = sparse(model_r)
model_Dr = sparse(model_Dr)
model_l = sparse(model_l)
model_Dl = sparse(model_Dl)
MVinv = inv(MM*VV)
# MVinv[MVinv .< 1e-10] .= 0.0 
# MVinv = sparse(MVinv)
MV = sparse(MM*VV)

model_Dr1 = -model_Dr[1:end-1,:]./ dr
model_r1 = model_r[1:end-1,:]./ 2
model_Dl1 = -model_Dl[2:end,:]./ dr
model_l1 = model_l[2:end,:]./ 2

MVinvAjl = Ajl .* MVinv 
MVinvAjr = Ajr .* MVinv 

f_muN_de = (us1,pp) -> myNN(reshape(Array(us1)[end,:],1,:), pp.p1, st)[1][:] .- 
    pp.K * (2/rr[end]*pp.β .+ 2*(us1[end-1,:] .-us1[end,:] .-dr*pp.β)/dr^2)

f_muN_in = (us1,pp) -> myNN(reshape(Array(us1)[end,:],1,:), pp.p1, st)[1][:] .- 
    pp.K * (-2/rr[end]*pp.β .+ 2*(us1[end-1,:] .-us1[end,:] .+dr*pp.β)/dr^2)

function get_phis_de!(dy, y, pp, t)   
    mu = get_mu(y,pp.p1) .- pp.K*vcat(6*(y[2]-y[1])/dr^2, D12C * y, 2/rr[end]* pp.β 
        + 2*(y[end-1]-y[end]-dr* pp.β)/dr^2)
    Fr = vcat((1 .- model_r1 * y) .* (model_r1 * y) .* 
        (model_Dr1*mu),-Iapp/(3*pp.theta2)*pp.theta1)
    Fl = vcat(0.0, (1 .- model_l1 * y) .* (model_l1 * y) .* (model_Dl1*mu))
    dy .= (MVinvAjl * Fl .- MVinvAjr * Fr) ./ pp.theta1 .*3600
end

function get_phis_in!(dy, y, pp, t)   
    mu = get_mu(y,pp.p1) .- pp.K*vcat(6*(y[2]-y[1])/dr^2, D12C * y, -2/rr[end]*pp.β 
        + 2*(y[end-1]-y[end]+dr*pp.β)/dr^2)
    Fr = vcat((1 .- model_r1 * y) .* (model_r1 * y) .* 
        (model_Dr1*mu),Iapp/(3*pp.theta2)*pp.theta1)
    Fl = vcat(0.0, (1 .- model_l1 * y) .* (model_l1 * y) .* (model_Dl1*mu))
    dy .= (MVinvAjl * Fl .- MVinvAjr * Fr) ./ pp.theta1 .*3600
end

