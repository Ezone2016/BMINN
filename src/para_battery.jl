using SparseArrays

R = 8.314                     # Gas ant [J.K-1.mol-1]
F = 96487                     # Faraday ant [C.mol-1]
T_ref = 25 + 273.15           # Reference temperature [K]
tplus = 0.2594                 # Transference number
Rc     = 20e-4             # Current collector resistance [ohm.m2]

As = 0.1027   # Electrode active surface area [m2] (assumed equal for anode and cathode)
Ce_init = 1000   #Electrode initial concentration [mol/m^3]
Rs1 = 5.86e-6                 # Anode solid particles' radius [m]
Rs3 =  5.22e-6                 # Cathode solid particles' radius [m]
L1 = 85.2e-6                  # Anode thickness [m]
L2 = 12e-6                    # Seperator thickness [m]
L3 = 75.6e-6                  # Cathode thickness [m]

eps1s = 0.75                # anode
b1 = 1.5
eps2e = 0.47               # separator
b2 = 1.5
eps3s = 0.665                # cathode
b3 = 1.5
a1 = 0.5          # anode symmetry
a3 = 0.5          # cathode symmetry

# Max Solid phase concentration [mol.m-3]
cs1_max = 29583
cs3_max = 51765

# Stoichiometry limits [-]
x1_soc0 = 0.0279              # at 0# soc in anode
x1_soc1 = 0.9014              # at 100# soc in anode
y3_soc0 = 0.9084              # at 0# soc in cathode
y3_soc1 = 0.2661          	# at 100# soc in cathode

# Diffusion coefficient of Li in active material [m2.s-1]
Ds1 = 3.3e-14     # Anode diffusion coeff @ ref temperature [m2.s-1]
Ds3 = 4e-15     # Cathode diffusion coeff @ ref temperature [m2.s-1]

Ea_Ds1 = 0         # Activation energy for Arrhenius' law [J/mol]
Ea_Ds3 = 0         # Activation energy for Arrhenius' law [J/mol]

# Reaction rate ant [m2.5 mol-0.5 s-1]
k1_ref = 6.716047e-12
k3_ref = 3.54458e-11

Ea_k1 = 35000      # Activation energy [J/mol]
Ea_k3 = 17800      # Activation energy [J/mol]

#  ce_avg = 1.0e3    # Average electrolyte concentration
L = L1 + L2 + L3;
as1 = (3*eps1s)/Rs1; # anode
as3 = (3*eps3s)/Rs3; # cathode
Kc = 0.18*eps3s^b3;  # Cathode conductivity [S.m-1]
Ka = 215*eps1s^b1;  # Anode conductivity [S.m-1]
coefL1 = 4*(1-eps1s)^(b1-1)/L1^2;
coefL2 = 4*eps2e^(b2-1)/L2^2;
coefL3 = 4*(1-eps3s)^(b3-1)/L3^2;
cons = 2*R*(1-tplus)/F

M = 2*Np
# Computational coordinates [Chebyshev nodes]

# Computing the Chebyshev differentiation matrices
xm,DM = chebdif(2*Np+1,2);
xr = xm[1:Np+1];
DM10 = DM[1,:,1];
DM2 = DM[1:Np,:,2];

P = Matrix{Float64}(I,Np,Np)[Np:-1:1,:]; # Permutation matrix, backward-identity

# Modified differentation matrices accounting for the solution symmetry
DN2 = DM2[:,1:Np] - DM2[:,Np+2:M+1]*P;
DN1 = DM10[1:Np] - P*DM10[Np+2:M+1];

A = DN2[2:Np,2:Np] + (DN2[2:Np,1]*DN1[2:Np]')/(1-DN1[1])
B = DN2[2:Np,1]/(1-DN1[1])
C = DN1[2:Np]/(1-DN1[1])'
D = 1/(1-DN1[1])

# Anode solid diffusion du/dt = Au + Bj & c = Cu + Dj
#  A1 = (1/Rs1^2)*A;
Di1 =(1/Rs1^2)*eigvals(A)[end:-1:1];
V1 = eigvecs(A)[:,end:-1:1];
V2 = inv(V1);
B2 = inv(V1)*B;
C1 = [C'; Diagonal(1 ./ (xm[2:Np]))];
D1 = [D; zeros(Np-1,1)];

#  A3 = (1/Rs3^2)*A;
Di3 = (1/Rs3^2)*eigvals(A)[end:-1:1];
V3 = eigvecs(A)[:,end:-1:1];
V4 = inv(V3);
C3 = [C'; Diagonal(1 ./ (xm[2:Np]))];
D3 = [D; zeros(Np-1,1)];

placeh,wm0 = clencurt(2*Np);
wn = [(wm0[1:Np]' + wm0[Np+2:2*Np+1]'*P) wm0[Np+1] ];
placeh,Wc0 = clencurt(Nc);  Ic = reverse(cumsummat(Nc+1),dims=1);
placeh,Ws0 = clencurt(Ns);  Is = reverse(cumsummat(Ns+1),dims=1);
placeh,Wa0 = clencurt(Na);  Ia = reverse(cumsummat(Na+1),dims=1);

Wc = Wc0
Ws = Ws0
Wa = Wa0
WcT = Wc'*L3/2;
WaT = Wa'*L1/2;
WsT = Ws'*L2/2;

x1,DeM01 = chebdif(Na+1,1);
DeM01 = DeM01 *2*(1-eps1s)^b1/L1;
x2,DeK01 = chebdif(Ns+1,1);
DeK01 = DeK01*2*eps2e^b2/L2;
x3,DeQ01 = chebdif(Nc+1,1);
DeQ01 = DeQ01*2*(1-eps3s)^b3/L3;

Ae = diagm(0 => [DeQ01[1,1]; -DeK01[1,1]+DeQ01[Nc+1,Nc+1]; -DeM01[1,1]+
       DeK01[Ns+1,Ns+1]; DeM01[Na+1,Na+1]]);
Ae[diagind(Ae, 1)] = [DeQ01[1,Nc+1]; -DeK01[1,Ns+1]; -DeM01[1,Na+1]];
Ae[diagind(Ae, -1)] = [DeQ01[Nc+1,1]; DeK01[Ns+1,1]; DeM01[Na+1,1]];

Be = [-DeQ01[1,2:Nc]' zeros(1, Ns+Na-2);
    -DeQ01[Nc+1,2:Nc]' DeK01[1,2:Ns]' zeros(1, Na-1);
    zeros(1, Nc-1) -DeK01[Ns+1,2:Ns]' DeM01[1,2:Na]';
    zeros(1, Ns+Nc-2) -DeM01[Na+1,2:Na]'];
Me = Ae\Be;

placeh,DeM1 = chebdif(Na+1,1); placeh,DeK1 = chebdif(Ns+1,1); placeh,DeQ1 = chebdif(Nc+1,1);
placeh,DeM2 = chebdif(Na+1,2); placeh,DeK2 = chebdif(Ns+1,2); placeh,DeQ2 = chebdif(Nc+1,2);
DM1 = DeM1[:,:,1]*2/L1;
DK1 = DeK1[:,:,1]*2/L2;
DQ1 = DeQ1[:,:,1]*2/L3;
DeM2 = DeM2[:,:,2]; DeK2 = DeK2[:,:,2]; DeQ2 = DeQ2[:,:,2];
DD3210 = [[DeQ2[2:Nc,2:Nc] zeros(Nc-1,Na+Ns-2)]; zeros(Na+Ns-2,Na+Ns+Nc-3)];
DD3210 += [[[zeros(Nc-1,Nc+Ns-2); zeros(Ns-1,Nc-1) DeK2[2:Ns,2:Ns]] zeros(Nc+Ns-2,Na-1)]; zeros(Na-1,Na+Ns+Nc-3)];
DD3210 += [zeros(Nc+Ns-2,Na+Ns+Nc-3); zeros(Na-1,Nc+Ns-2) DeM2[2:Na,2:Na]];

DD321 = DD3210 + hcat([DeQ2[2:Nc,1]; zeros(Na+Ns-2,1)],
               [DeQ2[2:Nc,Nc+1]; DeK2[2:Ns,1]; zeros(Na-1,1)],
               [zeros(Nc-1,1); DeK2[2:Ns,Ns+1]; DeM2[2:Na,1]],
               [zeros(Nc+Ns-2,1); DeM2[2:Na,Na+1]])*Me;

D3210 = [[DeQ1[2:Nc,2:Nc] zeros(Nc-1,Na+Ns-2)]; zeros(Na+Ns-2,Na+Ns+Nc-3)];
D3210 += [[[zeros(Nc-1,Nc+Ns-2); zeros(Ns-1,Nc-1) DeK1[2:Ns,2:Ns]] zeros(Nc+Ns-2,Na-1)]; zeros(Na-1,Na+Ns+Nc-3)];
D3210 += [zeros(Nc+Ns-2,Na+Ns+Nc-3); zeros(Na-1,Nc+Ns-2) DeM1[2:Na,2:Na]];

D321 = D3210 + hcat([DeQ1[2:Nc,1]; zeros(Na+Ns-2,1)],
               [DeQ1[2:Nc,Nc+1]; DeK1[2:Ns,1]; zeros(Na-1,1)],
               [zeros(Nc-1,1); DeK1[2:Ns,Ns+1]; DeM1[2:Na,1]],
               [zeros(Nc+Ns-2,1); DeM1[2:Na,Na+1]])*Me;

# phase-changing materials

kB = 1.38064852e-23;
NA = 6.02214076e23;
eV = 1.602e-19;
para_Rp = 1e-7; #particle radius [m]
para_h = 0.4/3.3 * para_Rp; #particle thickness [m]
para_cm = 1.379e28 #max concentration [m-3]
T = 298.2033; #[K]
k = 4e-3; #gradient energy coefficient [J/m] 4e-7
ome = 0.115*eV #enthalpy of mixing [J]
# ome = -0.0514*eV; #attraction
# ome = 2.57e-2*eV; #repulsive. high temperature
# para_D0 = 1.05e-15;
# para_D0 = 3.3e-14;
para_D0 = 1.0e-12;
para_alpha = 0.5;
para_I0 = para_Rp/para_cm/eV/para_D0*40; #reference current density [A/m2]
para_V0 = 0.13; #reference potential
para_ce0 = NA/para_cm #reference electrolyte concentration [m-3]
para_Omega = ome/kB/T
# para_Omega = 2
# para_Kappa = k/(para_cm*kB*T);
# para_Kappa = para_Kappa * 1e6;
para_Kappa = 3.13e9 * eV/(para_Rp^2*para_cm*kB*T)
para_beta = 0;
para_kT = kB * T;
para_eV = eV;
para_RT = kB*NA*T;
para_M = para_D0/para_kT

para_Tc = ome/kB/2 #Curie temperature [K]
para_lambda = sqrt(k/para_cm/ome); #interfacial thickness [m]
para_gamma_b = sqrt(k*ome*para_cm); #interfacial tension [J/m2]

# Phase seperating materials
SU = (x,xc,del) -> @. 0.5*(tanh((x-xc)/del)+1)
SD = (x,xc,del) -> @. 0.5*(-tanh((x-xc)/del)+1)
dr = 1/Npa
r = 0:dr:1
rr = r

MM = diagm(3/4*ones(Npa+1))
MM[diagind(MM, 1)] = 1/8*ones(Npa)
MM[diagind(MM, -1)] = 1/8*ones(Npa)
MM[2,1] = 1/4;  MM[end-1,end] = 1/4;

VV = diagm(4*pi*[dr^3/24; r[2:Npa].^2*dr .+ dr^3/12; r[Npa+1]^3/3-(r[Npa+1]-dr/2)^3/3])
# V = diagm(pi*para_h*[dr^2/4; 2*r[2:Npa]*dr; r[Npa+1]*dr-dr^2/4])
MVinv = sparse(inv(MM*VV))

D2C = diagm(-2*ones(Npa+1))
D2C[diagind(D2C, 1)] = ones(Npa)
D2C[diagind(D2C, -1)] = ones(Npa)
D1C = diagm(zeros(Npa+1))
D1C[diagind(D1C, 1)] = ones(Npa)
D1C[diagind(D1C, -1)] = -ones(Npa)
model_r = diagm(ones(Npa+1))
model_r[diagind(model_r, 1)] = ones(Npa)
model_l = diagm(ones(Npa+1))
model_l[diagind(model_l, -1)] = ones(Npa)
model_Dr = diagm(-ones(Npa+1))
model_Dr[diagind(model_Dr, 1)] = ones(Npa)
model_Dl = diagm(ones(Npa+1))
model_Dl[diagind(model_Dl, -1)] = -ones(Npa)

Ajr = [4*pi*(r[1:Npa].+dr/2).^2; 4*pi*r[Npa+1].^2]
Ajl = [0; 4*pi*(r[2:Npa+1].-dr/2).^2;]

# flakes
# Ajr = [2*pi*para_h*(r[1:Npa].+dr/2); 2*pi*para_h*r[Npa+1]]
# Ajl = [0; 2*pi*para_h*(r[2:Npa+1].-dr/2);]

function get_muh(c)
    c = abs.(c)
    mua = @. (-40*exp(-c/0.015) + 0.075*(tanh((c-0.17)/0.02-1))
                + (tanh((c-0.22)/0.04-1))) * SD(c,0.35,0.05)
    mub = -0.05 .* c.^(-0.85)
    muc = 10 .* SU(c,1,0.045)
    mud = @. 6.12 .* (0.4-c.^0.98) .* SD(c,0.49,0.045) .* SU(c,0.35,0.05)
    mue = @. (1.36 .* (0.74.-c) .+ 1.26) .* SU(c,0.5,0.02)
    muh = @. 0.18 + mua + mub + muc + mud + mue

    # muh = @. log(c/(1-c)) + para_Omega*(1-2*c)
end

# stoich = range(0,1,step=1e-2)
# plot(get_muh.(stoich))

# function get_muh(x1)
#     # muh  = @. 1.97938*2.7182818284*exp(-39.3631*x1) + 0.2482 -
#     #     0.0909*tanh(29.8538*(x1 - 0.1234)) - 0.04478*tanh(14.9159*(x1 - 0.2769)) -
#     #     0.0205*tanh(30.4444*(x1 - 0.6103))
#
#     # muh  = @. 0.7222 + 0.1387*x1 + 0.0290*x1.^(1/2) - 0.0172./x1 +
#     #  0.0019./(x1.^(1.5)) + 0.2808*exp(0.90-15*x1) -
#     #  0.7984*exp(0.4465*x1-0.4108)
#
#     muh  = @. 0.5*8.3145*T/F*log10((1-x1)/x1) + 3.5392*exp(-50.381*x1) - 0.13472 +
#         98.941*tanh(4.1465*(x1-0.33873)) + 102.43*tanh(3.8043*(x1-0.31895)) -
#         0.19988*tanh(22.515*(x1-0.11667)) - 200.87*tanh(3.9781*(x1-0.32969))
#
#     return .-muh
# end

function xc2xp(xc,domain::String)
   if  cmp(domain,"r1") .== 0
       xp = Rs1*xc
   elseif  cmp(domain,"r3") .== 0
       xp = Rs3*xc
   else
       error("ERROR: Domain must be r1 for anode | r3 for cathode")
   end
end

function xe2xep(xe,domain::String,)
   if  cmp(domain,"r1") .== 0
       xep = reverse(0.5*(xe .+1)*L1)
   elseif  cmp(domain,"r2") .== 0
       xep = reverse(0.5*(xe .+1)*L2) .+L1
   elseif  cmp(domain,"r3") .== 0
       xep = reverse(0.5*(xe .+1)*L3) .+L1 .+L2
   else
       error("ERROR: Domain must be r1 for anode, r2 for seperator | r3 for cathode")
   end
end

function get_OCP(x1,x3,T)

   # ANODE (Graphite - LiC6)
   # Open-circuit potential at reference temperature 25dC
   # V1_ref = similar(x1);
   @. V1_ref   = 0.7222 + 0.1387*x1 + 0.0290*x1^(1/2) - 0.0172/x1 +
       0.0019/(x1^(1.5)) + 0.2808*exp(0.90-15*x1) -
       0.7984*exp(0.4465*x1-0.4108);

   # Entropy coefficient
   c1_num = [-16515.05308;38379.18127;-37147.89470;19329.75490;
       -5812.27813;1004.91101;-91.79326;3.29927;0.00527];
   c1_den = [165705.85970;-385821.16070;374577.31520;-195881.64880;
       59431.30001;-10481.80419;1017.23480;-48.09287;1];
   dV1_dT_num = c1_num[1];
   dV1_dT_den = c1_den[1];
   for i = 2:7
       dV1_dT_num = dV1_dT_num.*x1 .+ c1_num[i];
       dV1_dT_den = dV1_dT_den.*x1 .+ c1_den[i];
   end
   dV1dT = 1e-3*dV1_dT_num./dV1_dT_den;

   # Open-circuit potential at temperature T
   V1 = V1_ref .+ (ones(size(x1,1),1)*T.-298.15).*dV1dT;

   # CATHODE (Cobalt oxide - LiCoO2)
   # Open-circuit potential at reference temperature 25dC
   V3_ref = zeros(size(x3));
   @. V3_ref = ( -4.656 + 88.669*x3^2 - 401.119*x3^4 + 342.909*x3^6 -
          462.471*x3^8 + 433.434*x3^10)/
       ( -1 + 18.933*x3^2 - 79.532*x3^4 + 37.311*x3^6 -
          73.083*x3^8 + 95.96*x3^10);

   # Entropy coefficient
   c3_num = [0;0.61154;-1.36455;0.92837;-0.19952];
   c3_den = [3.04876;-9.82431;11.47636;-5.66148;1];
   dV3dT_num = c3_num[1];
   dV3dT_den = c3_den[1];
   for i = 2:5
       dV3dT_num = dV3dT_num.*x3 .+ c3_num[i];
       dV3dT_den = dV3dT_den.*x3 .+ c3_den[i];
   end
   dV3dT = 1e-3*dV3dT_num./dV3dT_den;

   # Open-circuit potential at temperature T
   V3 = V3_ref .+ (ones(size(x3,1),1)*T.-298.15).*dV3dT;
   return V1, dV1dT, V3, dV3dT
end
