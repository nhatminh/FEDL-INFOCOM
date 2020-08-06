using Distributions
using HDF5

NumbDevs = 50 #1, 5, 10, 15, 20, 50
NUM_SIM = 1  #1, 100
### PROGRAM SETTING ###
Numb_Services = 3  #Number of simulations
REUSED_TRAFFIC = true
READ_RESULT = true       # Read results from result files

DEBUG = 0 #LEVEL 0, 1, 2, 3
HETEROGENEOUS = 2  # 1: heterogeneous, 2: reused  #Change setting paremeter


kaps = [0.5]
kaps_pareto = [0.05, 0.1, 0.5, 1, 2.5, 5, 7, 10]
kaps_s = [1., 1., 1.]
Numb_kaps = size(kaps)[1]
Numb_kaps1 = size(kaps_pareto)[1]

Numb_Kap_Vals = 8


if(READ_RESULT)
    REUSED_TRAFFIC = true
    HETEROGENEOUS = 2
end

# Numb_Iteration = 300 #1000 200
Numb_Iteration = 200
stop_epsilon1 = 1e-4 #5e-5
stop_epsilon2 = 1e-5 #5e-5
stop_epsilon3 = 2e-3

# RHO1 = 8e-1  #5e-1
# RHO2 = 8e-1
RHO1 = 450  #based on residual, we can control rho (high, faster convergence)
RHO2 = 200 #based on residual, we can control rho

if NumbDevs > 10
    RHO1 = 1000 # faster global convergence but can be fluctuate, slower decision variable convergence
    RHO2 = 10
    stop_epsilon1 = 4e-4 #4e-4
    stop_epsilon2 = 5e-5
    stop_epsilon3 = 1e-2
end

NU = RHO1 *(Numb_Services/(2-0.8) - 1) + 0.1 ## rho *(Numbservice/(2-alpha) - 1) (0< alpha <2)
#(slower NU faster convergence decreasing, increase alpha -> can be jumping around solution)
println("Chosen NU:",NU)

### LEARNING PARAMS ###
D_min   = 8* 10e6   #10 MB, datasize range (-> bits)
D_max   = 8* 20e6  #20 MB (-> bits)
D_avg   = (D_max + D_min)/2.
println("D_avg:",D_avg)  # 15 MB

# D_Total = 8* 20e6 #20 MB (-> bits)
# V_s     = 25e3 #weight params size (-> bits), and gradient => 25K nats (1bits/ln2) -> 4.5KB
# V_s     = 8* 500e3 #500KB
# kappa   = 0.001   #coeff of T_iter

### COMMUNICATIONS PARAMS ###
# Tx_Power = pow(10, 43. / 10) / 1000 #43 dBm -> Watts
Tx_Power = 10. #Watts
Tx_Power_BS = 40 #Watts
# Ptx_Max = 1.
# Ptx_Min = 0.2
N0  = 1e-10    #    -> Decrease BW -> Increase Theta
BW  = 20e6     #Mhz -> Increase BW -> Increase Theta
w_min = 0.01

Dist_min = 2  #2m
Dist_max = 50 #50m
Dist_avg   = (Dist_max + Dist_min)/2.
println("Dist_avg:",Dist_avg)  # 26m

### COMPUTING PARAMS ###
function save_setting(f_max, C_s, V_s, tau_extra, tau_mem)
    h5open("setting.h5", "w") do file
        # write(file,"f_min", f_min)
        write(file,"f_max", f_max)
        write(file,"C_s", C_s)
        write(file,"V_s", V_s)
        write(file,"tau_extra", tau_extra)
        write(file,"tau_mem", tau_mem)
    end
end

function read_setting()
    h5open("setting.h5", "r") do file
        # f_min =read(file,"f_min")
        f_max =read(file,"f_max")
        C_s = read(file,"C_s")
        V_s = read(file,"V_s")
        tau_extra = read(file,"tau_extra")
        tau_mem = read(file,"tau_mem")
        # T_avg = read(file,"T_avg")
        return f_max, C_s, V_s, tau_extra, tau_mem
    end
end

# if (HETEROGENEOUS == 0) # Homogeneous
#     cpu_max  = 2.  #GHz, cyles/1sec
#     # f_min = cpu_min*ones(NumbDevs)*1e9  #Hz
#     f_max = cpu_max*ones(NumbDevs)*1e9  #Hz
#     C_s   = 20*ones(Numb_Services)   #cycles/bits
#     V_s   =  5*ones(Numb_Services)*8e3
#     save_setting(f_max,C_s,V_s)

if (HETEROGENEOUS == 1) # Heterogeneous
    # cpu_min1  = .1  #GHz, cyles/1sec
    # cpu_min2  = .2  #GHz, cyles/1sec
    cpu_max1  = 1.  #GHz, cyles/1sec
    cpu_max2  = 2.  #GHz, cyles/1sec
    # f_min = rand(Uniform(cpu_min1,cpu_min2),NumbDevs)*1e9  #Hz
    f_max = rand(Uniform(cpu_max1,cpu_max2),NumbDevs)*1e9  #Hz
    # C_s   = rand(10:30,Numb_Services)   #cycles/bits
    C_s   = [50, 70, 90]      #cycles/bits #30, 50, 70
    V_s   = [100, 200, 300] *8e3   #bits  #300, 500, 700
    T_avg = [0.5, 1.,2.]
    tau_extra = zeros(Numb_Services,NumbDevs)
    tau_mem = zeros(Numb_Services,NumbDevs)
    for s =1:Numb_Services
        T_avg = zeros(Numb_Services)
        tau_extra[s,:] = rand(Uniform(0.,1.),NumbDevs) #extra time in second for communications
        tau_mem[s,:] = rand(Uniform(0.1,0.5),NumbDevs) #extra time in second for communications
    end
    save_setting(f_max,C_s,V_s,tau_extra, tau_mem)

elseif (HETEROGENEOUS == 2) # Reused params
    f_max, C_s, V_s, tau_extra, tau_mem = read_setting()
    C_s   = [50, 70, 90]      #cycles/bits
    V_s   = [100, 200, 300] *8e3   #bits
    T_avg = [0.5, 1.,2.]
end


# total_proportion_Vs = sum(V_s)
cpu_min  = .3 #GHz, cyles/1sec
f_min = cpu_min/Numb_Services*ones(Numb_Services)*1e9  #Hz

# println(C_s)
# println(f_min)
# println(f_max)

alpha    = 2e-28



#FEDL convergence analysis parameters
L=1 #L-Lipschitz
beta=0.5 #beta strongly convex
rho0 = L/beta
# gamma=0.5
c=1
C=c*rho0
thetas = [0.07,0.06,0.05]
K_l, K_g_const = zeros(Numb_Services), zeros(Numb_Services)
# A1s,A2s,B2s = zeros(Numb_Services), zeros(Numb_Services), zeros(Numb_Services)
As,Bs,Cs,Ds = zeros(Numb_Services), zeros(Numb_Services), zeros(Numb_Services), zeros(Numb_Services)
gamma = 1. #0.5


for s =1:Numb_Services
    K_l[s] = 2/gamma*log(C/thetas[s])
    # As[s] = 1
    As[s] = 1. / 2/rho0 #normalized 2*rho_0
    Bs[s] = (1+thetas[s])^2*rho0^2
    Cs[s] = 2*(thetas[s] -1)^2 - 2*(thetas[s]+1)*thetas[s]*rho0^2
    Ds[s] = rho0^2*(1+3thetas[s])*(1+thetas[s])
end

# println("Bs:",Bs)
println("Cs:",Cs)
# println("Ds:",Ds)
println("K_l:",K_l)
# println("theta:",thetas)
