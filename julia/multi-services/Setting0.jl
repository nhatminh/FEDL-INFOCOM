using Distributions
using HDF5

NumbDevs = 10 #1, 5, 10, 15, 20
NUM_SIM = 1
### PROGRAM SETTING ###
Numb_Services = 3  #Number of simulations
REUSED_TRAFFIC = true
READ_RESULT = true       # Read results from result files

DEBUG = 0 #LEVEL 0, 1, 2, 3
HETEROGENEOUS = 2  # 0: homogeneous, 1: heterogeneous, 2: reused  #Change setting paremeter
# if (REUSED_TRAFFIC == false)
#     HETEROGENEOUS = 1
# end

# # ByDataset = false
# NUMERICAL_RS = false      # Result for Section V: NUMERICAL RESULTS in the paper (50 devs)
# # NUMERICAL_RS = false   # Result for Section IV: Closed-Form solution in the paper (5 devs)

# kaps = [5e-5, 8e-5, 9e-5, 1e-4, 1.3e-4, 1.7e-4, 2e-4, 2.3e-4, 2.7e-4, 3e-4, 3.5e-4, 4e-4, 4.5e-4,
# 5e-4, 5.5e-4, 6e-4, 6.5e-4, 7e-4, 7.5e-4, 8e-4, 8.5e-4, 9e-4, 9.5e-4,
# 1e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3, 4e-3, 6e-3, 8e-3, 1e-2, 3e-2, 5e-2, 7e-2,
# 1e-1, 0.3, 0.5, 0.7, 0.85, 1.,1.5, 2.,2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8., 9., 1e1, 3e1, 5e1, 7e1, 1e2, 3e2, 5e2, 7e2, 1e3]

# kaps = [1e-1]
kaps = [0.2]


# if(NUMERICAL_RS)
#     kaps = [1e-4, 5e-4, 1e-3, 3e-3, 5e-3, 7e-3, 1e-2, 3e-2, 5e-2, 7e-2, 1e-1, 0.5, 1., 2., 3., 4., 5., 6., 7., 8., 9.,
#     1e1, 2e1, 3e1, 5e1, 7e1, 1e2, 3e2, 5e2, 7e2, 1e3, 3e3, 5e3, 7e3, 1e4]
#     NumbDevs = 50 #1, 5, 10, 15, 20
#     HETEROGENEOUS = 0
#     # D_ratios = collect(0.1:0.1:0.9)
#     D_ratios =[1.,0.5, 0.2, 0.01, 0.001]  #D_min/D_max
#     Numb_D = size(D_ratios)[1]
#     # Dis_ratios = collect(0.1:0.1:0.8) # => h_ratios
#     Dis_ratios = [1., 0.5, 0.1, 0.01, 0.001]
#     Numb_Dis = size(Dis_ratios)[1]
# end

if(READ_RESULT)
    REUSED_TRAFFIC = true
    HETEROGENEOUS = 2
end

Numb_Iteration = 300 #1000 200
stop_epsilon1 = 1e-4 #5e-5
stop_epsilon2 = 1e-5 #5e-5
stop_epsilon3 = 2e-3

# RHO1 = 8e-1  #5e-1
# RHO2 = 8e-1
RHO1 = 450  #based on residual, we can control rho (high, faster convergence but can be suboptimal)
RHO2 = 200 #based on residual, we can control rho

if NumbDevs > 10
    RHO1 = 1
    RHO2 = 1
    stop_epsilon1 = 5e-6
    stop_epsilon2 = 5e-6
end

NU = RHO1 *(Numb_Services/(2-0.8) - 1) + 0.1 ## rho *(Numbservice/(2-alpha) - 1) (0< alpha <2)
#(slower NU faster convergence decreasing increase alpha -> can be jumping around solution)
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
BW  = 10e6     #Mhz -> Increase BW -> Increase Theta
w_min = 0.01

Dist_min = 2  #2m
Dist_max = 50 #50m
Dist_avg   = (Dist_max + Dist_min)/2.
println("Dist_avg:",Dist_avg)  # 26m

### COMPUTING PARAMS ###
function save_setting(f_max, C_s, V_s)
    h5open("setting.h5", "w") do file
        # write(file,"f_min", f_min)
        write(file,"f_max", f_max)
        write(file,"C_s", C_s)
        write(file,"V_s", V_s)
    end
end

function read_setting()
    h5open("setting.h5", "r") do file
        # f_min =read(file,"f_min")
        f_max =read(file,"f_max")
        C_s = read(file,"C_s")
        V_s = read(file,"V_s")
        return f_max, C_s, V_s
    end
end

if (HETEROGENEOUS == 0) # Homogeneous
    cpu_max  = 2.  #GHz, cyles/1sec
    # f_min = cpu_min*ones(NumbDevs)*1e9  #Hz
    f_max = cpu_max*ones(NumbDevs)*1e9  #Hz
    C_s   = 20*ones(Numb_Services)   #cycles/bits
    V_s   =  5*ones(Numb_Services)*8e3
    save_setting(f_max,C_s,V_s)

elseif (HETEROGENEOUS == 1) # Heterogeneous
    # cpu_min1  = .1  #GHz, cyles/1sec
    # cpu_min2  = .2  #GHz, cyles/1sec
    cpu_max1  = 1.  #GHz, cyles/1sec
    cpu_max2  = 2.  #GHz, cyles/1sec
    # f_min = rand(Uniform(cpu_min1,cpu_min2),NumbDevs)*1e9  #Hz
    f_max = rand(Uniform(cpu_max1,cpu_max2),NumbDevs)*1e9  #Hz
    # C_s   = rand(10:30,Numb_Services)   #cycles/bits
    C_s   = [30, 50, 70]      #cycles/bits
    V_s   = [300, 400, 500] *8e3   #bits
    save_setting(f_max,C_s,V_s)

elseif (HETEROGENEOUS == 2) # Reused params
    f_max, C_s, V_s = read_setting()
    C_s   = [30, 50, 70]      #cycles/bits
    V_s   = [300, 400, 500] *8e3   #bits
end
cpu_min  = .3 #GHz, cyles/1sec
f_min = cpu_min/Numb_Services*ones(Numb_Services)*1e9  #Hz

println(C_s)
println(f_min)
println(f_max)

alpha    = 2e-28

Numb_kaps = size(kaps)[1]

#FEDL convergence analysis parameters
L=1 #L-Lipschitz
beta=0.5 #beta strongly convex
rho0 = L/beta
# gamma=0.5
c=1
C=c*rho0
thetas = [0.07,0.06,0.05]
K_l, K_g_const = zeros(Numb_Services), zeros(Numb_Services)
A1s,A2s,B2s = zeros(Numb_Services), zeros(Numb_Services), zeros(Numb_Services)
gamma = 1. #0.5

for s =1:Numb_Services
    K_l[s] = 2/gamma*log(C/thetas[s])
    A1s[s] = 1
    A2s[s] =-rho0*(1+thetas[s])^2
    K_g_const[s] = ((1-thetas[s])/rho0)^2 - thetas[s]*(1+thetas[s])
    B2s[s] = 2*rho0*K_g_const[s]
    # K_g= 2*eta*rho0 *( ((1-thetas[s])/rho0)^2 - thetas[s]*(1+thetas[s]) - (1+thetas[s])^2*eta/2 )
end
println("K_l:",K_l)
