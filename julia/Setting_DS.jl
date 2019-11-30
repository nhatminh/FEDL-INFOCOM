using Distributions
using HDF5

NumbDevs = 5 #1, 5, 10, 15, 20

### PROGRAM SETTING ###
Numb_SIMs = 1  #Number of simulations
REUSED_TRAFFIC = true
READ_RESULT = true       # Read results from result files

DEBUG = 0 #LEVEL 0, 1, 2, 3
HETEROGENEOUS = 2  # 0: homogeneous, 1: heterogeneous, 2: reused
ByDataset = true
# NUMERICAL_RS = true      # Result for Section V: NUMERICAL RESULTS in the paper (50 devs)
NUMERICAL_RS = false   # Result for Section IV: Closed-Form solution in the paper (5 devs)

kaps = [5e-5, 8e-5, 9e-5, 1e-4, 1.3e-4, 1.7e-4, 2e-4, 2.3e-4, 2.7e-4, 3e-4, 3.5e-4, 4e-4, 4.5e-4,
5e-4, 5.5e-4, 6e-4, 6.5e-4, 7e-4, 7.5e-4, 8e-4, 8.5e-4, 9e-4, 9.5e-4,
1e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3, 4e-3, 6e-3, 8e-3, 1e-2, 3e-2, 5e-2, 7e-2,
1e-1, 0.3, 0.5, 0.7, 0.85, 1.,1.5, 2.,2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8., 9., 1e1, 3e1, 5e1, 7e1, 1e2, 3e2, 5e2, 7e2, 1e3]

# kaps = [5e-5, 8e-5, 9e-5, 1e-4, 1.3e-4, 1.7e-4, 2e-4, 2.3e-4, 2.7e-4, 3e-4, 4e-4,
# 5e-4, 6e-4, 7e-4, 8e-4, 9e-4,
# 1e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3, 3.5e-3, 4e-3, 4.5e-3, 5e-3, 5.5e-3, 6e-3, 6.5e-3, 7e-3, 7.5e-3, 8e-3, 8.5e-3, 9e-3, 9.5e-3, 1e-2, 3e-2, 5e-2, 7e-2,
# 1e-1,0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9,  0.95 ,1.,1.25, 1.5, 1.75, 2., 2.25,2.5,2.75, 3., 3.5, 4., 4.5, 5., 6., 7., 8., 9., 1e1, 3e1, 5e1, 7e1, 1e2, 3e2, 5e2, 7e2, 1e3]

if(NUMERICAL_RS)
    kaps = [1e-4, 5e-4, 1e-3, 3e-3, 5e-3, 7e-3, 1e-2, 3e-2, 5e-2, 7e-2, 1e-1, 0.5, 1., 2., 3., 4., 5., 6., 7., 8., 9.,
    1e1, 2e1, 3e1, 5e1, 7e1, 1e2, 3e2, 5e2, 7e2, 1e3, 3e3, 5e3, 7e3, 1e4]
    NumbDevs = 50 #1, 5, 10, 15, 20
    HETEROGENEOUS = 0
    # D_ratios = collect(0.1:0.1:0.9)
    D_ratios =[1.,0.5, 0.2, 0.01, 0.001]  #D_min/D_max
    Numb_D = size(D_ratios)[1]
    # Dis_ratios = collect(0.1:0.1:0.8) # => h_ratios
    Dis_ratios = [1., 0.5, 0.1, 0.01, 0.001]
    Numb_Dis = size(Dis_ratios)[1]
end

if(READ_RESULT)
    REUSED_TRAFFIC = true
    HETEROGENEOUS = 2
end

### LEARNING PARAMS ###
D_min   = 8* 5e6   #5 MB, datasize range (-> bits)
D_max   = 8* 10e6  #10 MB (-> bits)
D_avg   = (D_max + D_min)/2.
D_DS    = 8* 0.94* 1e6 #~940KB
println("D_avg:",D_avg)  # 15 MB

D_Total = 8* 20e6 #20 MB (-> bits)
# S_n     = 25e3 #weight params size (-> bits), and gradient => 25K nats (1bits/ln2) -> 4.5KB
S_n     = 8* 1.6*1e6/log2(2) #(nats) = (bits/ln2)
# kappa   = 0.001   #coeff of T_iter

### COMMUNICATIONS PARAMS ###
Ptx_Max = 1.
Ptx_Min = 0.2
N0  = 1e-10    #    -> Decrease BW -> Increase Theta
BW  = 20e6     #Mhz -> Increase BW -> Increase Theta

Dist_min = 5  #5m
Dist_max = 40 #40m
Dist_avg   = (Dist_max + Dist_min)/2.
println("Dist_avg:",Dist_avg)  # 26m
### COMPUTING PARAMS ###
function save_setting(f_max, C_n)
    h5open("setting.h5", "w") do file
        # write(file,"f_min", f_min)
        write(file,"f_max", f_max)
        write(file,"C_n", C_n)
    end
end

function read_setting()
    h5open("setting.h5", "r") do file
        # f_min =read(file,"f_min")
        f_max =read(file,"f_max")
        C_n = read(file,"C_n")
        return f_max, C_n
    end
end

if (HETEROGENEOUS == 0) # Homogeneous
    cpu_max  = 2.  #GHz, cyles/1sec
    # f_min = cpu_min*ones(NumbDevs)*1e9  #Hz
    f_max = cpu_max*ones(NumbDevs)*1e9  #Hz
    C_n      = 20*ones(NumbDevs)   #cycles/bits
    save_setting(f_max,C_n)

elseif (HETEROGENEOUS == 1) # Heterogeneous
    # cpu_min1  = .1  #GHz, cyles/1sec
    # cpu_min2  = .2  #GHz, cyles/1sec
    cpu_max1  = 1.  #GHz, cyles/1sec
    cpu_max2  = 2.  #GHz, cyles/1sec
    # f_min = rand(Uniform(cpu_min1,cpu_min2),NumbDevs)*1e9  #Hz
    f_max = rand(Uniform(cpu_max1,cpu_max2),NumbDevs)*1e9  #Hz
    C_n   = rand(2000:2500,NumbDevs)   #cycles/bits
    save_setting(f_max,C_n)

elseif (HETEROGENEOUS == 2) # Reused params
    f_max, C_n = read_setting()
end
cpu_min  = .3 #GHz, cyles/1sec
f_min = cpu_min*ones(NumbDevs)*1e9  #Hz
println("CPU cycles/bits: ",C_n)
println("f_min: ",f_min)
println("f_max: ",f_max)

alpha    = 2e-28

Numb_kaps = size(kaps)[1]
