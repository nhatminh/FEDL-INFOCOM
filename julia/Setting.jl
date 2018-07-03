NumbDevs = 10

### PROGRAM SETTING ###
Numb_SIMs = 1  #Number of simulations
REUSED_TRAFFIC = true
DEBUG = 1 #LEVEL 0, 1, 2, 3
HETEROGENEOUS = true  #homogeneous:false, heterogeneous: true

### COMMUNICATIONS PARAMS ###
N0  = 1e-9    #
BW  = 2e6     #Mhz

Dist_min = 2  #2m
Dist_max = 50 #50m

### COMPUTING PARAMS ###
if HETEROGENEOUS
    cpu_min1  = .1  #GHz, cyles/1sec
    cpu_min2  = .2  #GHz, cyles/1sec
    cpu_max1  = 1.  #GHz, cyles/1sec
    cpu_max2  = 2.  #GHz, cyles/1sec
    f_min = rand(cpu_min1:cpu_min2,NumbDevs)*1e9  #Hz
    f_max = rand(cpu_max1:cpu_max2,NumbDevs)*1e9  #Hz
    C_n   = rand(50:150,NumbDevs)   #cycles/bits
else
    cpu_min  = .1 #GHz, cyles/1sec
    cpu_max  = 2.  #GHz, cyles/1sec
    f_min = cpu_min*ones(NumbDevs)*1e9  #Hz
    f_max = cpu_max*ones(NumbDevs)*1e9  #Hz
    C_n      = 100*ones(NumbDevs)   #cycles/bits

end

alpha    = 2e-28

### LEARNING PARAMS ###
D_min   = 8* 1e6 #1 MB, datasize range (-> bits)
D_max   = 8* 2e6 #2 MB (-> bits)
S_n     = 10e3 *8 #10KB, weight params size (-> bits), and gradient => 10K nats (1bits/ln2)
kappa   = 0.001   #coeff of T_iter
