NumbDevs = 10

### COMMUNICATIONS PARAMS ###
N0  = 1e-9    #
BW  = 2e6     #Mhz

Dist_min = 2  #2m
Dist_max = 50 #50m

Ptx_Min  = 0.1 # 0.1W -> 20dBm
Ptx_Max  = 1   # 1W   -> 30dBm

### COMPUTING PARAMS ###
cpu_min  = .1e9 #GHz, cyles/1sec
cpu_max  = 2.e9  #GHz, cyles/1sec
C_n      = 100   #cycles/bits
alpha    = 2e-28

### LEARNING PARAMS ###
D_min   = 8* 1e6 #1 MB, datasize range (-> bits)
D_max   = 8* 2e6 #2 MB (-> bits)
S_n     = 10e3 *8 #10KB, weight params size (-> bits), and gradient => 10K nats (1bits/ln2)
kappa   = 1.    #coeff of T_iter

### PROGRAM SETTING ###
REUSED_TRAFFIC = true
DEBUG = 1 #LEVEL 0, 1, 2, 3
