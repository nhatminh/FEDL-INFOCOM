
NumbDevs = 10

### COMMUNICATIONS PARAMS ###
N0  = 1e-9
BW  = 1       #Mhz

Dist_min = 2  #2m
Dist_max = 50 #50m

Ptx_Min  = 0.1 # 0.1W -> 20dBm
Ptx_Max  = 1   # 1W   -> 30dBm

### COMPUTING PARAMS ###
cpu_min  = 0.1 #GHz, cyles/1sec
cpu_max  = 2.  #GHz, cyles/1sec
C_n      = 1000   #cycles/bits
alpha    = 2e-28

### LEARNING PARAMS ###
D_min   = 10e6*8 #10 MB, datasize range (-> bits)
D_max   = 20e6*8 #20 MB (-> bits)
W       = 100e3 *8 #100KB, weight params size (-> bits)
kappa   = 1.    #coeff of T_iter

### PROGRAM SETTING ###
REUSED_TRAFFIC = False
DEBUG = 0 #LEVEL 0, 1, 2, 3