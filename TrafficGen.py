from Setting import *
import numpy as np
import h5py

def to_watt(x):
    return pow(10, x / 10)

def to_dBW(x):
    return 10*np.log10(x)

### Exponential distribution of channel_gain
def exp_gain(dist):
    g0 = -40 #dB
    g0_W = to_watt(g0)
    d0 = 1 #1 m
    mean = g0_W*pow((d0/dist),4)
    gain_h = np.random.exponential(mean);

    return gain_h #in Watts

def noise_per_gain(gain):  #N0/h
    return N0/gain

# def shanon_capacity():
#     return

def mobile_gen():
    if(REUSED_TRAFFIC):
        return simple_read_data()
    else:
        dist_list = np.random.randint(Dist_min,Dist_max,NumbDevs)
        gain_list = [exp_gain(dist_list[n]) for n in range(NumbDevs)]
        ratios    = [noise_per_gain(gain_list[n]) for n in range(NumbDevs)]
        D_n       = np.random.randint(D_min,D_max,NumbDevs)

        simple_save_data(dist_list, gain_list, ratios, D_n)

        return dist_list, gain_list, ratios, D_n

def simple_save_data(dist_list, gain_list, ratios, D_n):
    with h5py.File('data.h5', 'w') as hf:
        hf.create_dataset('dist_list', data=dist_list)
        hf.create_dataset('gain_list', data=gain_list)
        hf.create_dataset('ratios', data=ratios)
        hf.create_dataset('D_n', data=D_n)
        hf.close()

def simple_read_data():
    hf = h5py.File('data.h5', 'r')
    dist_list = np.array(hf.get('dist_list')[:])
    gain_list = np.array(hf.get('gain_list')[:])
    ratios = np.array(hf.get('ratios')[:])
    D_n = np.array(hf.get('D_n')[:])

    return dist_list, gain_list, ratios, D_n

def main():
    print(mobile_gen())


if __name__== "__main__":
  main()

