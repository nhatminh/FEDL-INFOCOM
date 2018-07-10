#include("Setting.jl")
using Distributions
using HDF5

function to_watt(x)
    return 10^(x / 10)
end

function to_dBW(x)
    return 10*log10(x)
end

### Exponential distribution of channel_gain
function exp_gain(dist)
    g0 = -40 #dB
    g0_W = to_watt(g0)
    d0 = 1 #1 m
    mean = g0_W*((d0/dist)^4)
    gain_h = rand(Exponential(mean))
    # println("here1: ", dist, " : ", mean, " : ",gain_h)
    return gain_h #in Watts
end

function noise_per_gain(gain)  #N0/h
    return N0/gain
end

# def shanon_capacity():
#     return

function mobile_gen()
    if(REUSED_TRAFFIC)
        return simple_read_data()
    else
        dist_list = zeros(Numb_SIMs,NumbDevs)
        gain_list = zeros(Numb_SIMs,NumbDevs)
        ratios    = zeros(Numb_SIMs,NumbDevs)
        D_n       = zeros(Numb_SIMs,NumbDevs)
        for s =1:Numb_SIMs
            dist_list[s,:] = rand(Uniform(Dist_min,Dist_max),NumbDevs)
            if(SCALE)
                D_n[s,:]   = D_Total/NumbDevs
            else
                D_n[s,:]   = rand(Uniform(D_min,D_max),NumbDevs)
            end

            for n=1:NumbDevs
                gain_list[s,n] = exp_gain(dist_list[n])
                ratios[s,n]    = noise_per_gain(gain_list[n])
            end

        end
        simple_save_data(dist_list, gain_list, ratios, D_n)
        return dist_list, gain_list, ratios, D_n
    end
end

function mobile_gen_sub1()
    if(REUSED_TRAFFIC)
        return simple_read_data("data_sub1.h5")
    else
        dist_list = zeros(Numb_D,NumbDevs)
        gain_list = zeros(Numb_D,NumbDevs)
        ratios    = zeros(Numb_D,NumbDevs)
        D_n       = zeros(Numb_D,NumbDevs)
        for d=1:Numb_D
            dist_list[d,:] = rand(Uniform(Dist_min,Dist_max),NumbDevs)

            global D_min = D_max * D_ratios[d]
            D_n[d,:]   = collect(D_min:((D_max-D_min)/(NumbDevs-1)):D_max)

            for n=1:NumbDevs
                gain_list[d,n] = exp_gain(dist_list[n])
                ratios[d,n]    = noise_per_gain(gain_list[n])
            end

        end
        D_n[:,NumbDevs] = D_max

        simple_save_data(dist_list, gain_list, ratios, D_n, "data_sub1.h5")
        return dist_list, gain_list, ratios, D_n
    end
end

function mobile_gen_sub2()
    if(REUSED_TRAFFIC)
        return simple_read_data("data_sub2.h5")
    else
        dist_list = zeros(Numb_Dis,NumbDevs)
        gain_list = zeros(Numb_Dis,NumbDevs)
        ratios    = zeros(Numb_Dis,NumbDevs)
        D_n       = zeros(Numb_Dis,NumbDevs)
        for d=1:Numb_Dis
            global Dist_max = Dist_min / D_ratios[d]
            dist_list[d,:] = collect(Dist_min:((Dist_max-Dist_min)/(NumbDevs-1)):Dist_max)

            # D_min = D_max * D_ratios[d]
            D_n[d,:]   = collect(D_min:((D_max-D_min)/(NumbDevs-1)):D_max)

            for n=1:NumbDevs
                gain_list[d,n] = exp_gain(dist_list[d,n])
                ratios[d,n]    = noise_per_gain(gain_list[d,n])
            end

        end
        # println("here: ", dist_list)

        simple_save_data(dist_list, gain_list, ratios, D_n, "data_sub2.h5")
        return dist_list, gain_list, ratios, D_n
    end
end

function simple_save_data(dist_list, gain_list, ratios, D_n, filename="data.h5")
    h5open(filename, "w") do file
        write(file,"dist_list", dist_list)
        write(file,"gain_list", gain_list)
        write(file,"ratios", ratios)
        write(file,"D_n", D_n)
    end
end

function simple_read_data(filename="data.h5")
    h5open("data.h5", "r") do file
        dist_list =read(file,"dist_list")
        gain_list =read(file,"gain_list")
        ratios =read(file,"ratios")
        D_n = read(file,"D_n")
        return dist_list, gain_list, ratios, D_n
    end
end

# mobile_gen()
