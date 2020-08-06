# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 13:19:14 2017

@author: Minh
"""

import h5py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_final():
    df_iters = read_files() #3 algorithms
    # df_iters = read_files2() # 2 algorithms

    # plt.figure(2, figsize=(7., 5.1))
    plt.figure(2,figsize=(8.7,5.8))
    #plt.figure(2,figsize=(4.5,5))
    #sns.set_style("whitegrid")
    sns.set_context("notebook", font_scale=1.5)
    sns.swarmplot(x="Algorithm",y="Iterations", data=df_iters)
    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    sns.boxplot(x="Algorithm", y="Iterations", data=df_iters, showcaps=True, boxprops={'facecolor': 'None'},
                showfliers=False, whiskerprops={'linewidth': 1}, linewidth=0.8)
    plt.xlabel('Algorithm',fontsize=20)
    # plt.ylim(0, 100)
    plt.savefig("figs\Iters_Stat_100.pdf")
    plt.show()


def read_files():
    NUM_SIM = 100
    filename = 'result_iter_'+str(NUM_SIM)+'.h5'
    
    df = h5py.File(filename, 'r')

    # stop1 = df['stop1'][:]
    stop2 = df['stop2'][:]
    stop3 = df['stop3'][:]
    stop4 = df['stop4'][:]
    rs_Objs = df['rs_Objs'][:]
    # print("Avg BCD:",np.average(stop1))
    print("Avg miADMM:",np.median(stop2))
    print("Avg JP-miADMM:",np.median(stop3))
    print("Avg JP-miADMM ES:", np.median(stop4))
    # print("Obj1:",rs_Objs[:,1])
    # print("Obj2:", rs_Objs[:, 2])
    print("Obj_ratio:", np.average(rs_Objs[:,2]/rs_Objs[:,0])*100)

    data = np.concatenate((stop3, stop4, stop2 ), axis=1)
    iters_cols = ['Decentralized', 'Decentralized-ES', 'miADMM']
    # data = np.concatenate((stop1, stop3), axis=1)
    # iters_cols =['Centralized','Decentralized']
    df_iters = pd.DataFrame(data, columns=iters_cols)
    df_iters = pd.melt(df_iters,var_name='Algorithm',value_name='Iterations')
    # print(df_iters)
    return df_iters


def read_files1():
    NUM_SIM = 100
    filename = 'result_iter_' + str(NUM_SIM) + '.h5'

    df = h5py.File(filename, 'r')

    stop1 = df['stop1'][:]
    stop2 = df['stop2'][:]
    stop3 = df['stop3'][:]
    print("Avg BCD:", np.average(stop1))
    print("Avg miADMM:", np.average(stop2))
    print("Avg JP-miADMM:", np.average(stop3))

    data = np.concatenate((stop1, stop3, stop2), axis=1)
    iters_cols = ['Centralized', 'Decentralized', 'miADMM']
    # data = np.concatenate((stop1, stop3), axis=1)
    # iters_cols =['Centralized','Decentralized']
    df_iters = pd.DataFrame(data, columns=iters_cols)
    df_iters = pd.melt(df_iters, var_name='Algorithm', value_name='Iterations')
    # print(df_iters)
    return df_iters


def read_files2():
    NUM_SIM = 100
    filename = 'result_iter_' + str(NUM_SIM) + '.h5'

    df = h5py.File(filename, 'r')

    stop1 = df['stop1'][:]
    stop2 = df['stop2'][:]
    stop3 = df['stop3'][:]
    print("Avg BCD:", np.average(stop1))
    print("Avg miADMM:", np.average(stop2))
    print("Avg JP-miADMM:", np.average(stop3))
    print(stop1.shape)
    print(stop3.shape)
    # data = np.concatenate((stop1, stop3), axis=1)
    # iters_cols = ['BCD', 'JP-miADMM']
    data = stop1
    iters_cols = ['BCD']

    df_iters = pd.DataFrame(data, columns=iters_cols)
    df_iters = pd.melt(df_iters, var_name='Algorithm', value_name='Iterations')
    # print(df_iters)
    return df_iters

plot_final()