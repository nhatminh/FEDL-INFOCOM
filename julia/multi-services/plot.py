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
    df_iters = read_files()
    
    plt.figure(2,figsize=(7.,5.1))
    sns.set_style("whitegrid")
    sns.set_context("notebook", font_scale=1.7)
    sns.swarmplot(x="Algorithm",y="Iterations", data=df_iters)
    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.xlabel('')
    plt.savefig("figs\Iters_Stat_100.pdf")
    plt.show()


def read_files():
    NUM_SIM = 100
    filename = 'result_iter_'+str(NUM_SIM)+'.h5'
    
    df = h5py.File(filename, 'r')

    stop1 = df['stop1'][:]
    stop2 = df['stop2'][:]
    stop3 = df['stop3'][:]
    print("Avg BCD:",np.average(stop1))
    print("Avg miADMM:",np.average(stop2))
    print("Avg JP-miADMM:",np.average(stop3))
    data = np.concatenate((stop1, stop2, stop3), axis=1)
    iters_cols =['BCD','miADMM','JP-miADMM']
    df_iters = pd.DataFrame(data, columns=iters_cols)
    df_iters = pd.melt(df_iters,var_name='Algorithm',value_name='Iterations')
    # print(df_iters)
    return df_iters

plot_final()