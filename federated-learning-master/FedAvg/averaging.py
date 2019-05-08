#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python version: 3.6

import copy
import torch
import numpy as np
from torch import nn


def average_weights(w):
    w_avg = copy.deepcopy(w[0])
    for k in w_avg.keys():
        for i in range(1, len(w)):
            w_avg[k] += w[i][k]
        w_avg[k] = torch.div(w_avg[k], len(w))
    return w_avg

def weighted_average_weights(w, num_samples):
    w_avg = copy.deepcopy(w[0]) # The first user
    total_size = np.sum(num_samples)

    for k in w_avg.keys():
        # print("w0:", w_avg[k])
        w_avg[k] *= num_samples[0] #Weighted for the first user
        # print("w0*Ds:", w_avg[k])
        # print("Numb of samples:", num_samples[0])
        # print("update weight key=", str(k))
        for i in range(1, len(w)): # The remanining users
            # print("Numb of samples:", num_samples[i])
            w_avg[k] += num_samples[i]*w[i][k]
        # print("w*Ds:", w_avg[k])
        w_avg[k] = torch.div(w_avg[k], total_size)
        # print("w_avg:", w_avg[k])
    return w_avg


def average_FSVRG_weights(w, ag_scalar, net, gpu=-1):
    """
    This method is for using FSVRG algo to update global parameters
    :param w: list of client's state_dict
    :param ag_scalar: simpilicity for A Matrix
    :param net: global net model
    :return: global state_dict
    """
    w_t = copy.deepcopy(net.state_dict())
    #print("=======================before==============================")
    #print(w_t)
    sg = {}
    total_size = np.array(np.sum([u[0] for u in w]))
    for key in w_t.keys():
        sg[key] = np.zeros(w_t[key].shape)
    for l in range(len(w)):
        for k in sg.keys():
            # += ag_scalar * w[l][0] * (w[l][1][k] - w_t[k]) / total_size
            if(gpu!= -1):
                tmp_w = (w[l][1][k] - w_t[k]).cpu()
                sg[k] = np.add(sg[k], w[l][0] * tmp_w)
            else:
                sg[k] = np.add(sg[k], w[l][0] * (w[l][1][k] - w_t[k]))#np.add(sg[k].long(), torch.div(ag_scalar * w[l][0] * (torch.add(w[l][1][k], -w_t[k])).long(), total_size.long()).long())
    for key in w_t.keys():
        if (gpu != -1):
            w_t[key] = np.add(w_t[key].cpu(), np.divide(ag_scalar * sg[key], total_size))
        else:
            w_t[key] = np.add(w_t[key], np.divide(ag_scalar * sg[key], total_size))
    #print('===========================after===================================')
    #print(w_t)
    return w_t
