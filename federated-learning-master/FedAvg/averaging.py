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


def average_FSVRG_weights(w, ag_scalar, net):
    """
    This method is for using FSVRG algo to update global parameters
    :param w: list of client's state_dict
    :param ag_scalar: simpilicity for A Matrix
    :param net: global net model
    :return: global state_dict
    """
    w_t = net.state_dict()
    sg = {}
    total_size = np.sum([u[0] for u in w])
    for key in w_t.keys():
        sg[key] = np.zeros(w_t[key].shape)
    for l in range(len(w)):
        for k in sg.keys():
            sg[k] = np.add(sg[k], np.divide(ag_scalar * w[l][0] * (np.subtract(w[l][1][k], w_t[k])), total_size))#+= ag_scalar * w[l][0] * (w[l][1][k] - w_t[k]) / total_size
    for key in w_t.keys():
        w_t[key] = np.add(w_t[key], sg[key])#+= sg[key]
    return w_t
