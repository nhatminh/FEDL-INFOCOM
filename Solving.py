from Setting import *
import cvxpy as cvx
import numpy as np

def Solving_global_prob(ratios, D_n):
    eps = 1e-4
    T_cmp = cvx.Variable()          # T computing
    f     = cvx.Variable(NumbDevs)  # CPU freq.

    T_com = cvx.Variable()          # T communications
    p     = cvx.Variable(NumbDevs)  # TX Power
    tau   = cvx.Variable(NumbDevs)  # Time slot of TDMA

    Theta = cvx.Variable()

    constraints = []
    constraints += [0 + eps <= T_cmp, T_cmp <= 1 - eps, 0 + eps <= T_com, T_com <= 1 - eps, 0 + eps <= Theta, Theta <= 1 - eps]
    constraints += [0 + eps <= tau, cpu_min <= f, f <= cpu_max, Ptx_Min <= p, p <= Ptx_Max]

    constraints += [cvx.sum_entries(tau) <= T_com, C_n*cvx.max(D_n*cvx.inv_pos(f)) <= T_cmp]

    for n in range(NumbDevs):
        constraints += [p[n] == ratios[n]*(pow(2,W*cvx.inv_pos(tau[n])))]

    E_com = np.transpose(tau) * p
    E_cmp = cvx.sum_entries(alpha/2*C_n*np.transpose(D_n*cvx.square(f)))
    T_iter = T_com + cvx.log(cvx.inv_pos(Theta))*T_cmp

    objective = cvx.Minimize(E_com+ E_cmp + kappa*T_iter)

    prob = cvx.Problem(objective, constraints)
    #
    # #    result = prob.solve(solver=cvx.GUROBI,verbose=True)
    #     result = prob.solve(solver=cvx.ECOS,verbose=False)
    # #     result = prob.solve(solver=cvx.CVXOPT,verbose=True,max_iters=10000)
    # #     result = prob.solve(solver=cvx.CVXOPT,verbose=True,kktsolver = "robust")
    # #     result = prob.solve(solver=cvx.SCS,verbose=True,eps = 1e-4, max_iters = 150000)

#     f_scaling = f.value.A1
#     bat_val = bat.value.A1
#     if (DEBUG > 0):
#         print("bat: ", bat_val)
#         print("f: ",f_scaling)

def Solving_sub_prob1(ratios, D_n):
    return

def Solving_sub_prob2(ratios, D_n):
    return

def Solving_sub_prob3(ratios, D_n):
    return

