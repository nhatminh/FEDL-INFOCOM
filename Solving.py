from Setting import *
import cvxpy as cvx
import numpy as np

eps = 1e-4
def Solving_global_prob(ratios, D_n):
    T_cmp = cvx.Variable()          # T computing
    f     = cvx.Variable(NumbDevs)  # CPU freq.

    T_com = cvx.Variable()          # T communications
    p     = cvx.Variable(NumbDevs)  # TX Power
    tau   = cvx.Variable(NumbDevs)  # Time slot of TDMA

    Theta = cvx.Variable()

    constraints = []
    constraints += [0 + eps <= T_com]
    constraints += [0 + eps <= Theta, Theta <= 1 - eps]
    constraints += [0 + eps <= tau, cpu_min <= f, f <= cpu_max, Ptx_Min <= p, p <= Ptx_Max]

    constraints += [cvx.sum_entries(tau) <= T_com, C_n*cvx.max_entries(D_n*cvx.inv_pos(f)) <= T_cmp]

    for n in range(NumbDevs):
        constraints += [p[n] == ratios[n]*(cvx.exp(W*cvx.inv_pos(tau[n])*np.log(2)) - 1)]

    E_com = np.transpose(tau) * p
    E_cmp = alpha/2*C_n*np.transpose(D_n)*cvx.square(f)
    T_iter = T_com + cvx.log(cvx.inv_pos(Theta))*T_cmp
    E_iter = E_com + cvx.log(cvx.inv_pos(Theta)) * E_cmp

    objective = cvx.Minimize(E_iter + kappa*T_iter)

    prob = cvx.Problem(objective, constraints)

    # result = prob.solve(solver=cvx.GUROBI,verbose=True)
    result = prob.solve(solver=cvx.ECOS,verbose=False)
#     result = prob.solve(solver=cvx.CVXOPT,verbose=True,max_iters=10000)
#     result = prob.solve(solver=cvx.CVXOPT,verbose=True,kktsolver = "robust")
#     result = prob.solve(solver=cvx.SCS,verbose=True,eps = 1e-4, max_iters = 150000)

def Solving_sub_prob1( D_n):
    print(D_n)
    T_cmp = cvx.Variable()          # T computing
    f     = cvx.Variable(NumbDevs)  # CPU freq.

    constraints =  [ T_cmp >= 0 ]
    constraints += [ cpu_min *1e-9 <= f, f <= cpu_max *1e-9]

    for n in range(NumbDevs):
        constraints += [C_n*(D_n[n]*1e-9*cvx.inv_pos(f[n])) <= T_cmp]
    # constraints += [ C_n*cvx.max_entries(D_n[:]*1e-9*cvx.inv_pos(f)) <= T_cmp]

    E_cmp = alpha/2 * C_n * np.transpose(D_n) * cvx.square(f) * 1e18

    objective = cvx.Minimize( E_cmp + kappa*T_cmp)

    prob = cvx.Problem(objective, constraints)

    # result = prob.solve(solver=cvx.GUROBI,verbose=True)
    result = prob.solve(solver=cvx.ECOS,verbose=False)
    # result = prob.solve(solver=cvx.CVXOPT,verbose=True,max_iters=10000)
#     result = prob.solve(solver=cvx.CVXOPT,verbose=True,kktsolver = "robust")
#     result = prob.solve(solver=cvx.SCS,verbose=True,eps = 1e-8, max_iters = 2000000)

    print(result)
    rs_T_cmp = T_cmp.value
    rs_f = f.value.A1
    rs_E_cmp = np.sum( alpha / 2 * C_n * D_n* np.square(rs_f)) * 1e18

    if (DEBUG > 0):
        print("Sub1 constraint: ", C_n*np.max(D_n*1e-9/rs_f))
        print("T_cmp: ", rs_T_cmp)
        print("f: ",rs_f)
        print("E_cmp: ", rs_E_cmp)

    return rs_T_cmp, rs_f, rs_E_cmp

def Solving_sub_prob2(ratios):
    T_com = cvx.Variable()          # T communications
    # p     = cvx.Variable(NumbDevs)  # TX Power
    tau   = cvx.Variable(NumbDevs)  # Time slot of TDMA

    constraints = []
    constraints += [0 + eps <= T_com, 0 + eps <= tau]
    # constraints += [Ptx_Min <= p, p <= Ptx_Max]

    constraints += [cvx.sum_entries(tau) <= T_com]

    E_com = 0

    for n in range(NumbDevs):
        # constraints += [p[n] == ratios[n]*(cvx.exp(W*cvx.inv_pos(tau[n]*BW)*np.log(2)) - 1)]
        # E_com += tau[n] * p[n]
        E_com += tau[n] * ratios[n]*(cvx.exp(W*cvx.inv_pos(tau[n])/BW*np.log(2)) - 1)
        constraints += [Ptx_Min <= ratios[n]*(cvx.exp(W*cvx.inv_pos(tau[n])/BW*np.log(2)) - 1), ratios[n]*(cvx.exp(W*cvx.inv_pos(tau[n])/BW*np.log(2)) - 1) <= Ptx_Max]
    # E_com = np.transpose(tau) * p

    objective = cvx.Minimize(E_com + kappa*T_com)

    prob = cvx.Problem(objective, constraints)

    # result = prob.solve(solver=cvx.GUROBI,verbose=True)
    result = prob.solve(solver=cvx.ECOS, verbose=False)
    # result = prob.solve(solver=cvx.CVXOPT,verbose=True,max_iters=10000)
    #     result = prob.solve(solver=cvx.CVXOPT,verbose=True,kktsolver = "robust")
    #     result = prob.solve(solver=cvx.SCS,verbose=True,eps = 1e-4, max_iters = 150000)

    rs_T_com = T_com.value
    rs_tau = tau.value.A1
    rs_p = ratios[n]*(pow(2,W/(rs_tau*BW)) - 1)
    rs_E_com = np.transpose(rs_tau) * rs_p

    if (DEBUG > 0):
        print("T_com: ", rs_T_com)
        print("p: ",rs_p)
        print("tau: ", rs_tau)
        print("E_com: ", rs_E_com)

    return rs_T_com, rs_p, rs_tau, rs_E_com

def Solving_sub_prob3( T_cmp, E_cmp, T_com, E_com):
    Theta = cvx.Variable()

    constraints = []
    constraints += [0 + eps <= Theta, Theta <= 1 - eps]

    # T_iter = T_com + cvx.log(cvx.inv_pos(Theta))*T_cmp
    # E_iter = E_com + cvx.log(cvx.inv_pos(Theta))*E_cmp
    T_iter = T_com - cvx.log(Theta)*T_cmp
    E_iter = E_com - cvx.log(Theta)*E_cmp

    # objective = cvx.Minimize(cvx.inv_pos(1-Theta)*(E_iter+ kappa*T_iter))
    objective = cvx.Minimize(cvx.inv_pos(1 - Theta) * (E_com - cvx.log(Theta)*E_cmp + kappa * (T_com - cvx.log(Theta)*T_cmp)))

    prob = cvx.Problem(objective, constraints)

    result = prob.solve(solver=cvx.GUROBI,verbose=True)
    #result = prob.solve(solver=cvx.ECOS, verbose=False)
    #result = prob.solve(solver=cvx.CVXOPT,verbose=True,max_iters=10000)
    #result = prob.solve(solver=cvx.CVXOPT,verbose=True,kktsolver = "robust")
    #result = prob.solve(solver=cvx.SCS,verbose=True,eps = 1e-4, max_iters = 150000)

    rs_Theta = Theta.value.A1
    if (DEBUG > 0):
        print("Theta: ", rs_Theta)

    return rs_Theta

def Solving_sub_prob4(T_cmp, E_cmp, T_com, E_com):
    Theta = cvx.Variable()

    constraints = []
    constraints += [0 + eps <= Theta, Theta <= 1 - eps]

    # T_iter = T_com + cvx.log(cvx.inv_pos(Theta))*T_cmp
    # E_iter = E_com + cvx.log(cvx.inv_pos(Theta))*E_cmp
    T_iter = T_com - cvx.log(Theta)*T_cmp
    E_iter = E_com - cvx.log(Theta)*E_cmp

    # objective = cvx.Minimize(cvx.inv_pos(1-Theta)*(E_iter+ kappa*T_iter))
    # objective = cvx.Minimize(cvx.inv_pos(1 - Theta) * (E_com - cvx.log(Theta)*E_cmp + kappa * (T_com - cvx.log(Theta)*T_cmp)))
    objective = cvx.Minimize(Theta*Theta)

    prob = cvx.Problem(objective, constraints)

    result = prob.solve(solver=cvx.GUROBI,verbose=True)
    #result = prob.solve(solver=cvx.ECOS, verbose=False)
    #result = prob.solve(solver=cvx.CVXOPT,verbose=True,max_iters=10000)
    #result = prob.solve(solver=cvx.CVXOPT,verbose=True,kktsolver = "robust")
    #result = prob.solve(solver=cvx.SCS,verbose=True,eps = 1e-4, max_iters = 150000)

    rs_Theta = Theta.value
    if (DEBUG > 0):
        print("Theta: ", rs_Theta)

    return rs_Theta
