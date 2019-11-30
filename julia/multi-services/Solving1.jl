# using Convex
using ECOS
using Gurobi
using SCS
using PyPlot
using HDF5
using JuMP
using Ipopt
using LambertW
using Roots
using DataStructures
using LinearAlgebra
#include("Setting.jl")

eps = 1e-4
function Solving_global_prob(D_n,capacity)
    up_capacity = capacity[1]
    down_capacity = capacity[2]
    tau_dl=zeros(Numb_Services) # Downlink transmission time
    for s =1:Numb_Services
        tau_dl[s]= maximum(V_s[s]/down_capacity[:])
    end
    println("\n===== Solving Global: Solver =====\n")
    # println("ratios: ", ratios)

    # Model(with_optimizer(CPLEX.Optimizer, CPX_PARAM_EPGAP=1e-2))
    prob = Model(with_optimizer(Ipopt.Optimizer,tol=1e-7, max_iter=10000, print_level =0))
    # prob = Model(with_optimizer(ECOS.Optimizer))
    # prob = Model(with_optimizer(Gurobi.Optimizer))
    # prob = Model(solver=BonminNLSolver())
    # prob = Model(solver=OsilBonminSolver())

    @variable(prob, T_cmp[1:Numb_Services] >= 0 + eps)         # T computing
    @variable(prob, f[1:Numb_Services,1:NumbDevs]>= 0 + eps)   # CPU freq.

    @variable(prob, T_com[1:Numb_Services] >= 0 + eps)         # T communications
    @variable(prob, w[1:NumbDevs]>= w_min)                     # BW fraction


    @variable(prob, 0 + eps <= Theta[1:Numb_Services] <= 1 - eps)    #Local error


    for s =1:Numb_Services
        for n =1:NumbDevs
            # @constraint(prob, f[n] <= f_max[n] *1e-9)
            @constraint(prob, f[s,n] >= f_min[s] *1e-9)
            @NLconstraint(prob, C_s[s]*(D_n[s,n]*1e-9/f[s,n]) <= T_cmp[s])
        end
    end
    for n =1:NumbDevs
        @constraint(prob, sum( f[s,n] for s =1:Numb_Services) == f_max[n]*1e-9)
    end

    @constraint(prob, sum(w[n] for n=1:NumbDevs) == 1 )

    for s =1:Numb_Services
        for n=1:NumbDevs
            @NLconstraint(prob, V_s[s]/(w[n]*up_capacity[n]) + tau_dl[s] <= T_com[s])
        end
    end

    # @NLobjective(prob, Min, sum(log(1/Theta[s])/(1-Theta[s]) * (sum( alpha/2 * C_s[s] * D_n[s,n] *f[s,n]^2 * 1e18 for n =1:NumbDevs) + kappa*T_cmp[s]) for s =1:Numb_Services))
    # @NLobjective(prob, Min, sum(1/(1-Theta[s]) * (sum( Tx_Power*V_s[s]/(w[n]*capacity[n]) for n =1:NumbDevs) + kappa*T_com[s]) for s =1:Numb_Services))
    # @NLobjective(prob, Min, sum(1/(1 - Theta[s]) * (E_com[s] - log(Theta[s])*E_cmp[s] + kappa * (T_com[s] - log(Theta[s])*T_cmp[s])) for s=1:Numb_Services) )

    @NLobjective(prob, Min, sum(1/(1-Theta[s])  * (sum( Tx_Power*V_s[s]/(w[n]*up_capacity[n]) for n =1:NumbDevs) -
                            log(Theta[s])*sum( alpha/2 * C_s[s] * D_n[s,n] *f[s,n]^2 * 1e18 for n =1:NumbDevs) +
                            kappa * (T_com[s] - log(Theta[s])*T_cmp[s])) for s=1:Numb_Services))

    optimize!(prob)
    println("Solve Status: ",termination_status(prob))

    rs_T_cmp = value.(T_cmp)
    rs_f = value.(f)
    rs_E_cmp = zeros(Numb_Services)
    for s =1:Numb_Services
        rs_E_cmp[s] = alpha / 2 * sum(C_s[s]* D_n[s,:].*(rs_f[s,:].^2)) * 1e18
    end

    rs_T_com = value.(T_com)
    rs_w = value.(w)
    rs_tau =  zeros(Numb_Services,NumbDevs)
    rs_E_com = zeros(Numb_Services)

    for s =1:Numb_Services
        for n=1:NumbDevs
            rs_tau[s,n] = V_s[s]/(rs_w[n]*up_capacity[n])
        end
        rs_E_com[s] =sum(Tx_Power*rs_tau[s,:])
    end

    rs_Theta = value.(Theta)
    Obj = 0
    for s =1:Numb_Services
        Obj += 1/(1 - rs_Theta[s]) * (rs_E_com[s] - log(rs_Theta[s])*rs_E_cmp[s] + kappa * (rs_T_com[s] - log(rs_Theta[s])*rs_T_cmp[s]))
    end
    println("Theta: ", rs_Theta)
    println("w: ",rs_w)
    println("f: ",rs_f)
    println("Obj-Global:",Obj)
    # println("rs_T_com:",rs_T_com)
    # println("rs_T_cmp:",rs_T_cmp)
    println("Computed-Obj:",compute_obj(rs_f, rs_w, rs_Theta, D_n, capacity))


    if (DEBUG > 1)
        println("T_cmp: ", rs_T_cmp)
        println("f: ",rs_f)
        println("E_cmp: ", rs_E_cmp)
        println("T_com: ", rs_T_com)
        println("w: ",rs_w)
        println("tau: ", rs_tau)
        println("E_com: ", rs_E_com)
        println("Theta: ", rs_Theta)
    end

    return Obj, rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, rs_Theta, rs_w, rs_f
end


function compute_obj(rs_f, rs_w, rs_Theta, D_n, capacity)
    up_capacity = capacity[1]
    down_capacity = capacity[2]
    tau_dl=zeros(Numb_Services) # Downlink transmission time
    for s =1:Numb_Services
        tau_dl[s]= maximum(V_s[s]/down_capacity[:])
    end

    rs_E_cmp = zeros(Numb_Services)
    rs_T_cmp = zeros(Numb_Services)
    for s =1:Numb_Services
        rs_T_cmp[s] = maximum(C_s[s]*(D_n[s,:].*1e-9./rs_f[s,:]))
        rs_E_cmp[s] = alpha / 2 * sum(C_s[s]* D_n[s,:].*(rs_f[s,:].^2)) * 1e18
    end

    rs_tau =  zeros(Numb_Services,NumbDevs)
    rs_E_com = zeros(Numb_Services)
    rs_T_com = zeros(Numb_Services)

    for s =1:Numb_Services
        for n=1:NumbDevs
            rs_tau[s,n] = V_s[s]/(rs_w[n]*up_capacity[n])
        end
        rs_T_com[s] = maximum(rs_tau[s,:]) + tau_dl[s]
        rs_E_com[s] =sum(Tx_Power*rs_tau[s,:])
    end

    Obj = 0
    for s =1:Numb_Services
        Obj += 1/(1 - rs_Theta[s]) * (rs_E_com[s] - log(rs_Theta[s])*rs_E_cmp[s] + kappa * (rs_T_com[s] - log(rs_Theta[s])*rs_T_cmp[s]))
    end
    # println("rs_T_com:",rs_T_com)
    # println("rs_T_cmp:",rs_T_cmp)
    return rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, Obj
end

function Solving_sub_prob1(Theta, D_n)
    if (DEBUG > 0) println("\n===== Solving Sub1: Solver =====\n") end
    # println("D_n: ",D_n)
    prob = Model(with_optimizer(Ipopt.Optimizer, tol=1e-7, max_iter=100000, print_level =0))
    # prob = Model(with_optimizer(ECOS.Optimizer))
    # prob = Model(with_optimizer(Gurobi.Optimizer))

    @variable(prob, T_cmp[1:Numb_Services] >= 0 + eps)    # T computing
    @variable(prob, f[1:Numb_Services,1:NumbDevs]>= 0 + eps)   # CPU freq.

    for s =1:Numb_Services
        for n =1:NumbDevs
            # @constraint(prob, f[n] <= f_max[n] *1e-9)
            @constraint(prob, f[s,n] >= f_min[s] *1e-9)
            @NLconstraint(prob, C_s[s]*(D_n[s,n]*1e-9/f[s,n]) <= T_cmp[s])
        end
    end
    for n =1:NumbDevs
        @constraint(prob, sum( f[s,n] for s =1:Numb_Services) == f_max[n]*1e-9)
    end

    @NLobjective(prob, Min, sum(log(1/Theta[s])/(1-Theta[s]) * (sum( alpha/2 * C_s[s] * D_n[s,n] *f[s,n]^2 * 1e18 for n =1:NumbDevs) + kappa*T_cmp[s]) for s =1:Numb_Services))

    optimize!(prob)
    # println("Solve Status: ",termination_status(prob))

    rs_T_cmp = value.(T_cmp)
    rs_f = value.(f)
    rs_E_cmp = zeros(Numb_Services)
    for s =1:Numb_Services
        rs_E_cmp[s] = alpha / 2 * sum(C_s[s]* D_n[s,:].*(rs_f[s,:].^2)) * 1e18
    end

    if (DEBUG > 1)
        println("T_cmp: ", rs_T_cmp)
        println("f: ",rs_f)
        println("E_cmp: ", rs_E_cmp)
        println("Objective: ", rs_E_cmp + kappa*rs_T_cmp )
    end

    return rs_T_cmp, rs_f, rs_E_cmp
end

function Solving_isub_prob1(Theta, D_n, s_idx, f_k, y, jpadmm = false) ## jpadmm  for parallel update
    numb_serv = 1
    if (DEBUG > 0) println("\n===== Solving iSub1: Solver =====\n") end
    # println("D_n: ",D_n)
    prob = Model(with_optimizer(Ipopt.Optimizer, tol=1e-7, max_iter=100000, print_level =0))
    # prob = Model(with_optimizer(ECOS.Optimizer))
    # prob = Model(with_optimizer(Gurobi.Optimizer))

    @variable(prob, T_cmp[1:numb_serv] >= 0 + eps)    # T computing
    @variable(prob, f[1:numb_serv,1:NumbDevs]>= 0 + eps)   # CPU freq.
    # @variable(prob, y[1:numb_serv,1:NumbDevs]>= 0 + eps)   # dual variable

    for s =1:numb_serv
        for n =1:NumbDevs
            # @constraint(prob, f[n] <= f_max[n] *1e-9)
            @constraint(prob, f[s,n] >= f_min[s_idx] *1e-9)
            @NLconstraint(prob, C_s[s_idx]*(D_n[s_idx,n]*1e-9/f[s,n]) <= T_cmp[s])
        end
    end
    # obj_func = sum(log(1/Theta[s_idx])/(1-Theta[s_idx]) *
    #                             (sum( alpha/2 * C_s[s_idx] * D_n[s_idx,n] *f[s,n]^2 * 1e18 for n =1:NumbDevs) + kappa*T_cmp[s])
    #                             + sum(y[n]*f[s,n] for n =1:NumbDevs) for s =1:numb_serv)
    #                             + RHO1/2*sum((sum(f_k[s,n] for s =1:Numb_Services)-f_k[s_idx,n] + f[1,n] - f_max[n]*1e-9)^2 for n =1:NumbDevs)
    if(jpadmm)
        @NLobjective(prob, Min, sum(log(1/Theta[s_idx])/(1-Theta[s_idx]) *
                                    (sum( alpha/2 * C_s[s_idx] * D_n[s_idx,n] *f[s,n]^2 * 1e18 for n =1:NumbDevs) + kappa*T_cmp[s])
                                    + sum(y[n]*f[s,n] for n =1:NumbDevs) for s =1:numb_serv)
                                    + RHO1/2*sum((sum(f_k[s,n] for s =1:Numb_Services)-f_k[s_idx,n] + f[1,n] - f_max[n]*1e-9)^2 for n =1:NumbDevs)
                                    + NU/2*sum((f[1,n]-f_k[s_idx,n])^2 for n =1:NumbDevs) )

    else
        @NLobjective(prob, Min,  sum(log(1/Theta[s_idx])/(1-Theta[s_idx]) *
                                    (sum( alpha/2 * C_s[s_idx] * D_n[s_idx,n] *f[s,n]^2 * 1e18 for n =1:NumbDevs) + kappa*T_cmp[s])
                                    + sum(y[n]*f[s,n] for n =1:NumbDevs) for s =1:numb_serv)
                                    + RHO1/2*sum((sum(f_k[s,n] for s =1:Numb_Services)-f_k[s_idx,n] + f[1,n] - f_max[n]*1e-9)^2 for n =1:NumbDevs) )
    end

    # @NLobjective(prob, Min, obj_func)

    optimize!(prob)
    # println("Solve Status: ",termination_status(prob))

    rs_T_cmp = value.(T_cmp)
    rs_f = value.(f)
    rs_E_cmp = zeros(numb_serv)

    for s =1:numb_serv
        rs_E_cmp[s] = alpha / 2 * sum(C_s[s_idx]* D_n[s_idx,:].*(rs_f[s,:].^2)) * 1e18
    end

    if (DEBUG > 1)
        println("T_cmp: ", rs_T_cmp)
        println("f: ",rs_f)
        println("E_cmp: ", rs_E_cmp)
        println("Objective: ", rs_E_cmp + kappa*rs_T_cmp )
    end

    return rs_T_cmp[1], rs_f[1,:], rs_E_cmp[1]
end

function Solving_sub_prob2(Theta, capacity)
    if (DEBUG > 0) println("\n===== Solving Sub2: Solver =====\n") end
    up_capacity = capacity[1]
    down_capacity = capacity[2]
    tau_dl=zeros(Numb_Services) # Downlink transmission time
    for s =1:Numb_Services
        tau_dl[s]= maximum(V_s[s]/down_capacity[:])
    end
    prob = Model(with_optimizer(Ipopt.Optimizer,tol=1e-10, max_iter=100000, print_level =0))

    @variable(prob, T_com[1:Numb_Services] >= 0 + eps)         # T communications
    @variable(prob, w[1:NumbDevs]>= w_min)                     # BW fraction

    @constraint(prob, sum(w[n] for n=1:NumbDevs) == 1 )

    for s =1:Numb_Services
        for n=1:NumbDevs
            @NLconstraint(prob, V_s[s]/(w[n]*up_capacity[n]) + tau_dl[s] <= T_com[s])
        end
    end

    @NLobjective(prob, Min, sum(1/(1-Theta[s]) * (sum( Tx_Power*V_s[s]/(w[n]*up_capacity[n]) for n =1:NumbDevs) + kappa*T_com[s]) for s =1:Numb_Services))


    optimize!(prob)
    # println("Solve Status: ",termination_status(prob))

    rs_T_com = value.(T_com)
    rs_w = value.(w)
    rs_tau =  zeros(Numb_Services,NumbDevs)
    rs_E_com = zeros(Numb_Services)

    for s =1:Numb_Services
        for n=1:NumbDevs
            rs_tau[s,n] = V_s[s]/(rs_w[n]*up_capacity[n])
        end
        rs_E_com[s] =sum(Tx_Power*rs_tau[s,:])
    end


    if (DEBUG > 1)
        println("T_com: ", rs_T_com)
        println("w: ",rs_w)
        println("tau: ", rs_tau)
        println("E_com: ", rs_E_com)
        # println("Objective: ", rs_E_com + kappa*rs_T_com )
        # println("Objective: ", JuMP.objective_value(prob) )
    end
    # for n = 1:NumbDevs
    #     x = V_s/(rs_tau[n]*BW)
    #     println("Check 1st derivative = 0: ", ratios[n]*(e^x - 1 - x*e^x) + kappa )
    # end

    return rs_T_com, rs_w, rs_tau, rs_E_com
end

function Solving_isub_prob2(Theta, capacity, s_idx, v, z)
    numb_serv = 1
    if (DEBUG > 0) println("\n===== Solving iSub2: Solver =====\n") end
    up_capacity = capacity[1]
    down_capacity = capacity[2]
    tau_dl=zeros(Numb_Services) # Downlink transmission time
    for s =1:Numb_Services
        tau_dl[s]= maximum(V_s[s]/down_capacity[:])
    end
    prob = Model(with_optimizer(Ipopt.Optimizer,tol=1e-10, max_iter=100000, print_level =0))

    @variable(prob, T_com[1:numb_serv] >= 0 + eps)         # T communications
    @variable(prob, w[1:numb_serv,1:NumbDevs]>= w_min)                     # BW fraction

    @constraint(prob, sum(sum(w[s,n] for n=1:NumbDevs) for s=1:numb_serv) == 1 )

    for s =1:numb_serv
        for n=1:NumbDevs
            @NLconstraint(prob, V_s[s_idx]/(w[s,n]*up_capacity[n]) + tau_dl[s] <= T_com[s])
        end
    end

    @NLobjective(prob, Min, sum(1/(1-Theta[s_idx]) *
                                (sum( Tx_Power*V_s[s_idx]/(w[s,n]*up_capacity[n]) for n =1:NumbDevs) + kappa*T_com[s])
                                + sum(v[s_idx,n]*(w[s,n]-z[n]) for n =1:NumbDevs)
                                + RHO2/2*sum((w[s,n]-z[n])^2 for n =1:NumbDevs) for s =1:numb_serv) )
    optimize!(prob)
    # println("Solve Status: ",termination_status(prob))

    rs_T_com = value.(T_com)
    rs_w = value.(w)
    rs_tau =  zeros(numb_serv,NumbDevs)
    rs_E_com = zeros(numb_serv)

    for s =1:numb_serv
        for n=1:NumbDevs
            rs_tau[s,n] = V_s[s_idx]/(rs_w[s,n]*up_capacity[n])
        end
        rs_E_com[s] =sum(Tx_Power*rs_tau[s,:])
    end


    if (DEBUG > 1)
        println("T_com: ", rs_T_com)
        println("w: ",rs_w)
        println("tau: ", rs_tau)
        println("E_com: ", rs_E_com)
        # println("Objective: ", rs_E_com + kappa*rs_T_com )
        # println("Objective: ", getobjectivevalue(prob) )
    end
    # for n = 1:NumbDevs
    #     x = V_s/(rs_tau[n]*BW)
    #     println("Check 1st derivative = 0: ", ratios[n]*(e^x - 1 - x*e^x) + kappa )
    # end

    return rs_T_com[1], rs_w[1,:], rs_tau[1,:], rs_E_com[1]
end

function Solving_sub_prob3( T_cmp, E_cmp, T_com, E_com)
    if (DEBUG > 0) println("\n===== Solving Sub3: Solver =====\n") end

    prob = Model(with_optimizer(Ipopt.Optimizer,tol=1e-10, max_iter=1000000, print_level =1))

    @variable(prob, 0 + eps <= Theta[1:Numb_Services] <= 1 - eps)    #Local error

    @NLobjective(prob, Min, sum(1/(1 - Theta[s]) * (E_com[s] - log(Theta[s])*E_cmp[s] + kappa * (T_com[s] - log(Theta[s])*T_cmp[s])) for s=1:Numb_Services) )

    optimize!(prob)
    # println("Solve Status: ",termination_status(prob))

    rs_Theta = value.(Theta)
    Obj = 0
    E_Obj = zeros(Numb_Services)
    T_Obj = zeros(Numb_Services)

    for s =1:Numb_Services
        E_Obj[s] =  1/(1 - rs_Theta[s]) * (E_com[s] - log(rs_Theta[s])*E_cmp[s] )
        T_Obj[s] =  1/(1 - rs_Theta[s]) * (T_com[s] - log(rs_Theta[s])*T_cmp[s] )
        Obj += E_Obj[s] + kappa * T_Obj[s]
    end

    if (DEBUG > 0)
        println("Theta: ", rs_Theta)
        println("Obj: ", Obj)
        println("Obj1: ", JuMP.objective_value(prob))
    end

    return rs_Theta, E_Obj, T_Obj, Obj
end

function Solving_isub_prob3( T_cmp, E_cmp, T_com, E_com)
    if (DEBUG > 0) println("\n===== Solving Sub3: Solver =====\n") end
    numb_serv = 1
    prob = Model(with_optimizer(Ipopt.Optimizer,tol=1e-10, max_iter=1000000, print_level =1))

    @variable(prob, 0 + eps <= Theta[1:numb_serv] <= 1 - eps)    #Local error

    @NLobjective(prob, Min, sum(1/(1 - Theta[s]) * (E_com[s] - log(Theta[s])*E_cmp[s] + kappa * (T_com[s] - log(Theta[s])*T_cmp[s])) for s=1:numb_serv) )

    optimize!(prob)
    # println("Solve Status: ",termination_status(prob))

    rs_Theta = value.(Theta)
    Obj = 0
    for s =1:numb_serv
        Obj += 1/(1 - rs_Theta[s]) * (E_com[s] - log(rs_Theta[s])*E_cmp[s] + kappa * (T_com[s] - log(rs_Theta[s])*T_cmp[s]))
    end

    if (DEBUG > 0)
        println("Theta: ", rs_Theta)
        println("Obj: ", Obj)
        println("Obj1: ", JuMP.objective_value(prob))
    end

    return rs_Theta[1], Obj
end

############ ############ ############
############ CLOSED FORM ############
############ ############ ############

function Solving_isub3( T_cmp, E_cmp, T_com, E_com)
    if (DEBUG > 0) println("\n===== Solving iSub3: Closed Form =====\n") end
    eta   = (E_cmp + kappa*T_cmp)/(E_cmp + E_com + kappa*(T_cmp + T_com))
    # println("Eta: ", eta)
    fx(x)  = log(x) + 1/x - 1/eta

    rs_Theta = find_zero(fx,1e-7)
    # Theta = find_zeros(fx, 0+0.00001, 1-0.00001)
    # println("Roots: ", Theta)

    Obj = 1/(1 - rs_Theta) * (E_com - log(rs_Theta)*E_cmp + kappa * (T_com - log(rs_Theta)*T_cmp))
    Obj_E = 1/(1 - rs_Theta) * (E_com - log(rs_Theta)*E_cmp)
    Obj_T = 1/(1 - rs_Theta) *(T_com - log(rs_Theta)*T_cmp)

    if (DEBUG > 0) & (NumbDevs <10)
        println("fx: ", fx(rs_Theta))
        println("Theta: ", rs_Theta)
        println("Obj: ", Obj)
    end

    # println("Test sub3:", abs(log(rs_Theta) + 1/rs_Theta - 1/eta))

    return rs_Theta, Obj, Obj_E, Obj_T, 1/eta
end
