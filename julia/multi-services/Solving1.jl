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


    @variable(prob, eta[1:Numb_Services] >=0 + eps)             #Learning parameter
    @variable(prob, 0 + eps <= Theta[1:Numb_Services] <=1 + eps)             #Learning parameter


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
        @NLconstraint(prob, Theta[s] == ( (Cs[s]*eta[s] - Ds[s]*eta[s]^2)/(Bs[s]*eta[s]^2 +1)/(2*rho0) ))
        for n=1:NumbDevs
            @NLconstraint(prob, V_s[s]/(w[n]*up_capacity[n]) + tau_dl[s] <= T_com[s])
        end
    end

    # @NLobjective(prob, Min, sum(log(1/eta[s])/(1-eta[s]) * (sum( alpha/2 * C_s[s] * D_n[s,n] *f[s,n]^2 * 1e18 for n =1:NumbDevs) + kappa*T_cmp[s]) for s =1:Numb_Services))
    # @NLobjective(prob, Min, sum(1/(1-eta[s]) * (sum( Tx_Power*V_s[s]/(w[n]*capacity[n]) for n =1:NumbDevs) + kappa*T_com[s]) for s =1:Numb_Services))
    # @NLobjective(prob, Min, sum(1/(1 - eta[s]) * (E_com[s] - log(eta[s])*E_cmp[s] + kappa * (T_com[s] - log(eta[s])*T_cmp[s])) for s=1:Numb_Services) )

    @NLobjective(prob, Min, sum(As[s]/Theta[s] * (sum( Tx_Power*V_s[s]/(w[n]*up_capacity[n]) for n =1:NumbDevs) +
                            K_l[s]*sum( alpha/2 * C_s[s] * D_n[s,n] *f[s,n]^2 * 1e18 for n =1:NumbDevs) +
                            kappa * (T_com[s] + K_l[s]*T_cmp[s])) for s=1:Numb_Services))

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

    rs_eta = value.(eta)
    K_g =zeros(Numb_Services)
    Obj = 0
    for s =1:Numb_Services
        K_g[s] = 2*rho0*As[s]*(Bs[s]*rs_eta[s]^2 +1) /(Cs[s]*rs_eta[s] - Ds[s]*rs_eta[s]^2)
        Obj += K_g[s] * (rs_E_com[s] +K_l[s]*rs_E_cmp[s] + kappa * (rs_T_com[s] +K_l[s]*rs_T_cmp[s]))
    end
    println("eta: ", rs_eta)
    println("w: ",rs_w)
    println("f: ",rs_f)
    println("Obj-Global:",Obj)
    # println("rs_T_com:",rs_T_com)
    # println("rs_T_cmp:",rs_T_cmp)
    println("Computed-Obj:",compute_obj(rs_f, rs_w, rs_eta, D_n, capacity))


    if (DEBUG > 1)
        println("T_cmp: ", rs_T_cmp)
        println("f: ",rs_f)
        println("E_cmp: ", rs_E_cmp)
        println("T_com: ", rs_T_com)
        println("w: ",rs_w)
        println("tau: ", rs_tau)
        println("E_com: ", rs_E_com)
        println("eta: ", rs_eta)
    end

    return Obj, rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, rs_eta, rs_w, rs_f
end


function compute_obj(rs_f, rs_w, rs_eta, D_n, capacity)
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

    K_g =zeros(Numb_Services)
    Obj_E, Obj_T=zeros(Numb_Services), zeros(Numb_Services)
    Obj = 0
    for s =1:Numb_Services
        # K_g[s] = A1s[s]/(A2s[s]*rs_eta[s]^2 + B2s[s]*rs_eta[s])
        K_g[s] = 2*rho0*As[s]*(Bs[s]*rs_eta[s]^2 +1)/ (Cs[s]*rs_eta[s] - Ds[s]* rs_eta[s]^2)
        Obj_E[s] = rs_E_com[s] +K_l[s]*rs_E_cmp[s]
        Obj_T[s] = rs_T_com[s] +K_l[s]*rs_T_cmp[s]
        Obj += K_g[s] * (Obj_E[s]  + kappa * Obj_T[s])
    end
    # println("rs_tau:",rs_tau)
    # println("rs_E_com:",rs_E_com)
    # println("K_l*rs_E_cmp:",K_l.*rs_E_cmp)
    # println("rs_T_com[s]:",rs_T_com)
    # println("K_l*rs_T_cmp:",K_l.*rs_T_cmp)

    # println("rs_T_cmp:",rs_T_cmp)
    return Obj_E, Obj_T, Obj
end

function compute_service_cost(T_cmp, E_cmp, T_com, E_com)
    C_s =zeros(Numb_Services)
    for s =1:Numb_Services
        C_s[s] = E_com[s] +K_l[s]*E_cmp[s] + kappa * (T_com[s] +K_l[s])*T_cmp[s]
    end
    return C_s
end

function Solving_sub_prob1(K_g, D_n)
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

    @NLobjective(prob, Min, sum(K_l[s]*K_g[s]* (sum( alpha/2 * C_s[s] * D_n[s,n] *f[s,n]^2 * 1e18 for n =1:NumbDevs) + kappa*T_cmp[s]) for s =1:Numb_Services))

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

function Solving_isub_prob1(K_g, D_n, s_idx, f_k, y, jpadmm = false) ## jpadmm  for parallel update
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

    if(jpadmm)
        @NLobjective(prob, Min, sum(K_l[s_idx]*K_g[s_idx] *
                                    (sum( alpha/2 * C_s[s_idx] * D_n[s_idx,n] *f[s,n]^2 * 1e18 for n =1:NumbDevs) + kappa*T_cmp[s])
                                    + sum(y[n]*f[s,n] for n =1:NumbDevs) for s =1:numb_serv)
                                    + RHO1/2*sum((sum(f_k[s,n] for s =1:Numb_Services)-f_k[s_idx,n] + f[1,n] - f_max[n]*1e-9)^2 for n =1:NumbDevs)
                                    + NU/2*sum((f[1,n]-f_k[s_idx,n])^2 for n =1:NumbDevs) )

    else
        @NLobjective(prob, Min,  sum(K_l[s_idx]*K_g[s_idx] *
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

function Solving_sub_prob2(K_g, capacity)
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

    @NLobjective(prob, Min, sum(K_g[s] * (sum( Tx_Power*V_s[s]/(w[n]*up_capacity[n]) for n =1:NumbDevs) + kappa*T_com[s]) for s =1:Numb_Services))


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

function Solving_isub_prob2(K_g, capacity, s_idx, v, z)
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

    @NLobjective(prob, Min, sum(K_g[s_idx] *
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

function Solving_sub_prob3( C_s)
    if (DEBUG > 0) println("\n===== Solving Sub3: Solver =====\n") end

    prob = Model(with_optimizer(Ipopt.Optimizer,tol=1e-10, max_iter=1000000, print_level =1))

    @variable(prob, eta[1:Numb_Services] >= 0 + eps)    #Local error
    @variable(prob, 0 + eps <= Theta[1:Numb_Services] <= 1 - eps )    #Local error
    for s =1:Numb_Services
        @NLconstraint(prob, Theta[s] == ( (Cs[s]*eta[s] - Ds[s]*eta[s]^2)/(Bs[s]*eta[s]^2 +1)/(2*rho0) ))
    end

    @NLobjective(prob, Min, sum(As[s]*C_s[s]/Theta[s] for s=1:Numb_Services) )

    optimize!(prob)
    # println("Solve Status: ",termination_status(prob))

    rs_eta = value.(eta)
    Obj = 0
    E_Obj = zeros(Numb_Services)
    T_Obj = zeros(Numb_Services)

    for s =1:Numb_Services
        # E_Obj[s] =  1/(1 - rs_eta[s]) * (E_com[s] - log(rs_eta[s])*E_cmp[s] )
        # T_Obj[s] =  1/(1 - rs_eta[s]) * (T_com[s] - log(rs_eta[s])*T_cmp[s] )
        # Obj += E_Obj[s] + kappa * T_Obj[s]
        Obj += 2*rho0*As[s]*C_s[s]*(Bs[s]*rs_eta[s]^2 +1) /(Cs[s]*rs_eta[s] - Ds[s]*rs_eta[s]^2)
    end

    if (DEBUG > 0)
        println("rs_eta: ", rs_eta)
        println("Obj: ", Obj)
        # println("Obj1: ", JuMP.objective_value(prob))
    end

    return rs_eta, E_Obj, T_Obj, Obj
end

function Solving_isub_prob3(s_idx, C_s)
    if (DEBUG > 0) println("\n===== Solving iSub3: Solver =====\n") end
    numb_serv = 1
    prob = Model(with_optimizer(Ipopt.Optimizer,tol=1e-10, max_iter=1000000, print_level =1))

    @variable(prob, eta[1:numb_serv] >= 0 + eps)    #Local error
    @variable(prob, 0 + eps <= Theta[1:numb_serv] <= 1 - eps )    #Local error
    @NLconstraint(prob, Theta[1] == ( (Cs[s_idx]*eta[1] - Ds[s_idx]*eta[1]^2)/(Bs[s_idx]*eta[1]^2 + 1)/(2*rho0) ))

    # @NLobjective(prob, Min, sum(A1s[s_idx]*C_s/(A2s[s_idx]*eta[i]^2 + B2s[s_idx]*eta[i]) for i=1:numb_serv) )
    @NLobjective(prob, Min, sum(As[s_idx]*C_s/Theta[i] for i=1:numb_serv) )

    optimize!(prob)
    # println("Solve Status: ",termination_status(prob))

    rs_eta = value.(eta)
    Obj = 0
    for i =1:numb_serv
        Obj += 2*rho0*As[s_idx]*C_s *(Bs[s_idx]*rs_eta[i]^2 +1) /(Cs[s_idx]*rs_eta[i] - Ds[s_idx]*rs_eta[i]^2)
    end

    if (DEBUG > 0)
        println("rs_eta: ", rs_eta)
        println("Obj: ", Obj)
        # println("Obj1: ", JuMP.objective_value(prob))
    end

    return rs_eta[1], Obj
end

############ ############ ############
############ CLOSED FORM ############
############ ############ ############

function Solving_isub3( s_idx, C_s)
    if (DEBUG > 0) println("\n===== Solving iSub3: Closed Form =====\n") end

    rs_eta = (-Ds[s_idx] + sqrt(Ds[s_idx]^2 + Bs[s_idx]*Cs[s_idx]^2))/(Bs[s_idx]*Cs[s_idx])  # (-D+sqrt(D^2+BC^2))/BC

    Obj = 2*rho0*As[s_idx]*C_s*(Bs[s_idx]*rs_eta^2 +1) /(Cs[s_idx]*rs_eta - Ds[s_idx]*rs_eta^2)    # A*C_s*(B*eta^2 +1)/(C*eta - D*eta^2)
    Obj_E = 0
    Obj_T = 0
    # println("rs_eta: ", rs_eta)
    if (DEBUG > 0) & (NumbDevs <10)
        println("rs_eta: ", rs_eta)
        # println("fx: ", fx(rs_eta))
        println("Obj: ", Obj)
        println("Test constraint Theta:", (Cs[s_idx]*rs_eta - Ds[s_idx]*rs_eta^2)/(Bs[s_idx]*rs_eta^2 + 1)/(2*rho0))
    end



    return rs_eta, Obj, Obj_E, Obj_T
end
