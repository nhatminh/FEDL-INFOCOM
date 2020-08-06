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
function Solving_global_prob1(D_n,capacity)
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
            @NLconstraint(prob, V_s[s]/(w[n]*up_capacity[n])+ tau_extra[s,n] + tau_dl[s] <= T_com[s])
        end
    end

    # @NLobjective(prob, Min, sum(log(1/eta[s])/(1-eta[s]) * (sum( alpha/2 * C_s[s] * D_n[s,n] *f[s,n]^2 * 1e18 for n =1:NumbDevs) + kappa*T_cmp[s]) for s =1:Numb_Services))
    # @NLobjective(prob, Min, sum(1/(1-eta[s]) * (sum( Tx_Power*V_s[s]/(w[n]*capacity[n]) for n =1:NumbDevs) + kappa*T_com[s]) for s =1:Numb_Services))
    # @NLobjective(prob, Min, sum(1/(1 - eta[s]) * (E_com[s] - log(eta[s])*E_cmp[s] + kappa * (T_com[s] - log(eta[s])*T_cmp[s])) for s=1:Numb_Services) )

    @NLobjective(prob, Min, sum(As[s]/Theta[s] * (sum( Tx_Power*V_s[s]/(w[n]*up_capacity[n]) for n =1:NumbDevs) +
                            K_l[s]*sum( alpha/2 * C_s[s] * D_n[s,n] *f[s,n]^2 * 1e18 for n =1:NumbDevs) +
                            kaps_s[s] * (T_com[s]+T_avg[s] + K_l[s]*T_cmp[s])) for s=1:Numb_Services))

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
        Obj += K_g[s] * (rs_E_com[s] +K_l[s]*rs_E_cmp[s] + kaps_s[s] * (rs_T_com[s]+T_avg[s] +K_l[s]*rs_T_cmp[s]))
    end
    println("eta: ", rs_eta)
    println("w: ",rs_w[1:5])
    println("f: ",rs_f[:,1:5])
    println("Obj-Global:",Obj)
    # println("rs_T_com:",rs_T_com)
    # println("rs_T_cmp:",rs_T_cmp)
    println("Computed-Obj:",compute_obj1(rs_f, rs_w, rs_eta, D_n, capacity))


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


function compute_obj1(rs_f, rs_w, rs_eta, D_n, capacity)
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
        rs_T_com[s] = maximum(rs_tau[s,:]+ tau_extra[s,:]) + tau_dl[s]
        rs_E_com[s] =sum(Tx_Power*rs_tau[s,:])
    end

    K_g =zeros(Numb_Services)
    Obj_E, Obj_T=zeros(Numb_Services), zeros(Numb_Services)
    Obj = 0
    for s =1:Numb_Services
        # K_g[s] = A1s[s]/(A2s[s]*rs_eta[s]^2 + B2s[s]*rs_eta[s])
        K_g[s] = 2*rho0*As[s]*(Bs[s]*rs_eta[s]^2 +1)/ (Cs[s]*rs_eta[s] - Ds[s]* rs_eta[s]^2)
        # print("rs_E_com[s]:",rs_E_com[s])
        # print("rs_E_cmp[s]:",rs_E_cmp[s])
        Obj_E[s] = rs_E_com[s] +K_l[s]*rs_E_cmp[s]
        Obj_T[s] = rs_T_com[s]+T_avg[s] +K_l[s]*rs_T_cmp[s]
        Obj += K_g[s] * (Obj_E[s]  + kaps_s[s] * Obj_T[s])
    end
    # println("rs_tau:",rs_tau)
    # println("rs_E_com:",rs_E_com)
    # println("K_l*rs_E_cmp:",K_l.*rs_E_cmp)
    # println("rs_T_com[s]:",rs_T_com)
    # println("K_l*rs_T_cmp:",K_l.*rs_T_cmp)

    # println("rs_T_cmp:",rs_T_cmp)
    return Obj_E, Obj_T, Obj
end

function compute_service_cost1(T_cmp, E_cmp, T_com, E_com)
    C_s =zeros(Numb_Services)
    for s =1:Numb_Services
        C_s[s] = E_com[s] +K_l[s]*E_cmp[s] + kaps_s[s] * (T_com[s] +T_avg[s]+K_l[s])*T_cmp[s]
    end
    return C_s
end
