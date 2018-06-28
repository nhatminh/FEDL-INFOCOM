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
include("Setting.jl")

eps = 1e-4
function Solving_global_prob(D_n, ratios)
    println("\n===== Solving Global: Solver =====\n")
    # println("ratios: ", ratios)
    prob = Model(solver=IpoptSolver(tol=1e-7, max_iter=10000, print_level =0))
    # prob = Model(solver=ECOSSolver())
    # prob = Model(solver=BonminNLSolver())
    # prob = Model(solver=OsilBonminSolver())

    @variable(prob, T_cmp >= 0 + eps)                               # T computing
    @variable(prob, cpu_min *1e-9<= f[1:NumbDevs]<= cpu_max *1e-9)  # CPU freq.

    @variable(prob, T_com >= 0 + eps)                     # T communications
    @variable(prob, tau[1:NumbDevs]>= 0 + eps)            # Time slot of TDMA

    @variable(prob, 0 + eps <= Theta <= 1 - eps)

    @constraint(prob, sum(tau[n] for n=1:NumbDevs) <= T_com )

    for n =1:NumbDevs
        @NLconstraint(prob, C_n*(D_n[n]*1e-9/f[n]) <= T_cmp)
        @NLconstraint(prob, S_n/BW /log(Ptx_Max/ratios[n] +1) <= tau[n]  <= S_n/BW/ log(Ptx_Min/ratios[n] + 1))
    end

    @NLobjective(prob, Min, 1/(1 - Theta) * (sum(tau[n] * ratios[n]*(e^(S_n/(tau[n]*BW)) - 1) for n=1:NumbDevs) -
                            log(Theta)*sum( alpha/2 * C_n * D_n[n] *f[n]^2 * 1e18 for n =1:NumbDevs) +
                            kappa * (T_com - log(Theta)*T_cmp)))

    status = solve(prob)
    println("Solve Status: ",status)

    rs_T_cmp = getvalue(T_cmp)
    rs_f = getvalue(f)[:]
    rs_E_cmp = alpha / 2 * C_n * sum(D_n.*(rs_f.^2)) * 1e18

    rs_T_com = getvalue(T_com)
    rs_tau = getvalue(tau)[:]
    rs_p = ratios.*(e.^(S_n./(rs_tau*BW)) - 1)
    rs_E_com = rs_tau'rs_p

    rs_Theta = getvalue(Theta)

    if (DEBUG > 0)
        println("Sub1 constraint: ", maximum(C_n*(D_n*1e-9./rs_f)))
        println("T_cmp: ", rs_T_cmp)
        println("f: ",rs_f)
        println("E_cmp: ", rs_E_cmp)
        println("T_com: ", rs_T_com)
        println("p: ",rs_p)
        println("tau: ", rs_tau)
        println("E_com: ", rs_E_com)
        println("Theta: ", rs_Theta)
    end

    return [rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, rs_Theta]
end

function Solving_sub_prob1(D_n)
    println("\n===== Solving Sub1: Solver =====\n")
    println("D_n: ",D_n)
    prob = Model(solver=IpoptSolver(tol=1e-7, max_iter=10000, print_level =0))

    @variable(prob, T_cmp >= 0 + eps)                               # T computing
    @variable(prob, cpu_min *1e-9<= f[1:NumbDevs]<= cpu_max *1e-9)  # CPU freq.

    for n =1:NumbDevs
        @NLconstraint(prob, C_n*(D_n[n]*1e-9/f[n]) <= T_cmp)
    end

    @NLobjective(prob, Min, sum( alpha/2 * C_n * D_n[n] *f[n]^2 * 1e18 for n =1:NumbDevs) + kappa*T_cmp)

    status = solve(prob)
    println("Solve Status: ",status)

    rs_T_cmp = getvalue(T_cmp)
    rs_f = getvalue(f)[:]
    rs_E_cmp = alpha / 2 * C_n * sum(D_n.*(rs_f.^2)) * 1e18
    # println("here: " ,sum(D_n.*(rs_f.^2)))

    if (DEBUG > 0)
        println("Sub1 constraint: ", maximum(C_n*(D_n*1e-9./rs_f)))
        println("T_cmp: ", rs_T_cmp)
        println("f: ",rs_f)
        println("E_cmp: ", rs_E_cmp)
    end

    return rs_T_cmp, rs_f, rs_E_cmp
end

function Solving_sub_prob2(ratios)
    println("\n===== Solving Sub2: Solver =====\n")
    println("ratios: ", ratios)
    prob = Model(solver=IpoptSolver(tol=1e-10, max_iter=1000000, print_level =1))

    @variable(prob, T_com >= 0 + eps)                     # T communications
    @variable(prob, tau[1:NumbDevs]>= 0 + eps)            # Time slot of TDMA
    # @variable(prob, Ptx_Min<= p[1:NumbDevs] <= Ptx_Max )  # TX Power

    @constraint(prob, sum(tau[n] for n=1:NumbDevs) <= T_com )

    # E_com  = 0
    for n =1:NumbDevs
        # @NLconstraint(prob, p[n] == ratios[n]*(2^(W/(tau[n]*BW)) - 1))
        @NLconstraint(prob, S_n/BW /log(Ptx_Max/ratios[n] +1) <= tau[n]  <= S_n/BW/ log(Ptx_Min/ratios[n] + 1))
        # println(n)
        # println("low_tau: ", S_n/BW /log(Ptx_Max/ratios[n] +1), " < up_tau:  ", S_n/BW/ log(Ptx_Min/ratios[n] + 1))
        # E_com += tau[n] * p[n]
    end
    # E_com = np.transpose(tau) * p
    # E_com = tau'p

    # objective = E_com + kappa*T_com
    @NLobjective(prob, Min, sum(tau[n] * ratios[n]*(exp(S_n/(tau[n]*BW)) - 1) for n=1:NumbDevs) + kappa*T_com)

    status = solve(prob)
    println("Solve Status: ",status)

    rs_T_com = getvalue(T_com)
    rs_tau = getvalue(tau)[:]
    rs_p = ratios.*(e.^(S_n./(rs_tau*BW)) - 1)
    # rs_p1 = zeros(NumbDevs)
    # for n =1:NumbDevs
    #     println(n)
    #     rs_p1[n] = ratios[n]*(e^(S_n/(rs_tau[n]*BW)) - 1)
    #     println("rs_tau: ", S_n/BW /log2(rs_p1[n]/ratios[n] +1) )
    #     println("rs_tau1: ", rs_tau[n])
    # end

    rs_E_com = rs_tau'rs_p

    if (DEBUG > 0)
        println("T_com: ", rs_T_com)
        println("p: ",rs_p)
        # println("p1: ",rs_p1)
        println("tau: ", rs_tau)
        println("E_com: ", rs_E_com)
        # println("Objective: ", rs_E_com + kappa*rs_T_com )
        # println("Objective: ", getobjectivevalue(prob) )
    end
    # for n = 1:NumbDevs
    #     x = S_n/(rs_tau[n]*BW)
    #     println("Check 1st derivative = 0: ", ratios[n]*(e^x - 1 - x*e^x) + kappa )
    # end

    return rs_T_com, rs_p, rs_tau, rs_E_com
end

function Solving_sub_prob3( T_cmp, E_cmp, T_com, E_com)
    println("\n===== Solving Sub3: Solver =====\n")

    prob = Model(solver=IpoptSolver(tol=1e-10, max_iter=1000000, print_level =1))

    @variable(prob, 0 + eps <= Theta <= 1 - eps)

    # T_iter = T_com + cvx.log(cvx.inv_pos(Theta))*T_cmp
    # E_iter = E_com + cvx.log(cvx.inv_pos(Theta))*E_cmp
    # T_iter = T_com - log(Theta)*T_cmp
    # E_iter = E_com - log(Theta)*E_cmp
    #
    # objective = 1/(1-Theta)*(E_iter+ kappa*T_iter)

    @NLobjective(prob, Min, 1/(1 - Theta) * (E_com - log(Theta)*E_cmp + kappa * (T_com - log(Theta)*T_cmp)))

    status = solve(prob)
    println("Solve Status: ",status)

    rs_Theta = getvalue(Theta)
    if (DEBUG > 0)
        println("Theta: ", rs_Theta)
        # println("Obj: ", 1/(1 - rs_Theta) * (E_com - log(rs_Theta)*E_cmp + kappa * (T_com - log(rs_Theta)*T_cmp)))
    end

    return rs_Theta
end

############ ############ ############
############ CLOSED FORM ############
############ ############ ############

function compute_T_cmp(D_n, list)
    tmp = 0
    for n in list
        tmp += alpha*(C_n*D_n[n])^3/kappa
    end
    return tmp ^(1/3)
end

function Solving_sub1(D_n)
    println("\n===== Solving Sub1: Closed Form =====\n")
    println("D_n: ",D_n)

    rs_T_cmp = 0
    rs_f = zeros(NumbDevs)
    f_max = cpu_max*ones(NumbDevs)
    f_min = cpu_min*ones(NumbDevs)

    if(kappa >= sum(alpha * f_max.^3))
        println("*** Sub1: CASE 1 ***")
        rs_f = f_max
        rs_T_cmp =  maximum(C_n*D_n./f_max)
    elseif (kappa <= sum(alpha * f_min.^3))
        println("*** Sub1: CASE 2 ***")
        rs_f = f_min
        rs_T_cmp =  maximum(C_n*D_n./f_min)
    else
        println("*** Sub1: CASE 3 ***")
        not_sastify  = []
        rs_T_cmp =  (sum(alpha*(C_n*D_n).^3 / kappa))^(1/3)

        # sorted_UEs = SortedDict(Dict{Float64,Float64}(), Base.Reverse)
        UEs = Dict{Int32,Float64}()
        sorted_UEs = OrderedDict{Int32,Float64}()
        for n = 1:NumbDevs
            UEs[n] = C_n*D_n[n]/cpu_min
        end

        sorted_UEs_array = sort(collect(UEs), by=x->x[2], rev=true)
        # println(sorted_UEs_array)
        for n=1:NumbDevs
            k,v = sorted_UEs_array[n]
            sorted_UEs[k] = v
        end

        println("N: ",sorted_UEs)
        i = NumbDevs
        for (k,v) in sorted_UEs
            if(sorted_UEs[k]<compute_T_cmp(D_n, keys(sorted_UEs)))
                delete!(sorted_UEs,k)
                rs_f[k] = cpu_min
            end
            i -= 1
        end

        println("N':",sorted_UEs)
        rs_T_cmp = compute_T_cmp(D_n,keys(sorted_UEs))
        for (k,v) in sorted_UEs
            rs_f[k] = C_n*D_n[k]/rs_T_cmp
        end
        # rs_T_cmp =  (sum(alpha*(C_n*D_n).^3 / kappa))^(1/3)
        # rs_f = C_n*D_n/rs_T_cmp
    end

    rs_f = rs_f * 1e-9
    rs_E_cmp = alpha / 2 * C_n * sum(D_n.*(rs_f.^2)) * 1e18
    # println("here: " ,sum(D_n.*(rs_f.^2)))

    if (DEBUG > 0)
        println("Sub1 constraint: ", maximum(C_n*(D_n*1e-9./rs_f)))
        println("T_cmp: ", rs_T_cmp)
        println("f: ",rs_f)
        println("E_cmp: ", rs_E_cmp)
    end

    return rs_T_cmp, rs_f, rs_E_cmp
end

function inv_g_func(tau, ratio)
    W_lam = (S_n/BW/tau - 1)      #S_n/BW = 0.04
    C = W_lam * exp(W_lam)                #inverse of lambert function
    # kap = C/((ratio - 1)/2 * log(2))
    kap = (C *e + 1) /ratio
    return kap
end

function g_func(kap, ratio)
    W_lam = lambertw((kap *ratio - 1)/e)
    tau = S_n/BW/(1 + W_lam)
    return tau
end

function Solving_sub2(ratios)
    println("\n===== Solving Sub2: Closed Form =====\n")
    println("ratios: ", ratios)

    rs_T_com = 0
    rs_tau = zeros(NumbDevs)
    # tau_max = zeros(NumbDevs)
    # tau_min = zeros(NumbDevs)

    for n =1:NumbDevs
        tau_max = S_n/BW /log(Ptx_Min/ratios[n] +1)
        tau_min = S_n/BW /log(Ptx_Max/ratios[n] +1)
        inv_g_max  = inv_g_func(tau_min, 1/ratios[n])
        inv_g_min  = inv_g_func(tau_max, 1/ratios[n])

        if(kappa >= inv_g_max)
            # println("*** Sub2: CASE 1 ***")
            rs_tau[n] = tau_min
        elseif (kappa <= inv_g_min)
            # println("*** Sub2: CASE 2 ***")
            rs_tau[n] = tau_max
        else
            println("Dev ", n)
            println("*** Sub2: CASE 3 ***")
            rs_tau[n] = g_func(kappa,1/ratios[n])
        end
        x = S_n/(rs_tau[n]*BW)
        # println("Check 1st derivative = 0: ", ratios[n]*(e^x - 1 - x*e^x) + kappa )
    end

    rs_T_com = sum(rs_tau)
    rs_p = ratios.*(e.^(S_n./(rs_tau*BW)) - 1)
    rs_E_com = rs_tau'rs_p

    if (DEBUG > 0)
        println("T_com: ", rs_T_com)
        println("p: ",rs_p)
        println("tau: ", rs_tau)
        println("E_com: ", rs_E_com)
        # println("Objective: ", rs_E_com + kappa*rs_T_com )
    end

    return rs_T_com, rs_p, rs_tau, rs_E_com
end

function Solving_sub3( T_cmp, E_cmp, T_com, E_com)
    println("\n===== Solving Sub3: Closed Form =====\n")
    eta   = (E_cmp + kappa*T_cmp)/(E_cmp + E_com + kappa*(T_cmp + T_com))
    println("Eta: ", eta)
    fx(x)  = log(x) + 1/x - 1/eta

    rs_Theta = find_zero(fx,0.001)
    Thetas = find_zeros(fx, 0+eps, 1-eps)
    println("Roots: ", Thetas)
    # @NLobjective(prob, Min, 1/(1 - Theta) * (E_com - log(Theta)*E_cmp + kappa * (T_com - log(Theta)*T_cmp)))
    #
    # status = solve(prob)
    # println("Solve Status: ",status)

    if (DEBUG > 0)
        println("fx: ", fx(rs_Theta))
        println("Theta: ", rs_Theta)
        # println("Obj: ", 1/(1 - rs_Theta) * (E_com - log(rs_Theta)*E_cmp + kappa * (T_com - log(rs_Theta)*T_cmp)))
    end

    return rs_Theta
end
