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
#include("Setting.jl")

eps = 1e-4
function Solving_global_prob(D_n, ratios)
    println("\n===== Solving Global: Solver =====\n")
    # println("ratios: ", ratios)
    prob = Model(with_optimizer(Ipopt.Optimizer,tol=1e-7, max_iter=10000, print_level =0))
    # prob = Model(solver=ECOSSolver())
    # prob = Model(solver=BonminNLSolver())
    # prob = Model(solver=OsilBonminSolver())

    @variable(prob, T_cmp >= 0 + eps)                   # T computing
    @variable(prob, f[1:NumbDevs] >= 0)                 # CPU freq.

    @variable(prob, T_com >= 0 + eps)                   # T communications
    @variable(prob, tau[1:NumbDevs]>= 0 + eps)          # Time slot of TDMA

    @variable(prob, 0 + eps <= Theta <= 1 - eps)

    @constraint(prob, sum(tau[n] for n=1:NumbDevs) <= T_com )

    for n =1:NumbDevs
        @constraint(prob, f[n] <= f_max[n] *1e-9)
        @constraint(prob, f[n] >= f_min[n] *1e-9)
        @NLconstraint(prob, C_n[n]*(D_n[n]*1e-9/f[n]) <= T_cmp)
        @NLconstraint(prob, S_n/BW /log(Ptx_Max/ratios[n] +1) <= tau[n]  <= S_n/BW/ log(Ptx_Min/ratios[n] + 1))
    end

    @NLobjective(prob, Min, 1/(1 - Theta) * (sum(tau[n] * ratios[n]*(e^(S_n/(tau[n]*BW)) - 1) for n=1:NumbDevs) -
                            log(Theta)*sum( alpha/2 * C_n[n] * D_n[n] *f[n]^2 * 1e18 for n =1:NumbDevs) +
                            kappa * (T_com - log(Theta)*T_cmp)))

    optimize!(prob)
    println("Solve Status: ",termination_status(prob))

    rs_T_cmp = value.(T_cmp)
    rs_f = value.(f)
    rs_E_cmp = alpha / 2 * sum(C_n .* D_n.*(rs_f.^2)) * 1e18

    rs_T_com = value.(T_com)
    rs_tau = value.(tau)
    rs_p = ratios.*(e.^(S_n./(rs_tau*BW)) - 1)
    rs_E_com = rs_tau'rs_p

    rs_Theta = value.(Theta)

    if (DEBUG > 0)
        println("Sub1 constraint: ", maximum(C_n.*(D_n*1e-9./rs_f)))
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

    prob = Model(with_optimizer(Ipopt.Optimizer,tol=1e-7, max_iter=10000, print_level =0))

    @variable(prob, T_cmp >= 0 + eps)                   # T computing
    @variable(prob, f[1:NumbDevs] >= 0)                 # CPU freq.

    for n =1:NumbDevs
        @constraint(prob, f[n] <= f_max[n] *1e-9)
        @constraint(prob, f[n] >= f_min[n] *1e-9)
        @NLconstraint(prob, C_n[n]*(D_n[n]*1e-9/f[n]) <= T_cmp)
    end

    @NLobjective(prob, Min, sum( alpha/2 * C_n[n] * D_n[n] *f[n]^2 * 1e18 for n =1:NumbDevs) + kappa*T_cmp)

    optimize!(prob)
    println("Solve Status: ",termination_status(prob))

    rs_T_cmp = value.(T_cmp)
    rs_f = value.(f)[:]
    rs_E_cmp = alpha / 2 * sum(C_n.* D_n.*(rs_f.^2)) * 1e18
    # println("here: " ,sum(D_n.*(rs_f.^2)))

    if (DEBUG > 0)
        println("Sub1 constraint: ", maximum(C_n.*(D_n*1e-9./rs_f)))
        println("T_cmp: ", rs_T_cmp)
        println("f: ",rs_f)
        println("E_cmp: ", rs_E_cmp)
        println("Objective: ", rs_E_cmp + kappa*rs_T_cmp )
    end

    return rs_T_cmp, rs_f, rs_E_cmp
end

function Solving_sub_prob2(ratios)
    println("\n===== Solving Sub2: Solver =====\n")
    println("ratios: ", ratios)
    prob = Model(with_optimizer(Ipopt.Optimizer,tol=1e-7, max_iter=10000, print_level =0))

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

    optimize!(prob)
    println("Solve Status: ",termination_status(prob))

    rs_T_com = value.(T_com)
    rs_tau = value.(tau)[:]
    rs_p = ratios.*(Base.MathConstants.e .^(S_n./(rs_tau*BW)) .- 1)
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

# function Solving_sub_prob3_old( T_cmp, E_cmp, T_com, E_com)
#     println("\n===== Solving Sub3: Solver =====\n")
#
#     prob = Model(with_optimizer(Ipopt.Optimizer,tol=1e-10, max_iter=1000000, print_level =1))
#     # prob = Model(solver=IpoptSolver(tol=1e-10, max_iter=1000000, print_level =1))
#
#     @variable(prob, 0 + eps <= Theta <= 1 - eps)
#
#     # T_iter = T_com + cvx.log(cvx.inv_pos(Theta))*T_cmp
#     # E_iter = E_com + cvx.log(cvx.inv_pos(Theta))*E_cmp
#     # T_iter = T_com - log(Theta)*T_cmp
#     # E_iter = E_com - log(Theta)*E_cmp
#     #
#     # objective = 1/(1-Theta)*(E_iter+ kappa*T_iter)
#
#     @NLobjective(prob, Min, 1/(1 - Theta) * (E_com - log(Theta)*E_cmp + kappa * (T_com - log(Theta)*T_cmp)))
#
#     optimize!(prob)
#     println("Solve Status: ",termination_status(prob))
#
#     rs_Theta = value.(Theta)
#     Obj = 1/(1 - rs_Theta) * (E_com - log(rs_Theta)*E_cmp + kappa * (T_com - log(rs_Theta)*T_cmp))
#     if (DEBUG > 0)
#         println("Theta: ", rs_Theta)
#         println("Obj: ", Obj)
#         println("Obj1: ", JuMP.objective_value(prob))
#     end
#
#     return rs_Theta, Obj
# end

function Solving_sub_prob3( T_cmp, E_cmp, T_com, E_com)
    println("\n===== Solving Sub3: Solver =====\n")

    prob = Model(with_optimizer(Ipopt.Optimizer,tol=1e-10, max_iter=1000000, print_level =1))
    # prob = Model(solver=IpoptSolver(tol=1e-10, max_iter=1000000, print_level =1))

    @variable(prob, 0 + eps <= Theta <= 1 - eps) #Upper_Theta
    @variable(prob, 0 + eps <= theta <= 1 - eps) #Small_theta
    @variable(prob, eta >= 0+ eps) #eta

    # T_iter = T_com + cvx.log(cvx.inv_pos(Theta))*T_cmp
    # E_iter = E_com + cvx.log(cvx.inv_pos(Theta))*E_cmp
    # T_iter = T_com - log(Theta)*T_cmp
    # E_iter = E_com - log(Theta)*E_cmp
    #
    # objective = 1/(1-Theta)*(E_iter+ kappa*T_iter)

    @NLconstraint(prob, Theta == 2*eta*L/beta *( ((1-theta)*beta/L)^2 - theta*(1+theta) - (1+theta)^2*eta/2 ) )

    @NLobjective(prob, Min, 1/Theta * ( E_com + (1/gamma*(log(C) - log(theta))*E_cmp + kappa * (T_com + 1/gamma*(log(C) - log(theta))*T_cmp))) )

    optimize!(prob)
    println("Solve Status: ",termination_status(prob))

    rs_Theta = value.(Theta)
    rs_theta = value.(theta)
    rs_eta = value.(eta)
    rs_Theta1 =  2*rs_eta*L/beta *( ((1-rs_theta)*beta/L)^2 - rs_theta*(1+rs_theta) - (1+rs_theta)^2*rs_eta/2 )
    # Obj = 1/ rs_Theta * (E_com - log(rs_theta)*E_cmp + kappa * (T_com - log(rs_theta)*T_cmp))
    Obj = 1/rs_Theta * ( E_com + (1/gamma*(log(C) - log(rs_theta))*E_cmp + kappa * (T_com + 1/gamma*(log(C) - log(rs_theta))*T_cmp)))

    if (DEBUG > 0)
        println("Theta: ", rs_Theta)
        println("Theta1: ", rs_Theta1)
        println("theta: ", rs_theta)
        println("eta: ", rs_eta)
        println("Obj: ", Obj)
        println("Obj1: ", JuMP.objective_value(prob))
    end

    return rs_Theta, rs_theta, rs_eta, Obj
end

function Solving_sub_prob3_search( T_cmp, E_cmp, T_com, E_com)
    x = collect(1.e-5:0.0003:0.99)
    y = collect(1.e-5:0.0003:0.99)

    obj   = zeros(size(x)[1],size(y)[1])

    min_obj = maxintfloat()
    min_theta = 0
    min_eta = 0
    min_Theta = 0
    for i=1:size(x)[1]
        for j=1:size(y)[1]
            Theta=( 2*y[j]*L/beta *( ((1-x[i])*beta/L)^2 - x[i]*(1+x[i]) - (1+x[i])^2*y[j]/2 ))
            obj[i,j] = 1/Theta * ( E_com + (1/gamma*(log(C) - log(x[i]))*E_cmp + kappa * (T_com + 1/gamma*(log(C) - log(x[i]))*T_cmp)))
            if(Theta<=0  || Theta >=1)
                obj[i,j] = 0
            else
                if min_obj >= obj[i,j]
                    min_obj = obj[i,j]
                    min_theta= x[i]
                    min_eta=y[j]
                    min_Theta = Theta
                end
            end
        end
    end

    return min_Theta, min_theta, min_eta, min_obj
end

############ ############ ############
############ CLOSED FORM ############
############ ############ ############

function T_N3k(D_n, list, k)
    tmp = 0
    N_C = []
    for n in list
        tmp += alpha*1e24*(C_n[n]*D_n[n]*1e-8)^3/kappa
        if (n!=k )
            # tmp += alpha*1e24*(C_n[n]*D_n[n]*1e-8)^3/kappa
            push!(N_C, n)
        end
    end
    return tmp^(1/3), N_C
end

function T_N3(D_n, list)
    tmp = 0
    for n in list
        tmp += alpha*1e24*(C_n[n]*D_n[n]*1e-8)^3/kappa
    end
    return tmp^(1/3)
end

function T_N1(D_n,list)
    tmp = 0
    for n in list
        tmp = max(tmp, C_n[n]*D_n[n]/f_max[n])
    end
    return tmp
end

function T_N2(D_n,list)
    tmp = 0
    for n in list
        tmp = max(tmp, C_n[n]*D_n[n]/f_min[n])
    end
    tmp1 = 0
    for n in list
        tmp1 = min(tmp, C_n[n]*D_n[n]/f_min[n])
    end
    return tmp
end

# function Solving_sub1(D_n)
#     println("\n===== Solving Sub1: Closed Form =====\n")
#     # println("D_n: ",D_n)
#
#     rs_T_cmp = 0
#     rs_f = zeros(NumbDevs)
#
#     #### ------------
#     UEs_min = Dict{Int32,Float64}()
#     UEs_max = zeros(NumbDevs)
#     sorted_UEs_min = OrderedDict{Int32,Float64}()
#
#     for n = 1:NumbDevs
#         UEs_min[n] = C_n[n]*D_n[n]/f_min[n]
#         UEs_max[n] = C_n[n]*D_n[n]/f_max[n]
#     end
#
#     N1=Int32[]
#     N3 = collect(1:NumbDevs)
#
#     sorted_UEs_min_arr = sort(collect(UEs_min), by=x->x[2])
#     # println(sorted_UEs_min_arr)
#
#     for n=1:NumbDevs
#         k,v = sorted_UEs_min_arr[n]
#         sorted_UEs_min[k] = v
#     end
#
#     i = NumbDevs
#     N2C = collect(keys(sorted_UEs_min))
#
#     N2= Int32[]
#     N2C_rs = copy(N2C)
#
#     for j in N2C
#         T_j, N2C_j = T_N3k(D_n, N3, j)
#         if (C_n[j] * D_n[j]/f_min[j] <= T_j)
#             # println("delete j:", j)
#             N3 = copy(N2C_j)
#             push!(N2,j)
#              rs_f[j] = f_min[j]
#         end
#     end
#
#     T_N3_init = T_N3(D_n, N3)
#     N1_max = maximum(UEs_max)
#     # println("HERE: ", N3)
#     if (N1_max >= T_N3_init) & (T_N3_init > 0)
#         N1 = find(a->a==N1_max, UEs_max)
#     end
#
#     for n in N1
#         rs_f[n] = f_max[n]
#     end
#
#     N3 = setdiff(N3,N1)
#
#     #### ------------
#     Tcmp_N1 = T_N1(D_n, N1)
#     Tcmp_N2 = T_N2(D_n, N2)
#     Tcmp_N3 = T_N3(D_n, N3)
#     rs_T_cmp = max(Tcmp_N1, Tcmp_N2, Tcmp_N3)
#
#     if (DEBUG > 0) & (NumbDevs <10)
#         println("N_min: ",sorted_UEs_min)
#         println("Order2: ", N2C)
#         println("N1: ", N1)
#         println("-> Tcmp1: ", Tcmp_N1)
#         println("N2: ", N2)
#         println("-> Tcmp2: ", Tcmp_N2)
#         println("N3: ", N3)
#         println("-> Tcmp3: ", Tcmp_N3)
#     end
#
#     for k in N3
#         rs_f[k] = C_n[k]*D_n[k]/rs_T_cmp
#     end
#
#     if (abs(Tcmp_N1-rs_T_cmp)<1e-6)
#         time = C_n.*D_n./f_max
#
#         for n =1:NumbDevs
#             if(abs(time[n]-rs_T_cmp)<1e-6)
#                 rs_f[n]=f_max[n]
#             else
#                 rs_f[n]= C_n[n]*D_n[n]/rs_T_cmp
#             end
#         end
#     elseif ((abs(Tcmp_N2-rs_T_cmp)<1e-6))
#         for n in N2
#             rs_f[n] = f_min[n]
#         end
#         for k in N3
#             rs_f[k] = max(C_n[k]*D_n[k]/rs_T_cmp, f_min[k])
#         end
#     end
#
#
#     rs_f = rs_f * 1e-9
#     rs_E_cmp = alpha / 2 * sum(C_n.* D_n.*(rs_f.^2)) * 1e18
#     # println("here: " ,sum(D_n.*(rs_f.^2)))
#
#     if (DEBUG > 0) & (NumbDevs <10)
#         println("Sub1 constraint: ", maximum(C_n.*(D_n*1e-9./rs_f)))
#         println("T_cmp: ", rs_T_cmp)
#         println("f: ",rs_f)
#         println("E_cmp: ", rs_E_cmp)
#         println("Objective: ", rs_E_cmp + kappa*rs_T_cmp )
#     end
#     min_N2, max_N2 = 100, 0
#     for n in N2
#         min_N2 = min(min_N2, (C_n[i]*D_n[i]/f_min[i])^3 )
#         max_N2 = max(max_N2, (C_n[i]*D_n[i]/f_min[i])^3 )
#     end
#     K3 = 0
#     for n in N3
#         K3  += alpha*1e27*((C_n[n]*D_n[n]*1e-9)^3)/maximum(C_n.*D_n./f_max)^3
#     end
#
#     if(kappa< minimum(alpha*(f_min.^3)))
#         push!(ZoneA,kappa)
#     # elseif(kappa < min_N2)
#     #     println("min K2 (Zone B): ", kappa)
#     #     push!(ZoneB1,kappa)
#     elseif(kappa < max_N2)
#         push!(ZoneB2,kappa)
#     elseif(kappa < K3)
#         push!(ZoneC,kappa)
#     else
#         push!(ZoneD,kappa)
#     end
#
#     return rs_T_cmp, rs_f, rs_E_cmp, Tcmp_N1, Tcmp_N2, Tcmp_N3, size(N1)[1], size(N2)[1], size(N3)[1]
# end

function Solving_sub1(D_n)
    println("\n===== Solving Sub1: Closed Form =====\n")
    # println("D_n: ",D_n)

    rs_T_cmp = 0
    rs_f = zeros(NumbDevs)

    #### ------------
    UEs_min = Dict{Int32,Float64}()
    UEs_max = zeros(NumbDevs)
    sorted_UEs_min = OrderedDict{Int32,Float64}()

    for n = 1:NumbDevs
        UEs_min[n] = C_n[n]*D_n[n]/f_min[n]
        UEs_max[n] = C_n[n]*D_n[n]/f_max[n]
    end

    N1=Int32[]
    N3 = collect(1:NumbDevs)

    sorted_UEs_min_arr = sort(collect(UEs_min), by=x->x[2])
    # println(sorted_UEs_min_arr)

    for n=1:NumbDevs
        k,v = sorted_UEs_min_arr[n]
        sorted_UEs_min[k] = v
    end

    i = NumbDevs
    N2C = collect(keys(sorted_UEs_min))

    N2= Int32[]
    N2C_rs = copy(N2C)

    for j in N2C
        T_N3_init = T_N3(D_n, N3)
        N1_max = maximum(UEs_max)

        if (N1_max >= T_N3_init) & (T_N3_init > 0)
            N1 = vcat(N1, findall(a->a==N1_max, UEs_max))
            N3 = setdiff(N3,N1)
        end

        T_j, N2C_j = T_N3k(D_n, N3, j)
        if (C_n[j] * D_n[j]/f_min[j] <= T_j)
            # println("delete j:", j)
            N3 = copy(N2C_j)
            push!(N2,j)
             rs_f[j] = f_min[j]
        end
    end

    for n in N1
        rs_f[n] = f_max[n]
    end


    #### ------------
    Tcmp_N1 = T_N1(D_n, N1)
    Tcmp_N2 = T_N2(D_n, N2)
    Tcmp_N3 = T_N3(D_n, N3)
    rs_T_cmp = max(Tcmp_N1, Tcmp_N2, Tcmp_N3)

    if (DEBUG > 0) & (NumbDevs <10)
        println("N_min: ",sorted_UEs_min)
        println("Order2: ", N2C)
        println("N1: ", N1)
        println("-> Tcmp1: ", Tcmp_N1)
        println("N2: ", N2)
        println("-> Tcmp2: ", Tcmp_N2)
        println("N3: ", N3)
        println("-> Tcmp3: ", Tcmp_N3)
    end

    for k in N3
        rs_f[k] = C_n[k]*D_n[k]/rs_T_cmp
    end

    if (abs(Tcmp_N1-rs_T_cmp)<1e-6)
        time = C_n.*D_n./f_max

        for n =1:NumbDevs
            if(abs(time[n]-rs_T_cmp)<1e-6)
                rs_f[n]=f_max[n]
            else
                rs_f[n]= C_n[n]*D_n[n]/rs_T_cmp
            end
        end
    elseif ((abs(Tcmp_N2-rs_T_cmp)<1e-6))
        for n in N2
            rs_f[n] = f_min[n]
        end
        for k in N3
            rs_f[k] = max(C_n[k]*D_n[k]/rs_T_cmp, f_min[k])
        end
    end


    rs_f = rs_f * 1e-9
    rs_E_cmp = alpha / 2 * sum(C_n.* D_n.*(rs_f.^2)) * 1e18
    # println("here: " ,sum(D_n.*(rs_f.^2)))

    if (DEBUG > 0) & (NumbDevs <10)
        println("Sub1 constraint: ", maximum(C_n.*(D_n*1e-9./rs_f)))
        println("T_cmp: ", rs_T_cmp)
        println("f: ",rs_f)
        println("E_cmp: ", rs_E_cmp)
        println("Objective: ", rs_E_cmp + kappa*rs_T_cmp )
    end
    min_N2, max_N2 = 100, 0
    for n in N2
        min_N2 = min(min_N2, (C_n[i]*D_n[i]/f_min[i])^3 )
        max_N2 = max(max_N2, (C_n[i]*D_n[i]/f_min[i])^3 )
    end
    K3 = 0
    K2 = 0
    for n in N3
        K2 += alpha*1e27*((C_n[n]*D_n[n]*1e-9)^3)/min_N2
        K3 += alpha*1e27*((C_n[n]*D_n[n]*1e-9)^3)/maximum(C_n.*D_n./f_max)^3
    end

    if(kappa< minimum(alpha*(f_min.^3)))
        push!(ZoneA,kappa)
    # elseif(kappa < min_N2)
    #     println("min K2 (Zone B): ", kappa)
    #     push!(ZoneB1,kappa)
    # elseif(kappa < max_N2)
    #     push!(ZoneB2,kappa)
    elseif(kappa < max_N2)
        push!(ZoneB2,kappa)
    elseif(kappa < K3)
        push!(ZoneC,kappa)
    else
        push!(ZoneD,kappa)
    end

    return rs_T_cmp, rs_f, rs_E_cmp, Tcmp_N1, Tcmp_N2, Tcmp_N3, size(N1)[1], size(N2)[1], size(N3)[1]
end

function inv_g_func(tau, ratio)
    W_lam = (S_n/BW/tau - 1)      #S_n/BW = 0.04
    C = W_lam * exp(W_lam)                #inverse of lambert function
    # kap = C/((ratio - 1)/2 * log(2))
    kap = (C *Base.MathConstants.e  + 1) /ratio
    return kap
end

function g_func(kap, ratio)
    W_lam = lambertw((kap *ratio - 1)/Base.MathConstants.e )
    tau = S_n/BW/(1 + W_lam)
    return tau
end

function Solving_sub2(ratios)
    println("\n===== Solving Sub2: Closed Form =====\n")
    # println("ratios: ", ratios)

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
            # println("Dev ", n)
            # println("*** Sub2: CASE 3 ***")
            rs_tau[n] = g_func(kappa,1/ratios[n])
        end
        x = S_n/(rs_tau[n]*BW)
        # println("Check 1st derivative = 0: ", ratios[n]*(e^x - 1 - x*e^x) + kappa )
    end

    rs_T_com = sum(rs_tau)
    rs_p = ratios.*(Base.MathConstants.e .^(S_n./(rs_tau*BW)) .- 1)
    rs_E_com = rs_tau'rs_p

    if (DEBUG > 0) & (NumbDevs <10)
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

    rs_Theta = find_zero(fx,1e-6)
    # Thetas = find_zeros(fx, 0+0.00001, 1-0.00001)
    # println("Roots: ", Thetas)

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
