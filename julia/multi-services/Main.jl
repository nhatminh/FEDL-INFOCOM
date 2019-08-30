include("Setting.jl")
# include("Setting_DS.jl")
include("TrafficGen.jl")
include("Solving1.jl")
# include("Utils.jl")
using LinearAlgebra

# Result for Section IV: Closed-Form solution
include("Plots_Figs.jl")


function main()
    #Generate data
    dist_list, gain_list, capacity, D_n = mobile_gen()
    # BCD(dist_list, gain_list, capacity, D_n)
    ADMM(dist_list, gain_list, capacity, D_n)

end

function BCD(dist_list, gain_list, capacity, D_n)
    T_cmp, T_cmp1   = zeros(Numb_kaps,Numb_Services), zeros(Numb_kaps)
    E_cmp, E_cmp1   = zeros(Numb_kaps,Numb_Services), zeros(Numb_kaps)
    f, f1           = 1/Numb_Services *ones(Numb_kaps,Numb_Services,NumbDevs), zeros(Numb_kaps,NumbDevs)
    N1, N2, N3      = zeros(Numb_kaps), zeros(Numb_kaps), zeros(Numb_kaps)
    Tcmp_N1, Tcmp_N2, Tcmp_N3   = zeros(Numb_kaps), zeros(Numb_kaps), zeros(Numb_kaps)


    T_com, T_com1   = zeros(Numb_kaps,Numb_Services),zeros(Numb_kaps)
    E_com, E_com1   = zeros(Numb_kaps,Numb_Services), zeros(Numb_kaps)
    w, w1           = zeros(Numb_kaps,NumbDevs), zeros(Numb_kaps,NumbDevs)
    tau, tau1       = zeros(Numb_kaps,Numb_Services,NumbDevs), zeros(Numb_kaps,NumbDevs)

    Theta, Theta1 = 0.1*ones(Numb_kaps,Numb_Services),  0.1*ones(Numb_kaps,Numb_Services)
    Obj, Obj1       = zeros(Numb_kaps), zeros(Numb_kaps)
    Obj_E, Obj_T    = zeros(Numb_kaps), zeros(Numb_kaps)
    d_eta  = zeros(Numb_kaps)

    # global ZoneA = []
    # global ZoneB1 = []
    # global ZoneB2 = []
    # global ZoneC = []
    # global ZoneD = []
    # # println("Numb_kaps: ", Numb_kaps)

    for k=1:Numb_kaps
        Theta_old, w_old = zeros(Numb_Services), zeros(Numb_Services)
        global kappa = kaps[k]
        for t =1:Numb_Iteration
            print("--> Iteration: ",t )

            ### Sub1 ###
            T_cmp[k,:], f[k,:,:], E_cmp[k,:]     = Solving_sub_prob1(Theta,D_n)
            # T_cmp1[k], f1[k,:], E_cmp1[k], Tcmp_N1[k], Tcmp_N2[k], Tcmp_N3[k],N1[k], N2[k], N3[k]   = Solving_sub1(D_n[s,:])
            # println("\n---->> Check Sub1 Solution: ", check([T_cmp, f, E_cmp], [T_cmp1, f1, E_cmp1]))

            ### Sub2 ###
            T_com[k,:], w[k,:], tau[k,:,:], E_com[k,:]     = Solving_sub_prob2(Theta,capacity)
            # T_com1[k], p1[k,:], tau1[k,:], E_com1[k] = Solving_sub2(ratios[s,:])
            # println("\n---->> Check Sub2 Solution: ", check([T_com, p, tau, E_com], [T_com1, p1, tau1, E_com1]))

            ### Sub3 ###
            Theta[k,:], Obj[k]  = Solving_sub_prob3(T_cmp[k,:],E_cmp[k,:],T_com[k,:],E_com[k,:])
            # Theta1[k], Obj1[k], Obj_E[k], Obj_T[k], d_eta[k] = Solving_sub3(T_cmp1[k],E_cmp1[k],T_com1[k],E_com1[k])
            # println("\n---->> Check Sub3 Solution: ", check([Theta], [Theta1]))

            if ((norm(Theta[k, :] - Theta_old) <= stop_epsilon1) && (norm(w[k, :] - w_old) <= stop_epsilon2))
                println("Theta:", Theta[k,:])
                println("w:", w[k,:])
                println("f:", f[k,:,:])
                println("Obj:", Obj[k])
                break
            else
                Theta_old = Theta[k, :]
                w_old = w[k, :]
            end
        end

        ### Global ###
        rs2 = Solving_global_prob(D_n,capacity)
        # rs  = [T_cmp, E_cmp, T_com, E_com, Theta]
        # println("\n---->> Check Global Solution: ", check(rs, rs2))
   end

   # # save_result(Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta)
   #
   # plot_sub1_f(f1)
   # plot_sub1_T(T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3)
   # plot_sub1_N(N1, N2, N3)
   #
   # plot_sub2_p(p1)
   # plot_sub2_tau(tau1)
   #
   # plot_sub3_equation(Theta1, d_eta)
   # plot_sub3_cvx(Theta1, Obj1, T_cmp1, E_cmp1, T_com1, E_com1)
   # plot_sub3_kappa_theta(Theta1, d_eta)
   # # # plot_numerical_pareto(Theta1, T_cmp1, E_cmp1, T_com1, E_com1)
   # # println("K1: ",minimum(alpha*(f_min.^3)) )
   # # println("ZoneA: ", ZoneA)
   # # # println("ZoneB1: ", ZoneB1)
   # # println("ZoneB2: ", ZoneB2)
   # # println("ZoneC: ", ZoneC)
   # # println("ZoneD: ", ZoneD)
end

function ADMM(dist_list, gain_list, capacity, D_n)
    T_cmp, T_cmp1   = zeros(Numb_kaps,Numb_Services), zeros(Numb_kaps)
    E_cmp, E_cmp1   = zeros(Numb_kaps,Numb_Services), zeros(Numb_kaps)
    f, f1           = 1/Numb_Services *ones(Numb_kaps,Numb_Iteration+1, Numb_Services,NumbDevs), zeros(Numb_kaps,NumbDevs)
    N1, N2, N3      = zeros(Numb_kaps), zeros(Numb_kaps), zeros(Numb_kaps)
    Tcmp_N1, Tcmp_N2, Tcmp_N3   = zeros(Numb_kaps), zeros(Numb_kaps), zeros(Numb_kaps)


    T_com, T_com1   = zeros(Numb_kaps,Numb_Services),zeros(Numb_kaps)
    E_com, E_com1   = zeros(Numb_kaps,Numb_Services), zeros(Numb_kaps)
    w, w1           = w_min*ones(Numb_kaps,NumbDevs), w_min*ones(Numb_kaps,Numb_Services,NumbDevs)
    tau, tau1       = zeros(Numb_kaps,Numb_Services,NumbDevs), zeros(Numb_kaps,NumbDevs)

    Theta, Theta1 = 0.1*ones(Numb_kaps,Numb_Services),  0.1*ones(Numb_kaps,Numb_Services)
    Obj, Obj1       = zeros(Numb_kaps), zeros(Numb_kaps)
    Obj_E, Obj_T    = zeros(Numb_kaps), zeros(Numb_kaps)
    d_eta  = zeros(Numb_kaps)

    r1, y  = zeros(Numb_kaps,NumbDevs), zeros(Numb_kaps,NumbDevs)    #primal residual and dual variable for SUB1
    r2, v  = zeros(Numb_kaps,Numb_Services,NumbDevs), zeros(Numb_kaps,Numb_Services,NumbDevs)    #primal residual and dual variable for SUB2
    z      = zeros(Numb_kaps,NumbDevs) # consensus variable for SUB2

    for k=1:Numb_kaps
        Theta_old, w_old = zeros(Numb_Services), zeros(Numb_Services)
        global kappa = kaps[k]

        for t =1:Numb_Iteration
            print("--> Iteration: ",t )
            f[k,t+1,:,:] = f[k,t,:,:] ## trick for cyclic update in primal ADMM

            # ### Sub1 ###
            # T_cmp[k,:], f[k,t+1,:,:], E_cmp[k,:]     = Solving_sub_prob1(Theta,D_n)

            for s=1:Numb_Services

                ### Sub1 ###
                T_cmp[k,s], f[k,t+1,s,:], E_cmp[k,s]  = Solving_isub_prob1(Theta,D_n,s,f[k,t+1,:,:],y[k,:])
                # T_cmp1[k], f1[k,:], E_cmp1[k], Tcmp_N1[k], Tcmp_N2[k], Tcmp_N3[k],N1[k], N2[k], N3[k]   = Solving_sub1(D_n[s,:])
                # println("\n---->> Check Sub1 Solution: ", check([T_cmp, f, E_cmp], [T_cmp1, f1, E_cmp1]))

                ### Sub2 ###
                T_com[k,s], w1[k,s,:], tau[k,s,:], E_com[k,s] = Solving_isub_prob2(Theta,capacity,s, v[k,:,:], z[k,:])
                # # T_com1[k], p1[k,:], tau1[k,:], E_com1[k] = Solving_sub2(ratios[s,:])
                # # println("\n---->> Check Sub2 Solution: ", check([T_com, p, tau, E_com], [T_com1, p1, tau1, E_com1]))
            end
            # ### Sub2 ###
            # T_com[k,:], w[k,:], tau[k,:,:], E_com[k,:] = Solving_sub_prob2(Theta,capacity)

            ### Sub3 ###
            Theta[k,:], Obj[k]  = Solving_sub_prob3(T_cmp[k,:],E_cmp[k,:],T_com[k,:],E_com[k,:])
            # Theta1[k], Obj1[k], Obj_E[k], Obj_T[k], d_eta[k] = Solving_sub3(T_cmp1[k],E_cmp1[k],T_com1[k],E_com1[k])
            # println("\n---->> Check Sub3 Solution: ", check([Theta], [Theta1]))

            ### Dual Update ###
            for n=1:NumbDevs
                z[k,n] = 1/Numb_Services * (sum(w1[k,:,n] + (1/RHO2)*v[k,:,n]))
                r1[k,n]   = sum(f[k,t+1,:,n]) - f_max[n]*1e-9
                y[k,n]    = y[k,n] + RHO1*r1[k,n]
                r2[k,:,n] = w1[k,:,n] - z[k,n]*ones(Numb_Services)
                v[k,:,n]    = v[k,:,n] + RHO2*r2[k,:,n]
            end
            println("Primal Residual1:",r1)
            println("Primal Residual2:",r2)

            if ((norm(Theta[k, :] - Theta_old) <= stop_epsilon1) && (norm(w[k, :] - w_old) <= stop_epsilon2) )
                println("Theta:", Theta[k,:])
                println("z:", z[k,:])
                println("f:", f[k,t,:,:])
                println("Obj:", Obj[k])
                break
            else
                Theta_old = Theta[k, :]
                w_old = w[k, :]
            end

        end

        ### Global ###
        rs2 = Solving_global_prob(D_n,capacity)
        # rs  = [T_cmp, E_cmp, T_com, E_com, Theta]
        # println("\n---->> Check Global Solution: ", check(rs, rs2))
   end

   # # save_result(Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta)
   #
   # plot_sub1_f(f1)
   # plot_sub1_T(T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3)
   # plot_sub1_N(N1, N2, N3)
   #
   # plot_sub2_p(p1)
   # plot_sub2_tau(tau1)
   #
   # plot_sub3_equation(Theta1, d_eta)
   # plot_sub3_cvx(Theta1, Obj1, T_cmp1, E_cmp1, T_com1, E_com1)
   # plot_sub3_kappa_theta(Theta1, d_eta)
   # # # plot_numerical_pareto(Theta1, T_cmp1, E_cmp1, T_com1, E_com1)
   # # println("K1: ",minimum(alpha*(f_min.^3)) )
   # # println("ZoneA: ", ZoneA)
   # # # println("ZoneB1: ", ZoneB1)
   # # println("ZoneB2: ", ZoneB2)
   # # println("ZoneC: ", ZoneC)
   # # println("ZoneD: ", ZoneD)
end


 # Result for Section IV: Closed-Form solution in the paper (5 devs)
if READ_RESULT
    dist_list, gain_list, ratios, D_n = mobile_gen()

    Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3,
    E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1,
    d_eta = read_result(string("result",NumbDevs,".h5"))

    plot_sub1_f(f1)
    plot_sub1_T(T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3)
    plot_sub1_N(N1, N2, N3)

    plot_sub2_p(p1)
    plot_sub2_tau(tau1)

    plot_sub3_equation(Theta1, d_eta)
    plot_sub3_cvx(Theta1, Obj1, T_cmp1, E_cmp1, T_com1, E_com1)
    plot_sub3_kappa_theta(Theta1, d_eta)
    plot_numerical_pareto(Theta1, T_cmp1, E_cmp1, T_com1, E_com1)

else
    main()
end
