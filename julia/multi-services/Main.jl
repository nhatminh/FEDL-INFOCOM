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

    ### Global ###
    rs2 = Solving_global_prob(D_n,capacity)
    # rs  = [T_cmp, E_cmp, T_com, E_com, Theta]
    # println("\n---->> Check Global Solution: ", check(rs, rs2))

end

function BCD(dist_list, gain_list, capacity, D_n)
    println("**** BCD Method ****")
    T_cmp   = zeros(Numb_kaps,Numb_Services)
    E_cmp   = zeros(Numb_kaps,Numb_Services)
    f       = 1/Numb_Services *ones(Numb_kaps,Numb_Services,NumbDevs)

    T_com   = zeros(Numb_kaps,Numb_Services)
    E_com   = zeros(Numb_kaps,Numb_Services)
    w       = zeros(Numb_kaps,NumbDevs)
    tau     = zeros(Numb_kaps,Numb_Services,NumbDevs)

    Theta   = 0.1*ones(Numb_kaps,Numb_Services)
    Obj     = zeros(Numb_kaps)
    Obj_E, Obj_T    = zeros(Numb_kaps), zeros(Numb_kaps)

    for k=1:Numb_kaps
        Theta_old, w_old = zeros(Numb_Services), zeros(NumbDevs)
        global kappa = kaps[k]
        for t =1:Numb_Iteration
            print("--> Iteration: ",t )

            ### Sub1 ###
            T_cmp[k,:], f[k,:,:], E_cmp[k,:]     = Solving_sub_prob1(Theta,D_n)

            ### Sub2 ###
            T_com[k,:], w[k,:], tau[k,:,:], E_com[k,:]     = Solving_sub_prob2(Theta,capacity)

            ### Sub3 ###
            Theta[k,:], Obj[k]  = Solving_sub_prob3(T_cmp[k,:],E_cmp[k,:],T_com[k,:],E_com[k,:])

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
   end

   # # save_result(Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta)
end

function ADMM(dist_list, gain_list, capacity, D_n)
    println("**** ADMM Method ****")
    T_cmp   = zeros(Numb_kaps,Numb_Services)
    E_cmp   = zeros(Numb_kaps,Numb_Services)
    f       = 1/Numb_Services *ones(Numb_kaps,Numb_Iteration+1, Numb_Services,NumbDevs)

    T_com   = zeros(Numb_kaps,Numb_Services)
    E_com   = zeros(Numb_kaps,Numb_Services)
    w, w1   = w_min*ones(Numb_kaps,NumbDevs), w_min*ones(Numb_kaps,Numb_Services,NumbDevs)  #w1 is for each service, while w is original variable
    tau     = zeros(Numb_kaps,Numb_Services,NumbDevs)

    Theta, Theta1 = 0.1*ones(Numb_kaps,Numb_Services),  0.1*ones(Numb_kaps,Numb_Services)
    Obj, Obj1       = zeros(Numb_kaps), zeros(Numb_kaps,Numb_Services)
    Obj_E, Obj_T    = zeros(Numb_kaps,Numb_Services), zeros(Numb_kaps,Numb_Services)
    d_eta  = zeros(Numb_kaps,Numb_Services)

    r1, y  = zeros(Numb_kaps,NumbDevs), zeros(Numb_kaps,NumbDevs)    #primal residual and dual variable for SUB1
    r2, v  = zeros(Numb_kaps,Numb_Services,NumbDevs), zeros(Numb_kaps,Numb_Services,NumbDevs)    #primal residual and dual variable for SUB2
    z      = zeros(Numb_kaps,NumbDevs) # consensus variable for SUB2

    for k=1:Numb_kaps
        Theta_old, w_old = zeros(Numb_Services), zeros(NumbDevs)
        global kappa = kaps[k]

        for t =1:Numb_Iteration
            print("--> Iteration: ",t )
            f[k,t+1,:,:] = f[k,t,:,:] ## trick for cyclic update in primal ADMM

            # ### Sub1 ###
            # T_cmp[k,:], f[k,t+1,:,:], E_cmp[k,:]     = Solving_sub_prob1(Theta,D_n)

            for s=1:Numb_Services

                ### Sub1 ###
                T_cmp[k,s], f[k,t+1,s,:], E_cmp[k,s]  = Solving_isub_prob1(Theta,D_n,s,f[k,t+1,:,:],y[k,:]) # trick for cyclic update in primal ADMM
                # T_cmp[k,s], f[k,t+1,s,:], E_cmp[k,s]  = Solving_isub_prob1(Theta,D_n,s,f[k,t,:,:],y[k,:]) # parallel update in primal ADMM
                # T_cmp1[k], f1[k,:], E_cmp1[k], Tcmp_N1[k], Tcmp_N2[k], Tcmp_N3[k],N1[k], N2[k], N3[k]   = Solving_sub1(D_n[s,:])

                ### Sub2 ###
                T_com[k,s], w1[k,s,:], tau[k,s,:], E_com[k,s] = Solving_isub_prob2(Theta,capacity,s, v[k,:,:], z[k,:])
                # # T_com1[k], p1[k,:], tau1[k,:], E_com1[k] = Solving_sub2(ratios[s,:])

                ### Sub3 ###
                # Theta[k,s], Obj1[k,s]  = Solving_isub_prob3(T_cmp[k,s],E_cmp[k,s],T_com[k,s],E_com[k,s])
                Theta[k,s], Obj1[k,s], Obj_E[k,s], Obj_T[k,s], d_eta[k,s] = Solving_isub3(T_cmp[k,s],E_cmp[k,s],T_com[k,s],E_com[k,s])
            end
            # ### Sub2 ###
            # T_com[k,:], w[k,:], tau[k,:,:], E_com[k,:] = Solving_sub_prob2(Theta,capacity)

            # ### Sub3 ###
            # Theta[k,:], Obj[k]  = Solving_sub_prob3(T_cmp[k,:],E_cmp[k,:],T_com[k,:],E_com[k,:])
            # # Theta1[k], Obj1[k], Obj_E[k], Obj_T[k], d_eta[k] = Solving_sub3(T_cmp1[k],E_cmp1[k],T_com1[k],E_com1[k])

            ### Dual Update ###
            for n=1:NumbDevs
                z[k,n] = 1/Numb_Services * (sum(w1[k,:,n] + (1/RHO2)*v[k,:,n]))
                r1[k,n]   = sum(f[k,t+1,:,n]) - f_max[n]*1e-9
                y[k,n]    = y[k,n] + RHO1*r1[k,n]
                r2[k,:,n] = w1[k,:,n] - z[k,n]*ones(Numb_Services)
                v[k,:,n]    = v[k,:,n] + RHO2*r2[k,:,n]
            end
            # println("Primal Residual1:",r1)
            # println("Primal Residual2:",r2)

            if ((norm(Theta[k, :] - Theta_old) <= stop_epsilon1) && (norm(z[k, :] - w_old) <= stop_epsilon2) )
                println("Theta:", Theta[k,:])
                println("z:", z[k,:])
                println("f:", f[k,t,:,:])
                # println("Obj:", Obj[k])
                println("Obj1:", sum(Obj1[k,:]))
                break
            else
                Theta_old = Theta[k, :]
                w_old = z[k, :]
            end

        end
   end
   # # save_result(Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta)
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
