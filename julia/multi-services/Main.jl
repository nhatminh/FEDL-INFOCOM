include("Setting.jl")
# include("Setting_DS.jl")
include("TrafficGen.jl")
include("Solving1.jl")
# include("Utils.jl")
using LinearAlgebra

# Result for Section IV: Closed-Form solution
include("Plots_Figs.jl")

function BCD(dist_list, gain_list, capacity, D_n)
    println("**** BCD Method ****")
    T_cmp   = zeros(Numb_kaps,Numb_Services)
    E_cmp   = zeros(Numb_kaps,Numb_Services)
    f       = ones(Numb_kaps,Numb_Services,NumbDevs,Numb_Iteration+1)
    for n=1:NumbDevs
        f[:,:,n,:] =f_max[n]/Numb_Services*1e-9*ones(Numb_kaps, Numb_Services,1,Numb_Iteration+1)
    end

    T_com   = zeros(Numb_kaps,Numb_Services)
    E_com   = zeros(Numb_kaps,Numb_Services)
    w       = 1/NumbDevs*ones(Numb_kaps,NumbDevs,Numb_Iteration+1)
    tau     = zeros(Numb_kaps,Numb_Services,NumbDevs)

    Theta   = 0.1*ones(Numb_kaps,Numb_Services,Numb_Iteration+1)
    Obj     = zeros(Numb_kaps, Numb_Iteration+1)
    Obj_E, Obj_T    = zeros(Numb_kaps), zeros(Numb_kaps)
    stop_k  = zeros(Int32,Numb_kaps)
    Heuristic_Obj = zeros(Numb_kaps)

    for k=1:Numb_kaps
        global kappa = kaps[k]
        t_cmp,e_cmp,t_com,e_com, Obj[k,1] = compute_obj(f[k,:,:,1], w[k,:,1], Theta[k,:,1], D_n, capacity)
        _, Heuristic_Obj[k]  = Solving_sub_prob3(t_cmp,e_cmp,t_com,e_com)

        for t =1:Numb_Iteration
            print("--> Iteration: ",t )

            ### Sub1 ###
            T_cmp[k,:], f[k,:,:,t+1], E_cmp[k,:]     = Solving_sub_prob1(Theta[k,:,t],D_n)

            ### Sub2 ###
            T_com[k,:], w[k,:,t+1], tau[k,:,:], E_com[k,:]     = Solving_sub_prob2(Theta[k,:,t],capacity)

            ### Sub3 ###
            Theta[k,:,t+1], Obj[k,t+1]  = Solving_sub_prob3(T_cmp[k,:],E_cmp[k,:],T_com[k,:],E_com[k,:])

            if ((norm(Theta[k, :, t+1] - Theta[k, :, t] ) <= stop_epsilon1) && (norm(w[k, :,t+1] - w[k, :,t] ) <= stop_epsilon2))
                println("Theta:", Theta[k, :, t+1])
                println("w:", w[k,:,t+1])
                println("f:", f[k,:,:,t+1])
                println("Obj:", Obj[k,t])
                stop_k[k]= t
                break
            end
        end
   end

   # save_BCD_result()
   # # save_result(Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta)
   return Obj, Theta, w, f, stop_k, Heuristic_Obj
end

function ADMM(dist_list, gain_list, capacity, D_n, jpadmm=false)
    println("**** ADMM Method ****")
    T_cmp   = zeros(Numb_kaps,Numb_Services)
    E_cmp   = zeros(Numb_kaps,Numb_Services)
    f       = ones(Numb_kaps, Numb_Services,NumbDevs,Numb_Iteration+1)
    for n=1:NumbDevs
        f[:,:,n,:] =f_max[n]/Numb_Services*1e-9*ones(Numb_kaps, Numb_Services,1,Numb_Iteration+1)
    end

    T_com   = zeros(Numb_kaps,Numb_Services)
    E_com   = zeros(Numb_kaps,Numb_Services)
    w, w1   = 1/NumbDevs*ones(Numb_kaps,NumbDevs), 1/NumbDevs*ones(Numb_kaps,Numb_Services,NumbDevs,Numb_Iteration+1)  #w1 is for each service, while w is original variable
    tau     = zeros(Numb_kaps,Numb_Services,NumbDevs)

    Theta           = 0.1*ones(Numb_kaps,Numb_Services, Numb_Iteration+1)
    Obj, Obj1       = zeros(Numb_kaps, Numb_Iteration+1), zeros(Numb_kaps,Numb_Services)
    Obj_E, Obj_T    = zeros(Numb_kaps,Numb_Services), zeros(Numb_kaps,Numb_Services)
    d_eta  = zeros(Numb_kaps,Numb_Services)

    r1, y  = zeros(Numb_kaps,NumbDevs,Numb_Iteration+1), zeros(Numb_kaps,NumbDevs)    #primal residual and dual variable for SUB1
    r2, v  = zeros(Numb_kaps,Numb_Services,NumbDevs,Numb_Iteration+1), zeros(Numb_kaps,Numb_Services,NumbDevs)    #primal residual and dual variable for SUB2
    z      = zeros(Numb_kaps,NumbDevs,Numb_Iteration+1) # consensus variable for SUB2
    stop_k  = Numb_Iteration*ones(Int32,Numb_kaps)

    z[:,:,1] = w[:,:]

    for k=1:Numb_kaps
        w_old = zeros(NumbDevs)
        global kappa = kaps[k]
        _,_,_,_,Obj[k,1] = compute_obj(f[k,:,:,1], w[k,:], Theta[k,:,1], D_n, capacity)

        for t =1:Numb_Iteration
            print("--> Iteration: ",t )
            f[k,:,:,t+1] = f[k,:,:,t] ## trick for cyclic update in primal ADMM

            # ### Sub1 ###
            # T_cmp[k,:], f[k,t+1,:,:], E_cmp[k,:]     = Solving_sub_prob1(Theta,D_n)

            for s=1:Numb_Services
                ### Sub1 ###
                if(jpadmm)
                    T_cmp[k,s], f[k,s,:,t+1], E_cmp[k,s]  = Solving_isub_prob1(Theta[k,:,t],D_n,s,f[k,:,:,t],y[k,:], jpadmm) # parallel update in primal JPADMM
                else
                    T_cmp[k,s], f[k,s,:,t+1], E_cmp[k,s]  = Solving_isub_prob1(Theta[k,:,t],D_n,s,f[k,:,:,t+1],y[k,:],jpadmm) # trick for cyclic update in primal ADMM
                end

                ### Sub2 ###
                T_com[k,s], w1[k,s,:,t+1], tau[k,s,:], E_com[k,s] = Solving_isub_prob2(Theta[k,:,t],capacity,s, v[k,:,:], z[k,:,t])

                ### Sub3 ###
                # Theta[k,s], Obj1[k,s]  = Solving_isub_prob3(T_cmp[k,s],E_cmp[k,s],T_com[k,s],E_com[k,s])
                Theta[k,s,t+1], Obj1[k,s], Obj_E[k,s], Obj_T[k,s], d_eta[k,s] = Solving_isub3(T_cmp[k,s],E_cmp[k,s],T_com[k,s],E_com[k,s])
                Obj[k,t+1] += Obj1[k,s]
            end
            # ### Sub2 ###
            # T_com[k,:], w[k,:], tau[k,:,:], E_com[k,:] = Solving_sub_prob2(Theta,capacity)

            # ### Sub3 ###
            # Theta[k,:], Obj[k]  = Solving_sub_prob3(T_cmp[k,:],E_cmp[k,:],T_com[k,:],E_com[k,:])

            ### Dual Update ###
            for n=1:NumbDevs
                z[k,n,t+1]      = 1/Numb_Services * (sum(w1[k,:,n,t+1] + (1/RHO2)*v[k,:,n]))
                r1[k,n,t]   = sum(f[k,:,n,t+1]) - f_max[n]*1e-9
                y[k,n]      = y[k,n] + RHO1*r1[k,n,t]
                r2[k,:,n,t] = w1[k,:,n,t+1] - z[k,n,t+1]*ones(Numb_Services)
                v[k,:,n]    = v[k,:,n] + RHO2*r2[k,:,n,t]
            end
            # println("Primal Residual1:",r1)
            # println("Primal Residual2:",r2)

            if ((norm(Theta[k, :,t+1] - Theta[k, :,t]) <= stop_epsilon1) && (norm(z[k, :,t+1] - z[k, :,t]) <= stop_epsilon2) )
                println("Theta:", Theta[k, :, t+1])
                println("z:", z[k,:,t+1])
                println("f:", f[k,:,:,t+1])
                # println("Obj:", Obj[k])
                println("Obj1:", Obj[k,t])
                stop_k[k]= t
                break
            end

        end
   end
   # # save_result(Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta)
   return Obj, r1, r2, Theta, z, w1, f, stop_k
end



function main()
    jpadmm=false
    #Generate data
    dist_list, gain_list, capacity, D_n = mobile_gen()
    Obj1, Theta1, w1, f1, stop1, Heuristic_Obj = BCD(dist_list, gain_list, capacity, D_n)
    Obj2, r1, r2, Theta2, w2, ws2, f2, stop2 = ADMM(dist_list, gain_list, capacity, D_n,jpadmm)

    ### Global ###
    rs_Obj, rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, rs_Theta, rs_w, rs_f = Solving_global_prob(D_n,capacity)

    println("Heuristic_Obj:",Heuristic_Obj[1])
    plot_convergence(Obj1, Obj2, rs_Obj, r1, r2, Theta1, Theta2, rs_Theta, w1, w2, ws2, rs_w, f1, f2, rs_f, stop1, stop2)

end

function run_multiple_time()
    global REUSED_TRAFFIC = false
    stop1 =zeros(Numb_kaps,NUM_SIM)
    stop2 =zeros(Numb_kaps,NUM_SIM)
    stop3=zeros(Numb_kaps,NUM_SIM)
    for i =1:NUM_SIM
        #Generate data
        println("---- Simulation ",i)
        dist_list, gain_list, capacity, D_n = mobile_gen()
        _, _, _, _, stop1[:,i], _ = BCD(dist_list, gain_list, capacity, D_n)
        _, _, _, _, _, _, _, stop2[:,i] = ADMM(dist_list, gain_list, capacity, D_n,false)
        _, _, _, _, _, _, _, stop3[:,i] = ADMM(dist_list, gain_list, capacity, D_n,true)
    end
    save_result_iteration(stop1,stop2,stop3)

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
    run_multiple_time()
    # main()
end
