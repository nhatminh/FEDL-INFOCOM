include("Setting.jl")
# include("Setting_DS.jl")
include("TrafficGen.jl")
include("Solving.jl")
include("Solving_hete_kap.jl")
# include("Utils.jl")
using LinearAlgebra

# Result for Section IV: Closed-Form solution
include("Plots_Figs.jl")

function pareto(D_n,capacity)
    rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, rs_eta = zeros(Numb_kaps1,Numb_Services),zeros(Numb_kaps1,Numb_Services),
    zeros(Numb_kaps1,Numb_Services),zeros(Numb_kaps1,Numb_Services),zeros(Numb_kaps1,Numb_Services)
    for k=1:Numb_kaps1
        global kappa = kaps_pareto[k]
        _, rs_T_cmp[k,:], rs_E_cmp[k,:], rs_T_com[k,:], rs_E_com[k,:], rs_eta[k,:], _, _ = Solving_global_prob(D_n,capacity)
        # rs_Obj, rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, rs_eta, rs_w, rs_f = Solving_global_prob(D_n,capacity)
    end
    return rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, rs_eta
end

function priority(D_n,capacity)
    rs_T_cmp, rs_T_com, rs_eta = zeros(Numb_Kap_Vals,Numb_Services), zeros(Numb_Kap_Vals,Numb_Services), zeros(Numb_Kap_Vals,Numb_Services)
    rs_E_cmp, rs_E_com =  zeros(Numb_Kap_Vals,Numb_Services), zeros(Numb_Kap_Vals,Numb_Services)
    for k =1:Numb_Kap_Vals
        kaps_s[3] = k
        _, rs_T_cmp[k,:], rs_E_cmp[k,:], rs_T_com[k,:], rs_E_com[k,:], rs_eta[k,:], _, _ = Solving_global_prob1(D_n,capacity)
    end
    return rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, rs_eta
end

function BCD(dist_list, gain_list, capacity, D_n)
    println("**** BCD Method ****")
    T_cmp   = ones(Numb_kaps,Numb_Services)
    E_cmp   = ones(Numb_kaps,Numb_Services)
    f       = ones(Numb_kaps,Numb_Services,NumbDevs,Numb_Iteration+1) # Local CPU equal allocation
    f1      = zeros(Numb_kaps,Numb_Services,NumbDevs,Numb_Iteration+1) # Local CPU Proportional allocation
    w       = 1/NumbDevs*ones(Numb_kaps,NumbDevs,Numb_Iteration+1)
    w1      = zeros(Numb_kaps,NumbDevs,Numb_Iteration+1)

    up_capacity = capacity[1]
    w_total_n = sum(up_capacity)
    sorted_up_capacity =  sort(up_capacity) #small_to_large

    for n=1:NumbDevs
        f[:,:,n,:]  =f_max[n]/Numb_Services*1e-9*ones(Numb_kaps, Numb_Services,1,Numb_Iteration+1)
        # Proportional allocation for local CPU and BW
        D_total_n = sum(D_n[:,n])
        # D_total_n = 0
        # for s=1:Numb_Services
        #     D_total_n += C_s[s]*D_n[s,n]
        # end
        idx = searchsortedfirst(sorted_up_capacity,up_capacity[n]) # return indexed in the sorted list

        w1[:,n,:] = sorted_up_capacity[NumbDevs-idx+1]/w_total_n*ones(Numb_kaps,1,Numb_Iteration+1)
        for s=1:Numb_Services
            # print(" fraction:",D_n[s,n]/D_total_n)
            f1[:,s,n,:] =D_n[s,n]/D_total_n * f_max[n]*1e-9*ones(Numb_kaps, 1 ,1,Numb_Iteration+1)
            # f1[:,s,n,:] =C_s[s]*D_n[s,n]/D_total_n * f_max[n]*1e-9*ones(Numb_kaps, 1 ,1,Numb_Iteration+1)
        end
    end
    # print("w1:",w1[1,:,1])
    T_com   = ones(Numb_kaps,Numb_Services)
    E_com   = ones(Numb_kaps,Numb_Services)
    tau     = zeros(Numb_kaps,Numb_Services,NumbDevs)

    eta   = 0.1*ones(Numb_kaps,Numb_Services,Numb_Iteration+1)
    Obj     = zeros(Numb_kaps, Numb_Iteration+1)
    Obj_E, Obj_T    = zeros(Numb_kaps, Numb_Services), zeros(Numb_kaps, Numb_Services)
    stop_k  = zeros(Int32,Numb_kaps)

    ### Local CPU equal allocation
    Heuristic_Obj = zeros(Numb_kaps)
    Heuristic_Obj_E, Heuristic_Obj_T, Heuristic_eta = zeros(Numb_kaps, Numb_Services), zeros(Numb_kaps, Numb_Services), zeros(Numb_Services)
    ### Local CPU proportional allocation
    Heuristic_Obj1 = zeros(Numb_kaps)
    Heuristic_Obj_E1, Heuristic_Obj_T1,Heuristic_eta1 = zeros(Numb_kaps, Numb_Services), zeros(Numb_kaps, Numb_Services), zeros(Numb_Services)
    for k=1:Numb_kaps
        global kappa = kaps[k]

        ### BW, Local CPU proportional allocation
        obj_e,obj_t, Obj[k,1] = compute_obj(f1[k,:,:,1], w1[k,:,1], eta[k,:,1], D_n, capacity)
        C_s = obj_e  .+ kappa .* obj_t
        Heuristic_eta1[:], _, _, _  = Solving_sub_prob3(C_s)
        Heuristic_Obj_E1[k,:], Heuristic_Obj_T1[k,:], Heuristic_Obj1[k] = compute_obj(f1[k,:,:,1], w1[k,:,1], Heuristic_eta1, D_n, capacity)

        ### Local CPU equal allocation
        obj_e,obj_t, Obj[k,1] = compute_obj(f[k,:,:,1], w[k,:,1], eta[k,:,1], D_n, capacity)
        C_s = obj_e  .+ kappa .* obj_t
        Heuristic_eta[:], _, _, _  = Solving_sub_prob3(C_s)
        Heuristic_Obj_E[k,:], Heuristic_Obj_T[k,:], Heuristic_Obj[k] = compute_obj(f[k,:,:,1], w[k,:,1], Heuristic_eta, D_n, capacity)

        K_g=zeros(Numb_Services)

        for t =1:Numb_Iteration
            print("--> Iteration: ",t )
            ### Sub3 ###
            C_s = compute_service_cost(T_cmp[k,:],E_cmp[k,:],T_com[k,:],E_com[k,:])
            eta[k,:,t+1], Obj_E[k,:], Obj_T[k,:], _  = Solving_sub_prob3(C_s)
            # println("eta:",eta[k,:,t+1])
            for s=1:Numb_Services
                # K_g[s] =  A1s[s]/(A2s[s]*eta[k,s,t+1]^2 + B2s[s]*eta[k,s,t+1])
                tmp_eta = eta[k,s,t+1]
                K_g[s] =  2*rho0*As[s]*(Bs[s]*tmp_eta^2 +1)/ (Cs[s]*tmp_eta - Ds[s]* tmp_eta^2) # Kg = 2*rho*A*(B*eta^2 +1)/(C*eta - D*eta^2)
            end

            ### Sub1 ###
            T_cmp[k,:], f[k,:,:,t+1], E_cmp[k,:]     = Solving_sub_prob1(K_g,D_n)

            ### Sub2 ###
            T_com[k,:], w[k,:,t+1], tau[k,:,:], E_com[k,:]     = Solving_sub_prob2(K_g,capacity)

            Obj_E[k,:],Obj_T[k,:], Obj[k,t+1] = compute_obj(f[k,:,:,t+1], w[k,:,t+1], eta[k,:,t+1], D_n, capacity)

            if ((norm(f[k, :, :, t+1] - f[k,:, :, t] ) <= stop_epsilon1) && (norm(w[k, :,t+1] - w[k, :,t] ) <= stop_epsilon2))
                println("\neta:", eta[k, :, t+1])
                println("w:", w[k,1:5,t+1])
                println("f:", f[k,:,1:5,t+1])
                println("Obj:", Obj[k,t])
                stop_k[k]= t
                break
            end
        end
   end

   # save_BCD_result()
   # # save_result(eta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta)
   return Obj, Obj_E, Obj_T, eta, w, f, stop_k, Heuristic_Obj, Heuristic_Obj_E, Heuristic_Obj_T, Heuristic_eta,
   Heuristic_Obj1, Heuristic_Obj_E1, Heuristic_Obj_T1, Heuristic_eta1
end

function ADMM(dist_list, gain_list, capacity, D_n, jpadmm=false, early_stop=false)
    if(jpadmm == true)
        println(string("**** JP-miADMM Method **** early_stop:",early_stop))
    else
        println("**** miADMM Method ****")
    end
    T_cmp   = ones(Numb_kaps,Numb_Services)
    E_cmp   = ones(Numb_kaps,Numb_Services)
    f       = ones(Numb_kaps, Numb_Services,NumbDevs,Numb_Iteration+1)
    for n=1:NumbDevs
        f[:,:,n,:] =f_max[n]/Numb_Services*1e-9*ones(Numb_kaps, Numb_Services,1,Numb_Iteration+1)
    end

    T_com   = ones(Numb_kaps,Numb_Services)
    E_com   = ones(Numb_kaps,Numb_Services)
    w, w1   = 1/NumbDevs*ones(Numb_kaps,NumbDevs), 1/NumbDevs*ones(Numb_kaps,Numb_Services,NumbDevs,Numb_Iteration+1)  #w1 is for each service, while w is original variable
    tau     = zeros(Numb_kaps,Numb_Services,NumbDevs)

    eta           = 0.1*ones(Numb_kaps,Numb_Services, Numb_Iteration+1)
    Obj, Obj1       = zeros(Numb_kaps, Numb_Iteration+1), zeros(Numb_kaps,Numb_Services)
    Obj_E, Obj_T    = zeros(Numb_kaps,Numb_Services), zeros(Numb_kaps,Numb_Services)
    d_eta  = zeros(Numb_kaps,Numb_Services)

    r1, y  = zeros(Numb_kaps,NumbDevs,Numb_Iteration+1), zeros(Numb_kaps,NumbDevs)    #primal residual and dual variable for SUB1
    r2, v  = zeros(Numb_kaps,Numb_Services,NumbDevs,Numb_Iteration+1), zeros(Numb_kaps,Numb_Services,NumbDevs)    #primal residual and dual variable for SUB2
    z      = zeros(Numb_kaps,NumbDevs,Numb_Iteration+1) # consensus variable for SUB2
    stop_k  = Numb_Iteration*ones(Int32,Numb_kaps)

    z[:,:,1] = w[:,:]
    K_g =zeros(Numb_Services)
    for k=1:Numb_kaps
        w_old = zeros(NumbDevs)
        global kappa = kaps[k]
        _,_,Obj[k,1] = compute_obj(f[k,:,:,1], w[k,:], eta[k,:,1], D_n, capacity)

        for t =1:Numb_Iteration
            print("--> Iteration: ",t )
            f[k,:,:,t+1] = f[k,:,:,t] ## trick for cyclic update in primal ADMM
            eta1=zeros(Numb_Services)
            if(t>2) #eta does not change
                eta[k,:,t+1] = eta[k,:,t]
            else
                for s=1:Numb_Services
                    ### Sub3 ###
                    C_s = E_com[k,s] +K_l[s]*E_cmp[k,s] + kappa * (T_com[k,s]+T_avg[s] +K_l[s])*T_cmp[k,s]
                    # eta[k,s], Obj1[k,s]  = Solving_isub_prob3(s, C_s)

                    eta[k,s,t+1], Obj1[k,s], Obj_E[k,s], Obj_T[k,s] = Solving_isub3(s, C_s)
                    eta1[s],_ = Solving_isub_prob3(s, C_s)
                    # K_g[s] =  A1s[s]/(A2s[s]*eta[k,s,t+1]^2 + B2s[s]*eta[k,s,t+1])
                    tmp_eta = eta[k,s,t+1]
                    K_g[s] =  2*rho0*As[s]*(Bs[s]*tmp_eta^2 +1)/ (Cs[s]*tmp_eta - Ds[s]* tmp_eta^2)  # Kg = 2*rho0*A*(B*eta^2 +1)/(C*eta - D*eta^2)
                end
            end
            # println("K_g:",K_g)
            # println("eta1:",eta1)
            # println("eta:",eta[k,:,t+1])

            for s=1:Numb_Services
                ### Sub1 ###
                if(jpadmm)
                    T_cmp[k,s], f[k,s,:,t+1], E_cmp[k,s]  = Solving_isub_prob1(K_g,D_n,s,f[k,:,:,t],y[k,:], jpadmm) # parallel update in primal JPADMM
                else
                    T_cmp[k,s], f[k,s,:,t+1], E_cmp[k,s]  = Solving_isub_prob1(K_g,D_n,s,f[k,:,:,t+1],y[k,:],jpadmm) # trick for cyclic update in primal ADMM
                end

                ### Sub2 ###
                T_com[k,s], w1[k,s,:,t+1], tau[k,s,:], E_com[k,s] = Solving_isub_prob2(K_g,capacity,s, v[k,:,:], z[k,:,t])
            end

            # ### Sub2 ###
            # T_com[k,:], w[k,:], tau[k,:,:], E_com[k,:] = Solving_sub_prob2(eta,capacity)

            # ### Sub3 ###
            # eta[k,:], Obj[k]  = Solving_sub_prob3(T_cmp[k,:],E_cmp[k,:],T_com[k,:],E_com[k,:])

            ### Dual Update ###
            for n=1:NumbDevs
                z[k,n,t+1]      = 1/Numb_Services * (sum(w1[k,:,n,t+1] + (1/RHO2)*v[k,:,n]))
                r1[k,n,t]   = sum(f[k,:,n,t+1]) - f_max[n]*1e-9
                y[k,n]      = y[k,n] + RHO1*r1[k,n,t]
                r2[k,:,n,t] = w1[k,:,n,t+1] - z[k,n,t+1]*ones(Numb_Services)
                v[k,:,n]    = v[k,:,n] + RHO2*r2[k,:,n,t]
            end

            Obj_E[k,:],Obj_T[k,:], Obj[k,t+1] = compute_obj(f[k,:,:,t+1], z[k,:,t+1], eta[k,:,t+1], D_n, capacity)
            # println("Primal Residual1:",r1)
            # println("Primal Residual2:",r2)

            if ( (norm(f[k,:, :,t+1] - f[k,:, :,t]) <= stop_epsilon1) && (norm(z[k, :,t+1] - z[k, :,t]) <= stop_epsilon2) ||
                 ((early_stop==true) &&  (norm(r1[k, :,t]) <= stop_epsilon3) && (norm(r2[k,:, :,t]) <= stop_epsilon3)) )
                println("\neta:", eta[k, :, t+1])
                println("z:", z[k,1:5,t+1])
                println("f:", f[k,:,1:5,t+1])
                # println("Obj:", Obj[k])
                println("Obj1:", Obj[k,t])
                stop_k[k]= t
                break
            end

        end
   end
   return Obj, r1, r2, eta, z, w1, f, stop_k
end



function main()
    jpadmm=true
    #Generate data
    dist_list, gain_list, capacity, D_n = mobile_gen()
    Obj1, Obj_E, Obj_T, eta1, w1, f1, stop1, Heuristic_Obj, Heuristic_Obj_E, Heuristic_Obj_T,Heuristic_eta,
    Heuristic_Obj1, Heuristic_Obj_E1, Heuristic_Obj_T1,Heuristic_eta1 = BCD(dist_list, gain_list, capacity, D_n)

    Obj2, r1, r2, eta2, w2, ws2, f2, stop2 = ADMM(dist_list, gain_list, capacity, D_n,jpadmm)

    ### Global ###
    rs_Obj, rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, rs_eta, rs_w, rs_f = Solving_global_prob(D_n,capacity)
    rs_T_cmp1, rs_E_cmp1, rs_T_com1, rs_E_com1, rs_eta1 = pareto(D_n,capacity)
    rs_T_cmp2, rs_E_cmp2, rs_T_com2, rs_E_com2, rs_eta2 = priority(D_n,capacity)

    save_result(rs_eta, rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, rs_Obj, rs_w, rs_f, Obj1, Obj2, r1, r2,
    eta1, eta2, w1, w2, ws2, f1, f2, stop1, stop2, Obj_E, Obj_T, Heuristic_Obj, Heuristic_Obj_E, Heuristic_Obj_T, Heuristic_eta,
    Heuristic_Obj1, Heuristic_Obj_E1, Heuristic_Obj_T1, Heuristic_eta1, rs_T_cmp1, rs_E_cmp1, rs_T_com1, rs_E_com1, rs_eta1,
    rs_T_cmp2, rs_E_cmp2, rs_T_com2, rs_E_com2, rs_eta2)
    plot_sub3(rs_eta, Obj_E, Obj_T)

    println("Heuristic_Obj:",Heuristic_Obj[1])
    println("Heuristic_Obj1:",Heuristic_Obj1[1])
    plot_convergence(Obj1, Obj2, rs_Obj, r1, r2, eta1, eta2, rs_eta, w1, w2, ws2, rs_w, f1, f2, rs_f, stop1, stop2)
    # plot_convergence1(Obj1, Obj2, rs_Obj, r1, r2, eta1, eta2, rs_eta, w1, w2, ws2, rs_w, f1, f2, rs_f, stop1, stop2)
    plot_comparison(rs_Obj, Obj_E[1,:], Obj_T[1,:], Heuristic_Obj[1,:],Heuristic_Obj_E[1,:], Heuristic_Obj_T[1,:],
    Heuristic_Obj1[1,:],Heuristic_Obj_E1[1,:], Heuristic_Obj_T1[1,:], rs_eta, Heuristic_eta, Heuristic_eta1)
    plot_data_distance(D_n,dist_list)
    plot_pareto(rs_T_cmp1, rs_E_cmp1, rs_T_com1, rs_E_com1, rs_eta1)
    plot_priority(rs_T_cmp2, rs_E_cmp2, rs_T_com2, rs_E_com2, rs_eta2)

end

function run_multiple_time()
    stop1 =zeros(Int64,Numb_kaps,NUM_SIM)
    stop2 =zeros(Int64,Numb_kaps,NUM_SIM)
    stop3 =zeros(Int64,Numb_kaps,NUM_SIM)
    stop4 =zeros(Int64,Numb_kaps,NUM_SIM)
    rs_Objs  =zeros(3,NUM_SIM)

    global REUSED_TRAFFIC = false
    for i =1:NUM_SIM
    # for i =1:1
        #Generate data
        println("---- Simulation ",i)
        dist_list, gain_list, capacity, D_n = mobile_gen()
        _, _, _, _, _, _, stop1[:,i], _, _, _ = BCD(dist_list, gain_list, capacity, D_n)
        Objs, _, _, _, _, _, _, stop2[:,i] = ADMM(dist_list, gain_list, capacity, D_n,false) #mi-ADMM

        rs_Objs[1,i]= Objs[stop2[1,i]]
        Objs, _, _, _, _, _, _, stop3[:,i] = ADMM(dist_list, gain_list, capacity, D_n,true)  #JP-miADMM
        rs_Objs[2,i]= Objs[stop3[1,i]]
        Objs, _, _, _, _, _, _, stop4[:,i] = ADMM(dist_list, gain_list, capacity, D_n,true,true)  #JP-miADMM, early_stop
        rs_Objs[3,i]= Objs[stop4[1,i]]
    end
    println("Stop1:",stop1)
    println("Stop2:",stop2)
    println("Stop3:",stop3)
    println("Stop4:",stop4)
    println("rs_Objs:",rs_Objs)
    save_result_iteration(stop1,stop2,stop3,stop4,rs_Objs)

end


 # Result for Section IV: Closed-Form solution in the paper (5 devs)
if READ_RESULT
    global REUSED_TRAFFIC = true
    dist_list, _, _, D_n = mobile_gen()

    rs_eta, rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, rs_Obj, rs_w, rs_f, Obj1, Obj2, r1, r2, eta1, eta2, w1, w2,
    ws2, f1, f2, stop1, stop2, Obj_E, Obj_T, Heuristic_Obj, Heuristic_Obj_E, Heuristic_Obj_T, Heuristic_eta,
    Heuristic_Obj1, Heuristic_Obj_E1, Heuristic_Obj_T1, Heuristic_eta1, rs_T_cmp1, rs_E_cmp1, rs_T_com1, rs_E_com1, rs_eta1,
    rs_T_cmp2, rs_E_cmp2, rs_T_com2, rs_E_com2, rs_eta2 = read_result(string("result",NumbDevs,".h5"))
    plot_sub3(rs_eta, Obj_E, Obj_T)

    println("Heuristic_Obj:",Heuristic_Obj[1])
    println("Heuristic_Obj1:",Heuristic_Obj1[1])
    println("rs_Obj:",rs_Obj[1])
    plot_convergence(Obj1, Obj2, rs_Obj, r1, r2, eta1, eta2, rs_eta, w1, w2, ws2, rs_w, f1, f2, rs_f, stop1, stop2)
    # plot_convergence1(Obj1, Obj2, rs_Obj, r1, r2, eta1, eta2, rs_eta, w1, w2, ws2, rs_w, f1, f2, rs_f, stop1, stop2)
    plot_comparison(rs_Obj, Obj_E[1,:], Obj_T[1,:], Heuristic_Obj[1,:],Heuristic_Obj_E[1,:], Heuristic_Obj_T[1,:],
    Heuristic_Obj1[1,:],Heuristic_Obj_E1[1,:], Heuristic_Obj_T1[1,:], rs_eta, Heuristic_eta, Heuristic_eta1)
    plot_data_distance(D_n,dist_list)
    plot_pareto(rs_T_cmp1, rs_E_cmp1, rs_T_com1, rs_E_com1, rs_eta1)
    plot_priority(rs_T_cmp2, rs_E_cmp2, rs_T_com2, rs_E_com2, rs_eta2)

    println("Distances:",dist_list)
    println("D_n:",D_n[:,26])
    println("D_n:",D_n)
    println("rs_w:",findmax(rs_w))
    println("rs_w_high:",rs_w[26],rs_w[24],rs_w[38])
    println("rs_f:",rs_f)
    println("rs_eta:",rs_eta)
    println("rs_T_cmp:",rs_T_cmp)
    println("rs_T_com:",rs_T_com)
    println("stop2:",stop2)

else
    if(NUM_SIM>1)
        run_multiple_time()
    else
        main()
    end
end
