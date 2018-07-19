include("Setting.jl")
include("TrafficGen.jl")
include("Solving1.jl")
if(NUMERICAL_RS)
    include("Plots_Figs1.jl")
else
    include("Plots_Figs.jl")
end

ck_e = 1e-2
function check(list, list1)
    rs = true
    for i =1:size(list)[1]
        rs = rs & (norm(list[i] - list1[i]) < ck_e)
        # println(i,":",rs)
    end

    if(rs)
        return "MATCH"
    else
        return "NOT MATCH"
    end
end

function main()
    #Generate data
    dist_list, gain_list, ratios, D_n = mobile_gen()
    T_cmp, T_cmp1   = zeros(Numb_kaps), zeros(Numb_kaps)
    E_cmp, E_cmp1   = zeros(Numb_kaps), zeros(Numb_kaps)
    f, f1           = zeros(Numb_kaps,NumbDevs), zeros(Numb_kaps,NumbDevs)
    N1, N2, N3      = zeros(Numb_kaps), zeros(Numb_kaps), zeros(Numb_kaps)
    Tcmp_N1, Tcmp_N2, Tcmp_N3   = zeros(Numb_kaps), zeros(Numb_kaps), zeros(Numb_kaps)


    T_com, T_com1   = zeros(Numb_kaps),zeros(Numb_kaps)
    E_com, E_com1   = zeros(Numb_kaps), zeros(Numb_kaps)
    p, p1           = zeros(Numb_kaps,NumbDevs), zeros(Numb_kaps,NumbDevs)
    tau, tau1       = zeros(Numb_kaps,NumbDevs), zeros(Numb_kaps,NumbDevs)

    Theta, Theta1   = zeros(Numb_kaps), zeros(Numb_kaps)
    Obj, Obj1       = zeros(Numb_kaps), zeros(Numb_kaps)
    Obj_E, Obj_T    = zeros(Numb_kaps), zeros(Numb_kaps)
    d_eta  = zeros(Numb_kaps)
    global ZoneA = []
    global ZoneB1 = []
    global ZoneB2 = []
    global ZoneC = []
    global ZoneD = []
    # println("Numb_kaps: ", Numb_kaps)
    for s =1:Numb_SIMs
        # for k=1:1
        for k=1:Numb_kaps
            global kappa = kaps[k]
            ### Sub1 ###
            T_cmp[k], f[k,:], E_cmp[k]     = Solving_sub_prob1(D_n[s,:])
            T_cmp1[k], f1[k,:], E_cmp1[k], Tcmp_N1[k], Tcmp_N2[k], Tcmp_N3[k],N1[k], N2[k], N3[k]   = Solving_sub1(D_n[s,:])
            # println("\n---->> Check Sub1 Solution: ", check([T_cmp, f, E_cmp], [T_cmp1, f1, E_cmp1]))

            ### Sub2 ###
            T_com[k], p[k,:], tau[k,:], E_com[k]     = Solving_sub_prob2(ratios[s,:])
            T_com1[k], p1[k,:], tau1[k,:], E_com1[k] = Solving_sub2(ratios[s,:])
            # println("\n---->> Check Sub2 Solution: ", check([T_com, p, tau, E_com], [T_com1, p1, tau1, E_com1]))

            ### Sub3 ###
            Theta[k], Obj[k]  = Solving_sub_prob3(T_cmp[k],E_cmp[k],T_com[k],E_com[k])
            Theta1[k], Obj1[k], Obj_E[k], Obj_T[k], d_eta[k] = Solving_sub3(T_cmp1[k],E_cmp1[k],T_com1[k],E_com1[k])
            # println("\n---->> Check Sub3 Solution: ", check([Theta], [Theta1]))

            # ### Global ###
            # rs2 = Solving_global_prob(D_n[s,:],ratios[s,:])
            # rs  = [T_cmp, E_cmp, T_com, E_com, Theta]
            # println("\n---->> Check Global Solution: ", check(rs, rs2))
        end
   end

   save_result(Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta)

   plot_sub1_f(f1)
   plot_sub1_T(T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3)
   plot_sub1_N(N1, N2, N3)

   plot_sub2_p(p1)
   plot_sub2_tau(tau1)

   plot_sub3_equation(Theta1, d_eta)
   plot_sub3_cvx(Theta1, Obj1, T_cmp1, E_cmp1, T_com1, E_com1)
   plot_sub3_kappa_theta(Theta1, d_eta)
   plot_numerical_pareto(Theta1, T_cmp1, E_cmp1, T_com1, E_com1)
   println("K1: ",minimum(alpha*(f_min.^3)) )
   println("ZoneA: ", ZoneA)
   # println("ZoneB1: ", ZoneB1)
   println("ZoneB2: ", ZoneB2)
   println("ZoneC: ", ZoneC)
   println("ZoneD: ", ZoneD)
end

function main_sub1()
    #Generate data
    dist_list, gain_list, ratios, D_n = mobile_gen_sub1()
    T_cmp1, T_cmp   = zeros(Numb_D,Numb_kaps), zeros(Numb_D,Numb_kaps)
    E_cmp1   = zeros(Numb_D,Numb_kaps)
    f1           = zeros(Numb_D,Numb_kaps,NumbDevs)
    N1, N2, N3      = zeros(Numb_D,Numb_kaps), zeros(Numb_D,Numb_kaps), zeros(Numb_D,Numb_kaps)
    Tcmp_N1, Tcmp_N2, Tcmp_N3   = zeros(Numb_D,Numb_kaps), zeros(Numb_D,Numb_kaps), zeros(Numb_D,Numb_kaps)

    T_com1   = zeros(Numb_D,Numb_kaps)
    E_com1   = zeros(Numb_D,Numb_kaps)
    p1           = zeros(Numb_D,Numb_kaps,NumbDevs)
    tau1       = zeros(Numb_D,Numb_kaps,NumbDevs)

    Theta1   = zeros(Numb_D,Numb_kaps)
    Obj1       = zeros(Numb_D,Numb_kaps)
    Obj_E, Obj_T    = zeros(Numb_D,Numb_kaps), zeros(Numb_D,Numb_kaps)
    d_eta  = zeros(Numb_D,Numb_kaps)

    t_ratios = 1./D_ratios*(f_min[1]/f_max[1])

    tau_max = minimum(S_n/BW./log.(Ptx_Min./ratios[:] .+1))
    tau_min = maximum(S_n/BW./log.(Ptx_Max./ratios[:] .+1))
    tau_ratios = tau_min/tau_max

    for s =1:Numb_D
        for k=1:Numb_kaps
            global kappa = kaps[k]
            ### Sub1 ###
            T_cmp[s,k], _, _    = Solving_sub_prob1(D_n[s,:])
            T_cmp1[s,k], f1[s,k,:], E_cmp1[s,k], Tcmp_N1[s,k], Tcmp_N2[s,k], Tcmp_N3[s,k],N1[s,k], N2[s,k], N3[s,k]   = Solving_sub1(D_n[s,:])

            ### Sub2 ###
            # T_com[s,k], p[s,k,:], tau[s,k,:], E_com[s,k]     = Solving_sub_prob2(ratios[s,:])
            T_com1[s,k], p1[s,k,:], tau1[s,k,:], E_com1[s,k] = Solving_sub2(ratios[:])
            # println("\n---->> Check Sub2 Solution: ", check([T_com, p, tau, E_com], [T_com1, p1, tau1, E_com1]))

            ### Sub3 ###
            # Theta[s,k], Obj[s,k]  = Solving_sub_prob3(T_cmp[k],E_cmp[k],T_com[k],E_com[k])
            Theta1[s,k], Obj1[s,k], Obj_E[s,k], Obj_T[s,k], d_eta[s,k] = Solving_sub3(T_cmp1[s,k],E_cmp1[s,k],T_com1[s,k],E_com1[s,k])
            # println("\n---->> Check Sub3 Solution: ", check([Theta], [Theta1]))
        end
   end
   filename = string("result",NumbDevs,"_sub1.h5")
   save_result(filename,Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta,t_ratios, tau_ratios)

   plot_sub1_f(f1[:,10,:], t_ratios)
   plot_sub1_T(T_cmp[:,10], T_cmp1[:,10], Tcmp_N1[:,10], Tcmp_N2[:,10], Tcmp_N3[:,10], t_ratios)
   plot_sub1_N(N1[:,10], N2[:,10], N3[:,10], t_ratios)

    # plot_sub3_equation(d_eta, t_ratios)
    # plot_sub3_cvx(Theta1, Obj1, T_cmp1, E_cmp1, T_com1, E_com1, t_ratios)
    plot_sub3_kappa_theta(Theta1, d_eta, t_ratios, 1)
    plot_numerical_pareto(Theta1, T_cmp1, E_cmp1, T_com1, E_com1, t_ratios, 1)
    plot_total_cost(Obj1, t_ratios, 1)
end

function main_sub2()
    #Generate data
    dist_list, gain_list, ratios, D_n = mobile_gen_sub2()

    T_cmp1, T_cmp   = zeros(Numb_D,Numb_kaps), zeros(Numb_D,Numb_kaps)
    E_cmp1   = zeros(Numb_Dis,Numb_kaps)
    f1           = zeros(Numb_Dis,Numb_kaps,NumbDevs)
    N1, N2, N3      = zeros(Numb_Dis,Numb_kaps), zeros(Numb_Dis,Numb_kaps), zeros(Numb_Dis,Numb_kaps)
    Tcmp_N1, Tcmp_N2, Tcmp_N3   = zeros(Numb_Dis,Numb_kaps), zeros(Numb_Dis,Numb_kaps), zeros(Numb_Dis,Numb_kaps)

    T_com1   = zeros(Numb_Dis,Numb_kaps)
    E_com1   = zeros(Numb_Dis,Numb_kaps)
    p1       = zeros(Numb_Dis,Numb_kaps,NumbDevs)
    tau1     = zeros(Numb_Dis,Numb_kaps,NumbDevs)

    Theta1   = zeros(Numb_Dis,Numb_kaps)
    Obj1       = zeros(Numb_Dis,Numb_kaps)
    Obj_E, Obj_T    = zeros(Numb_Dis,Numb_kaps), zeros(Numb_Dis,Numb_kaps)
    d_eta  = zeros(Numb_Dis,Numb_kaps)

    global kappa = 1.
    tau_ratios = zeros(Numb_Dis)
    t_ratios = maximum(D_n./f_max)/minimum(D_n./f_min)

    for s =1:Numb_Dis
        tau_max = minimum(S_n/BW./log.(Ptx_Min./ratios[s,:] .+1))
        tau_min = maximum(S_n/BW./log.(Ptx_Max./ratios[s,:] .+1))
        tau_ratios[s] = tau_min/tau_max
   end

   tau_ratios_idx = sortperm(tau_ratios)

   s = 1
   for t in tau_ratios_idx
       for k=1:Numb_kaps
           global kappa = kaps[k]
           ### Sub1 ###
           # T_cmp[s,k], _, _ = Solving_sub_prob1(D_n[s,:])
           T_cmp1[s,k], f1[s,k,:], E_cmp1[s,k], Tcmp_N1[s,k], Tcmp_N2[s,k], Tcmp_N3[s,k],N1[s,k], N2[s,k], N3[s,k]   = Solving_sub1(D_n[:])

           ### Sub2 ###
           # T_com[s,k], p[s,k,:], tau[s,k,:], E_com[s,k]     = Solving_sub_prob2(ratios[s,:])
           # println(size(ratios[t,:])[1])
           T_com1[s,k], p1[s,k,:], tau1[s,k,:], E_com1[s,k] = Solving_sub2(ratios[t,:])
           # println("\n---->> Check Sub2 Solution: ", check([T_com, p, tau, E_com], [T_com1, p1, tau1, E_com1]))

           ### Sub3 ###
           # Theta[s,k], Obj[s,k]  = Solving_sub_prob3(T_cmp[k],E_cmp[k],T_com[k],E_com[k])
           Theta1[s,k], Obj1[s,k], Obj_E[s,k], Obj_T[s,k], d_eta[s,k] = Solving_sub3(T_cmp1[s,k],E_cmp1[s,k],T_com1[s,k],E_com1[s,k])
           # println("\n---->> Check Sub3 Solution: ", check([Theta], [Theta1]))
       end
       s += 1
   end
   tau_ratios_sorted = sort(tau_ratios)
   # println("here3: ",tau_ratios_sorted)
   filename = string("result",NumbDevs,"_sub2.h5")
   save_result(filename,Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta,t_ratios,tau_ratios_sorted)

   plot_sub2_p(p1[:,10,:], tau_ratios_sorted)
   plot_sub2_tau(tau1[:,10,:], tau_ratios_sorted)
   plot_sub2_Tcom(T_com1[:,10],tau_ratios_sorted)
   # plot_sub3_equation(d_eta, tau_ratios_sorted)
   # plot_sub3_cvx(Theta1, Obj1, T_cmp1, E_cmp1, T_com1, E_com1, tau_ratios_sorted)
   plot_sub3_kappa_theta(Theta1, d_eta, tau_ratios_sorted, 2)
   plot_numerical_pareto(Theta1, T_cmp1, E_cmp1, T_com1, E_com1, tau_ratios_sorted, 2)
   plot_total_cost(Obj1, tau_ratios_sorted, 2)
   println("here1: ",tau_ratios)
   println("here2: ",tau_ratios_idx)
   println("here3: ",tau_ratios_sorted)
end

function main1()
    #Generate data
    dist_list, gain_list, ratios, D_n = mobile_gen()
    for s =1:Numb_SIMs
        ### Sub1 ###
        T_cmp, f, E_cmp     = Solving_sub_prob1(D_n[s,:])
        T_cmp1, f1, E_cmp1  = Solving_sub1(D_n[s,:])
        println("\n---->> Check Sub1 Solution: ", check([T_cmp, f, E_cmp], [T_cmp1, f1, E_cmp1]))


        # ### Sub2 ###
        # T_com, p, tau, E_com    = Solving_sub_prob2(ratios[s,:])
        # T_com1, p1, tau1, E_com1= Solving_sub2(ratios[s,:])
        # println("\n---->> Check Sub2 Solution: ", check([T_com, p, tau, E_com], [T_com1, p1, tau1, E_com1]))
        #
        # ### Sub3 ###
        # Theta  = Solving_sub_prob3(T_cmp,E_cmp,T_com,E_com)
        # Theta1 = Solving_sub3(T_cmp,E_cmp,T_com,E_com)
        # println("\n---->> Check Sub3 Solution: ", check([Theta], [Theta1]))
        #
        # ### Global ###
        # rs2 = Solving_global_prob(D_n[s,:],ratios[s,:])
        # rs  = [T_cmp, E_cmp, T_com, E_com, Theta]
        # println("\n---->> Check Global Solution: ", check(rs, rs2))
    end
end

if(NUMERICAL_RS)
    if READ_RESULT
        ### RATIO 1
        dist_list, gain_list, ratios, D_n = mobile_gen_sub1()

        # Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3,
        # E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1,
        # d_eta, levels = read_result(string("result",NumbDevs,"_sub1.h5"))
        Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3,
        E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1,
        d_eta,levels, tau_ratios = read_result(string("result",NumbDevs,"_sub1.h5"))

        id_kap = 10
        println("Figs of kappa = ",kaps[id_kap])
        plot_sub1_f(f1[:,id_kap,:], levels)
        plot_sub1_T(T_cmp[:,id_kap],T_cmp1[:,id_kap], Tcmp_N1[:,id_kap], Tcmp_N2[:,id_kap], Tcmp_N3[:,id_kap], levels)
        plot_sub1_N(N1[:,id_kap], N2[:,id_kap], N3[:,id_kap], levels)
        plot_sub3_kappa_theta(Theta1, d_eta, levels, 1)
        plot_numerical_pareto(Theta1, T_cmp1, E_cmp1, T_com1, E_com1, levels, 1)
        plot_total_cost(Obj1, levels, 1)

        ### RATIO 2
        dist_list, gain_list, ratios, D_n = mobile_gen_sub2()

        Theta1, Obj1, Obj_E, Obj_T, T_cmp2, T_cmp12, Tcmp_N1, Tcmp_N2, Tcmp_N3,
        E_cmp12, T_com12, E_com12, N1, N2, N3, f1, tau1, p1,
        d_eta, t_ratios, levels2 = read_result(string("result",NumbDevs,"_sub2.h5"))
        plot_sub2_p(p1[:,id_kap,:], levels2)
        plot_sub2_tau(tau1[:,id_kap,:], levels2)
        plot_sub2_Tcom(T_com12[:,id_kap],levels2)
        plot_sub3_kappa_theta(Theta1, d_eta, levels1, 2)
        plot_numerical_pareto(Theta1, T_cmp12, E_cmp12, T_com12, E_com12, levels1, 2)
        plot_total_cost(Obj1, levels2, 2)

        plot_ratios(T_cmp1, E_cmp1, T_com1, E_com1, T_cmp12, E_cmp12, T_com12, E_com12, levels, levels2)
        println("L_com_fix: ",tau_ratios)
        println("L_cmp_fix: ",t_ratios)
    else
        main_sub1()
        main_sub2()
    end
elseif READ_RESULT
    dist_list, gain_list, ratios, D_n = mobile_gen()

    Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3,
    E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1,
    d_eta = read_result(string("result",NumbDevs,".h5"))

    # Theta1, Obj1, Obj_E1, Obj_T1, T_cmp11, E_cmp11, T_com11, E_com11,
    # N11, N21, N31, f11, tau11, p11,
    # d_eta1 = read_result(string("result5_homo.h5"))

    plot_sub1_f(f1)
    plot_sub1_T(T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3)
    plot_sub1_N(N1, N2, N3)

    plot_sub2_p(p1)
    plot_sub2_tau(tau1)

    plot_sub3_equation(Theta1, d_eta)
    plot_sub3_cvx(Theta1, Obj1, T_cmp1, E_cmp1, T_com1, E_com1)
    plot_sub3_kappa_theta(Theta1, d_eta)
    plot_numerical_pareto(Theta1, T_cmp1, E_cmp1, T_com1, E_com1)

    # println("K1: ",minimum(alpha*(f_min.^3)) )
    # println("K2: ",sum(alpha*1e27*((C_n.*D_n*1e-9).^3))/maximum(C_n.*D_n./f_min)^3 )
    # println("old K2: ",sum(alpha*1e27*((f_min*1e-9).^3)) )
    # println("K3: ",sum(alpha*1e27*((C_n.*D_n*1e-9).^3))/maximum(C_n.*D_n./f_max)^3 )
    # println("old K3: ",sum(alpha*1e27*((f_max*1e-9).^3)) )

else
    main()
end
