include("Setting.jl")
include("TrafficGen.jl")
include("Solving1.jl")
include("Plots_Figs.jl")
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
    Numb_kaps = size(kaps)[1]
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

   save_result(Theta1, Obj1, Obj_E, Obj_T, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta)

   plot_sub1_T(T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3)
   plot_sub1_N(N1, N2, N3)

   if(NumbDevs <=5)
       plot_sub1_f(f1)
       plot_sub2_tau(tau1)
       plot_sub2_p(p1)
   end

   plot_sub3_equation(d_eta)
   plot_sub3_cvx(Theta1, Obj1, T_cmp1, E_cmp1, T_com1, E_com1)
   plot_sub3_kappa_theta(Theta1, d_eta)
   plot_numerical_pareto(Theta1, T_cmp1, E_cmp1, T_com1, E_com1)
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

if READ_RESULT
    Theta1, Obj1, Obj_E, Obj_T, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3,
    E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1,
    d_eta = read_result(string("result5.h5"))
    #
    # Theta1, Obj1, Obj_E1, Obj_T1, T_cmp11, E_cmp11, T_com11, E_com11,
    # N11, N21, N31, f11, tau11, p11,
    # d_eta1 = read_result(string("result5_homo.h5"))

    plot_sub1_T(T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3)
    plot_sub1_N(N1, N2, N3)

    if(NumbDevs <=5)
        plot_sub1_f(f1)
        plot_sub2_tau(tau1)
        plot_sub2_p(p1)
    end

    plot_sub3_equation(d_eta)
    plot_sub3_cvx(Theta1, Obj1, T_cmp1, E_cmp1, T_com1, E_com1)
    plot_sub3_kappa_theta(Theta1, d_eta)
    plot_numerical_pareto(Theta1, T_cmp1, E_cmp1, T_com1, E_com1)
else
    main()
end
