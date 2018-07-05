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
    # kaps = [5e-5, 8e-5, 1e-4, 3e-4, 5e-4, 7e-4, 1e-3, 3e-3, 5e-3, 7e-3, 1e-2, 3e-2, 5e-2, 7e-2, 1e-1, 0.3, 0.5, 0.7, 1.,3., 5., 7., 1e1, 5e1, 1e2]
    kaps = [5e-5, 8e-5, 1e-4, 3e-4, 5e-4, 7e-4, 1e-3, 3e-3, 5e-3, 7e-3, 1e-2, 3e-2, 5e-2, 7e-2, 1e-1, 0.3, 0.5, 0.7, 1.,3., 5., 7., 1e1, 5e1, 7e1, 1e2]
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

   plot_sub1_T(kaps, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3)
   plot_sub1_N(kaps, N1, N2, N3)
   plot_sub1_f(kaps, f, f1)
   plot_sub2_tau(kaps, tau, tau1)
   plot_sub2_p(kaps, p, p1)

   plot_sub3_equation(kaps, d_eta)
   plot_sub3_cvx(kaps, Theta, Theta1, Obj, Obj1, T_cmp1, E_cmp1, T_com1, E_com1)
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

main()
# main1()
