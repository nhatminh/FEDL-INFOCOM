include("TrafficGen.jl")
include("Solving.jl")

ck_e = 1e-4
function check(list, list1)
    rs = true
    for i =1:size(list)[1]
        rs = rs & (norm(list[i] - list1[i]) < ck_e)
    end

    if(rs)
        return "MATCH"
    else
        return "NOT MATCH"
    end
end

function main()
    dist_list, gain_list, ratios, D_n = mobile_gen()

    ### Sub1 ###
    T_cmp, f, E_cmp     = Solving_sub_prob1(D_n)
    T_cmp1, f1, E_cmp1  = Solving_sub1(D_n)
    println("\n---->> Check Sub1 Solution: ", check([T_cmp, f, E_cmp], [T_cmp1, f1, E_cmp1]))

    # ### Sub2 ###
    # T_com, p, tau, E_com    = Solving_sub_prob2(ratios)
    # T_com1, p1, tau1, E_com1= Solving_sub2(ratios)
    # println("\n---->> Check Sub2 Solution: ", check([T_com, p, tau, E_com], [T_com1, p1, tau1, E_com1]))
    #
    # ### Sub3 ###
    # Theta  = Solving_sub_prob3(T_cmp,E_cmp,T_com,E_com)
    # Theta1 = Solving_sub3(T_cmp,E_cmp,T_com,E_com)
    # println("\n---->> Check Sub3 Solution: ", check([Theta], [Theta1]))
    #
    # ### Global ###
    # rs2 = Solving_global_prob(D_n,ratios)
    # rs  = [T_cmp, E_cmp, T_com, E_com, Theta]
    # println("\n---->> Check Global Solution: ", check(rs, rs2))
end

main()
