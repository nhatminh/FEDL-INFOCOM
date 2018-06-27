from TrafficGen import *
from Solving import *

if __name__== "__main__":
    dist_list, gain_list, ratios, D_n = mobile_gen()
    # T_com, p, tau, E_com = Solving_sub_prob2(ratios)
    T_cmp, f, E_cmp     = Solving_sub_prob1(D_n)
    # T_com, p, tau, E_com= Solving_sub_prob2(ratios)
    Theta  = Solving_sub_prob4(T_cmp,E_cmp,T_cmp,E_cmp)