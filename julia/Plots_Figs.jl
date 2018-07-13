# using Plots
using PyPlot
using PyCall
@pyimport matplotlib.patches as patch

# fig_size = (7.,5.1)
fig_size = (6.,4.3)
fig_size1 = (5.5,4.)
label_fontsize = 18-1.5
legend_fontsize = label_fontsize - 3
patterns = ["","."]

label_fontsize1 = label_fontsize
marker_size=6
l_width=1.2
# stride = convert(Int, max_iters/10) +1
stride = 6

colors=["m","b","coral","g","k","r"]
algs = ["PBCD", "Consensus_BCD2", "JP-ADMM", "JP-ADMM_BCD4","IpOpt Solver","Exhaustive Search"]
markers = ["x","o",">","^", ".","s"]

folder = string("figs//")

function kappa_finding_sub1(f1)
    global kaps_draw = zeros(3)
    global kaps_draw_idx = zeros(Int32,3)
    # UEs_min = Numb_kaps * ones(NumbDevs, Numb_kaps)
    min_UEs1 = Numb_kaps
    min_UEs2 = 1
    max_UEs  = Numb_kaps

    for n =1:NumbDevs
        UEs_min  = find(a->abs(a-f_min[n]*1e-9)<5e-4, f1[:,n])
        UEs_max  = find(a->abs(a-f_max[n]*1e-9)<5e-4, f1[:,n])

        if size(UEs_min)[1] > 0
            min_UEs1 = min(min_UEs1, maximum(UEs_min))
            min_UEs2 = max(min_UEs2, maximum(UEs_min))
        end
        if size(UEs_max)[1] > 0
            max_UEs  = min(max_UEs, minimum(UEs_max))
        end
    end

    kaps_draw[1] = kaps[min_UEs1]
    kaps_draw[2] = kaps[min_UEs2+1]
    kaps_draw[3] = kaps[max_UEs]
    kaps_draw_idx[1] = min_UEs1
    kaps_draw_idx[2] = min_UEs2+1
    kaps_draw_idx[3] = max_UEs
    println("kaps_thresh1: ", kaps_draw)
end

function kappa_finding_sub2(p1)
    global kaps_draw2 = zeros(2)
    global kaps_draw_idx2 = zeros(Int32,2)
    # UEs_min = Numb_kaps * ones(NumbDevs, Numb_kaps)
    min_UEs1 = Numb_kaps
    max_UEs  = 1

    for n =1:NumbDevs
        UEs_min  = find(a->abs(a-Ptx_Min)<1e-4, p1[:,n])
        UEs_max  = find(a->abs(a-Ptx_Max)<1e-4, p1[:,n])

        if size(UEs_min)[1] > 0
            min_UEs1 = min(min_UEs1, maximum(UEs_min))
        end
        if size(UEs_max)[1] > 0
            max_UEs  = max(max_UEs, minimum(UEs_max))
        end
    end

    kaps_draw2[1] = kaps[min_UEs1]
    kaps_draw2[2] = kaps[max_UEs]
    kaps_draw_idx2[1] = min_UEs1
    kaps_draw_idx2[2] = max_UEs
    println("kaps_thresh2: ", kaps_draw2)
end

function plot_sub1_T(T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3)
    clf()
    cfig = figure(1,figsize=fig_size)
    ax = cfig[:add_subplot](1,1,1)
    ax[:tick_params]("both",labelsize=legend_fontsize-1)
    # plot(kaps,T_cmp,color=colors[1],linestyle="-",linewidth=l_width,label="Solver")
    plot(kaps,T_cmp1 + 0.02*maximum(T_cmp1),color=colors[6],linestyle="-",linewidth=l_width,label="\$T_{cmp}^*\$")
    plot(kaps,Tcmp_N1,color=colors[2],linestyle="--",linewidth=l_width+0.2,label="\$T_{\\mathcal{N}_1}\$")
    plot(kaps,Tcmp_N2,color=colors[4],linestyle="-.",linewidth=l_width+0.2,label="\$T_{\\mathcal{N}_2}\$")
    plot(kaps,Tcmp_N3,color=colors[5],linestyle="-",linewidth=l_width+0.2,label="\$T_{\\mathcal{N}_3}\$")

    r1 = patch.Rectangle([0,0],kaps_draw[1],1.05*maximum(T_cmp1), alpha=0.07,fc="k",ec="blue",linewidth=.7)
    r2 = patch.Rectangle([kaps_draw[1],0],kaps_draw[2] - kaps_draw[1],T_cmp1[kaps_draw_idx[1]]+ 0.02*maximum(T_cmp1), alpha=0.12,fc="k",ec="blue",linewidth=.7)
    r3 = patch.Rectangle([kaps_draw[2],0],kaps_draw[3] - kaps_draw[2],T_cmp1[kaps_draw_idx[2]]+ 0.02*maximum(T_cmp1), alpha=0.16,fc="k",ec="blue",linewidth=.7)
    r4 = patch.Rectangle([kaps_draw[3],0],maximum(kaps)- kaps_draw[3],T_cmp1[kaps_draw_idx[3]]+ 0.02*maximum(T_cmp1), alpha=0.2,fc="k",ec="blue",linewidth=.7)
    ax[:add_patch](r1)
    ax[:add_patch](r2)
    ax[:add_patch](r3)
    ax[:add_patch](r4)

    annotate("a", xy=[kaps_draw[1]/3;(T_cmp1[1]+1.)/2.], xycoords="data",size=19)
    annotate("b", xy=[kaps_draw[1] + (kaps_draw[2] - kaps_draw[1])/7; (T_cmp1[kaps_draw_idx[1]]+0.5)/2], xycoords="data",size=19)
    annotate("c", xy=[kaps_draw[2] + (kaps_draw[3] - kaps_draw[2])/30;(T_cmp1[kaps_draw_idx[2]]+0.5)/2], xycoords="data",size=19)
    annotate("d", xy=[kaps_draw[3] + (maximum(kaps)- kaps_draw[3])/10;(T_cmp1[kaps_draw_idx[3]]+0.5)/3.5], xycoords="data",size=19)

    legend(loc="best",fontsize=legend_fontsize+2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+2)
    xscale("log")
    ylabel("\$T_{cmp}\$ (sec)",fontsize=label_fontsize1+1)
    ylim(0, 1.065*T_cmp1[1])
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_T.pdf"))
end

function plot_sub1_N(N1, N2, N3)
    clf()
    cfig = figure(2,figsize=fig_size)
    ax = cfig[:add_subplot](1,1,1)
    ax[:tick_params]("both",labelsize=legend_fontsize-1)
    # plot(kaps,T_cmp,color=colors[1],linestyle="-",linewidth=l_width,label="Solver")
    step(kaps,N1,color=colors[4],linestyle="-",linewidth=l_width,label="\$\\mathcal{N}_1\$", where="post", marker=markers[2], markersize=marker_size, markevery=11)
    step(kaps,N2,color=colors[3],linestyle="-",linewidth=l_width,label="\$\\mathcal{N}_2\$", where="post", marker=markers[3], markersize=marker_size, markevery=11)
    step(kaps,N3,color=colors[2],linestyle="-",linewidth=l_width,label="\$\\mathcal{N}_3\$", where="post", marker=markers[6], markersize=marker_size, markevery=11)

    r1 = patch.Rectangle([0,0],kaps_draw[1],NumbDevs, alpha=0.07,fc="k",ec="blue",linewidth=.7)
    r2 = patch.Rectangle([kaps_draw[1],0],kaps_draw[2] - kaps_draw[1],NumbDevs, alpha=0.12,fc="k",ec="blue",linewidth=.7)
    r3 = patch.Rectangle([kaps_draw[2],0],kaps_draw[3] - kaps_draw[2],NumbDevs, alpha=0.16,fc="k",ec="blue",linewidth=.7)
    r4 = patch.Rectangle([kaps_draw[3],0],maximum(kaps)- kaps_draw[3],NumbDevs, alpha=0.2,fc="k",ec="blue",linewidth=.7)
    ax[:add_patch](r1)
    ax[:add_patch](r2)
    ax[:add_patch](r3)
    ax[:add_patch](r4)

    annotate("a", xy=[kaps_draw[1]/3; NumbDevs/2.], xycoords="data",size=19)
    annotate("b", xy=[kaps_draw[1] + (kaps_draw[2] - kaps_draw[1])/10;NumbDevs/2.], xycoords="data",size=19)
    annotate("c", xy=[kaps_draw[2] + (kaps_draw[3] - kaps_draw[2])/30;NumbDevs/2.], xycoords="data",size=19)
    annotate("d", xy=[kaps_draw[3] + (maximum(kaps)- kaps_draw[3])/14;NumbDevs/2.], xycoords="data",size=19)

    legend(loc=1,fontsize=legend_fontsize)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+2)
    xscale("log")
    ylabel("Three subsets by Alg.1",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_N.pdf"))
end

function plot_sub1_f(f1)
    kappa_finding_sub1(f1)

    clf()
    cfig = figure(3,figsize=fig_size)
    ax = cfig[:add_subplot](1,1,1)
    ax[:tick_params]("both",labelsize=legend_fontsize-1)
    plot(kaps,f_min[1]*ones(Numb_kaps)*1e-9,linestyle=":",color=colors[6])

    if (HETEROGENEOUS == 0) # Homogeneous
        plot(kaps,f_max[1]*ones(Numb_kaps)*1e-9,linestyle="--",color=colors[6])
    end

    for n = 1:5
        if (HETEROGENEOUS > 0)  & (abs(f_max[n]*1e-9 - maximum(f1[:,n])) < 1e-3)
            plot(kaps,f_max[n]*ones(Numb_kaps)*1e-9,linestyle="--",color=colors[n])
            # plot(kaps,f_min[n]*ones(Numb_kaps)*1e-9,linestyle=":",color=colors[n])
        end

        plot(kaps,f1[:,n],color=colors[n],linestyle="-",linewidth=l_width,label=string("UE ",n))
    end

    r1 = patch.Rectangle([0,0],kaps_draw[1],minimum(f1)+0.1, alpha=0.07,fc="k",ec="blue",linewidth=.7)
    r2 = patch.Rectangle([kaps_draw[1],0],kaps_draw[2] - kaps_draw[1],minimum(f1)+ 0.3, alpha=0.12,fc="k",ec="blue",linewidth=.7)
    r3 = patch.Rectangle([kaps_draw[2],0],kaps_draw[3] - kaps_draw[2],maximum(f1), alpha=0.16,fc="k",ec="blue",linewidth=.7)
    r4 = patch.Rectangle([kaps_draw[3],0],maximum(kaps)- kaps_draw[3],maximum(f1) + 0.15, alpha=0.2,fc="k",ec="blue",linewidth=.7)
    ax[:add_patch](r1)
    ax[:add_patch](r2)
    ax[:add_patch](r3)
    ax[:add_patch](r4)

    annotate("a", xy=[kaps_draw[1]/3;(minimum(f1)+0.1)/1.5], xycoords="data",size=19)
    annotate("b", xy=[kaps_draw[1] + (kaps_draw[2] - kaps_draw[1])/7;(minimum(f1)+ 0.3)/2], xycoords="data",size=19)
    annotate("c", xy=[kaps_draw[2] + (kaps_draw[3] - kaps_draw[2])/30;maximum(f1)/2], xycoords="data",size=19)
    annotate("d", xy=[kaps_draw[3] + (maximum(kaps)- kaps_draw[3])/14;(maximum(f1) + 0.15)/2], xycoords="data",size=19)

    # axvline(x=kaps_draw[1])
    # axvline(x=kaps_draw[2])
    # axvline(x=kaps_draw[3])

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+2)
    xscale("log")
    ylim(0,maximum(f1) + 0.2 )
    ylabel("f (GHz)",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_f.pdf"))
end

function plot_sub2_tau(tau1)
    clf()
    cfig = figure(4,figsize=fig_size)
    ax = cfig[:add_subplot](1,1,1)
    ax[:tick_params]("both",labelsize=legend_fontsize+1)

    for n = 1:5
        plot(kaps,tau1[:,n], color=colors[n], linestyle="-",linewidth=l_width,label=string("UE ",n))
    end

    max_tau = maximum(tau1[1,:])

    r1 = patch.Rectangle([0,0],kaps_draw2[1], 1.1*max_tau, alpha=0.09,fc="k",ec="blue",linewidth=.7)
    r2 = patch.Rectangle([kaps_draw2[1],0],kaps_draw2[2] - kaps_draw2[1],maximum(tau1[kaps_draw_idx2[1],:]), alpha=0.14,fc="k",ec="blue",linewidth=.7)
    r3 = patch.Rectangle([kaps_draw2[2],0],maximum(kaps)- kaps_draw2[2],maximum(tau1[kaps_draw_idx2[2],:]), alpha=0.2,fc="k",ec="blue",linewidth=.7)
    ax[:add_patch](r1)
    ax[:add_patch](r2)
    ax[:add_patch](r3)

    annotate("a", xy=[kaps_draw2[1]/7;(1.1*max_tau)/2], xycoords="data",size=19)
    annotate("b", xy=[kaps_draw2[1] + (kaps_draw2[2] - kaps_draw2[1])/50;maximum(tau1[kaps_draw_idx2[1],:])/2], xycoords="data",size=19)
    annotate("c", xy=[kaps_draw2[2] + (maximum(kaps)- kaps_draw2[2])/10;maximum(tau1[kaps_draw_idx2[2],:])/2.5], xycoords="data",size=19)

    legend(loc="best",fontsize=legend_fontsize-1)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    # yscale("log")
    xscale("log")
    ylim(0, 1.15*max_tau)
    ylabel("\$\\tau_n\$ (sec)",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub2_Tau.pdf"))
end

function plot_sub2_p(p1)
    kappa_finding_sub2(p1)

    clf()
    cfig = figure(5,figsize=fig_size)
    ax = cfig[:add_subplot](1,1,1)
    ax[:tick_params]("both",labelsize=legend_fontsize+1)

    plot(kaps,Ptx_Max*ones(Numb_kaps),linestyle=":",color=colors[6])
    plot(kaps,Ptx_Min*ones(Numb_kaps),linestyle=":",color=colors[6])
    for n = 1:5
        plot(kaps,p1[:,n],color=colors[n],linestyle="-",linewidth=l_width,label=string("UE ",n))
    end

    r1 = patch.Rectangle([0,0],kaps_draw2[1],minimum(p1)+0.1, alpha=0.09,fc="k",ec="blue",linewidth=.7)
    r2 = patch.Rectangle([kaps_draw2[1],0],kaps_draw2[2] - kaps_draw2[1],maximum(p1), alpha=0.14,fc="k",ec="blue",linewidth=.7)
    r3 = patch.Rectangle([kaps_draw2[2],0],maximum(kaps)- kaps_draw2[2],maximum(p1) + 0.1, alpha=0.2,fc="k",ec="blue",linewidth=.7)
    ax[:add_patch](r1)
    ax[:add_patch](r2)
    ax[:add_patch](r3)

    annotate("a", xy=[kaps_draw2[1]/7;(minimum(p1)+0.1)/1.5], xycoords="data",size=19)
    annotate("b", xy=[kaps_draw2[1] + (kaps_draw2[2] - kaps_draw2[1])/50;maximum(p1)/2], xycoords="data",size=19)
    annotate("c", xy=[kaps_draw2[2] + (maximum(kaps)- kaps_draw2[2])/10;(maximum(p1) + 0.1)/2], xycoords="data",size=19)

    legend(loc=2,fontsize=legend_fontsize-1)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    xscale("log")
    ylim(0,maximum(p1) + 0.15)
    ylabel("p (Watt)",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub2_p.pdf"))
end

function plot_sub3_cvx(Theta1, Obj1, T_cmp1, E_cmp1, T_com1, E_com1)
    clf()
    cfig = figure(6,figsize=fig_size1)
    ax = cfig[:add_subplot](1,1,1)
    ax[:tick_params]("both",labelsize=legend_fontsize+2)

    x = collect(1.e-5:0.001:0.99)
    obj   = zeros(size(x)[1])
    glob_cost_iter = zeros(size(x)[1])
    glob_numb_iter = zeros(size(x)[1])
    id = 34
    println("Convex for kappa: ",  kaps[id])
    for i=1:size(x)[1]
        obj[i] = 1/(1 - x[i])* (E_com1[id] - log(x[i])*E_cmp1[id] + kaps[id] * (T_com1[id] - log(x[i])*T_cmp1[id]))
        glob_cost_iter[i] = E_com1[id] - log(x[i])*E_cmp1[id] + kaps[id] * (T_com1[id] - log(x[i])*T_cmp1[id])
        glob_numb_iter[i] = 1/(1 - x[i])
        # obj[i]   = obj_E[i] + obj_T[i]
    end
    plot(x, obj,linestyle="-",color="k", label=string("SUB3 Obj: \$\\kappa\$ =", kaps[id]))
    plot(x, glob_cost_iter,linestyle="--",color=colors[2], label=string("\$E_{glob} + \\kappa * T_{glob}\$"))
    plot(x, glob_numb_iter,linestyle="--",color=colors[3], label=string("\$ K(\\theta)\$"))
    # println(x)
    plot(Theta1[id], Obj1[id],color="r", marker=markers[2], markersize=marker_size)

    legend(loc="best",fontsize=legend_fontsize+6)
    xlabel("\$\\theta\$",fontsize=label_fontsize1+3)
    # ylabel("Objective",fontsize=label_fontsize1+1)
    yscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub3_obj.pdf"))
    println("Theta: ", minimum(Theta1), " - ", maximum(Theta1))
end

function plot_sub3_kappa_theta(Theta, d_eta)
    clf()
    cfig = figure(10,figsize=fig_size1)
    ax = cfig[:add_subplot](1,1,1)
    ax[:tick_params]("both",labelsize=legend_fontsize+3.5)
    # plot(Numb_devs, Objs_E[:,11],linestyle="--",color=colors[1],marker=markers[1], markersize=marker_size, label=string("\$\\kappa\$ =", kaps[11]))
    # plot(kaps, 1./d_eta,linestyle="-",color=colors[3],label="Heterogeneous: \$\\eta\$")
    plot(kaps, Theta,linestyle="-",color=colors[2])
    # plot(kaps, Theta1,linestyle="-",color=colors[3],label="Homogeneous:\$\\kappa\$")
    # plot(kaps, 1./d_eta,linestyle="-",color=colors[3],label="Homogeneous")

    # legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+7)
    ylabel("\$\\theta^*\$",fontsize=label_fontsize1+4)
    xscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"sub3_kappa_theta.pdf"))
end

function plot_sub3_equation(Theta, d_eta)
    clf()
    x = collect(1.e-6:0.001:0.999)

    cfig = figure(7,figsize=fig_size1)
    ax = cfig[:add_subplot](1,1,1)
    ax[:tick_params]("both",labelsize=legend_fontsize+2)
    id1 = 18
    id2 = 28
    plot(x,d_eta[id1]*ones(size(x)),linestyle="-",color="b",label=string("\$\\kappa\$ = ",kaps[id1]))
    plot(x,d_eta[id2]*ones(size(x)),linestyle="-",color="g",label=string("\$\\kappa\$ = ",kaps[id2]))
    # hlines(y=d_eta[20],xmin=0, xmax=Theta[20], linestyle=":",color="k", zorder=1)
    # hlines(y=d_eta[30],xmin=0, xmax=Theta[30], linestyle=":",color="k", zorder=1)
    vlines(x=Theta[id1],ymin=0, ymax=d_eta[id1], linestyle=":",color="k", zorder=2)
    vlines(x=Theta[id2],ymin=0, ymax=d_eta[id2], linestyle=":",color="k", zorder=2)

    plot(x, 1./x + log.(x),linestyle="-",color=colors[6], label="\$\\log(e^{1/\\theta} \\theta)\$")
    # for k = 1:Numb_kaps
    #     plot(x,d_eta[k]*ones(size(x)),linestyle=":",color="k")
    # end

    annotate(string("(",round(Theta[id1],3),", 1/",round(1/d_eta[id1],3),")"), xy=[Theta[id1];1.05*d_eta[id1]], xycoords="data",size=18)
    annotate(string("(",round(Theta[id2],3),", 1/",round(1/d_eta[id2],3),")"), xy=[Theta[id2];1.2*d_eta[id2]], xycoords="data",size=18)
    scatter(Theta[id1], d_eta[id1],color="k")
    scatter(Theta[id2], d_eta[id2],color="k")

    legend(loc="best",fontsize=legend_fontsize+6)
    xlim(0, 0.3)
    ylim(0.98,maximum(d_eta)+0.1*maximum(d_eta))
    xlabel("\$\\theta\$",fontsize=label_fontsize1+3)
    ylabel("\$1/\\eta\$",fontsize=label_fontsize1+3)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub3_eq.pdf"))
end

function plot_numerical_pareto(Theta1, T_cmp1, E_cmp1, T_com1, E_com1)
    clf()
    figure(9,figsize=fig_size)

    E_obj   = zeros(Numb_kaps)
    T_obj   = zeros(Numb_kaps)

    for i=1:Numb_kaps
        E_obj[i] = 1/(1 - Theta1[i])* (E_com1[i] - log(Theta1[i])*E_cmp1[i])
        T_obj[i] = 1/(1 - Theta1[i])* (T_com1[i] - log(Theta1[i])*T_cmp1[i])
    end
    scatter(E_obj, T_obj)

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("Energy Cost",fontsize=label_fontsize1+1)
    ylabel("Time Cost",fontsize=label_fontsize1+1)
    # yscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"pareto.pdf"))
end

# function plot_scale_result()
#     Sims = size(Numb_devs)[1]
#     Thetas = zeros(Sims, Numb_kaps)
#     Objs   = zeros(Sims, Numb_kaps)
#     Objs_E = zeros(Sims, Numb_kaps)
#     Objs_T = zeros(Sims, Numb_kaps)
#
#     for i = 1:Sims
#         Thetas[i,:], Objs[i,:], Objs_E[i,:], Objs_T[i,:], T_cmp1, E_cmp1, T_com1, E_com1,
#         N1, N2, N3, f1, tau1, p1,
#         d_eta = read_result(string("result",Numb_devs[i],".h5"))
#     end
#
#     # clf()
#     # figure(8,figsize=fig_size)
#     # plot(Numb_devs, Objs[:,11],linestyle="--",color=colors[1],marker=markers[1], markersize=marker_size, label=string("Objective: \$\\kappa\$ =", kaps[11]))
#     # plot(Numb_devs, Objs[:,15],linestyle="--",color=colors[2],marker=markers[2], markersize=marker_size, label=string("Objective: \$\\kappa\$ =", kaps[15]))
#     # plot(Numb_devs, Objs[:,19],linestyle="--",color=colors[3],marker=markers[3], markersize=marker_size, label=string("Objective: \$\\kappa\$ =", kaps[19]))
#     # plot(Numb_devs, Objs[:,23],linestyle="--",color=colors[4],marker=markers[4], markersize=marker_size, label=string("Objective: \$\\kappa\$ =", kaps[23]))
#     #
#     # legend(loc="best",fontsize=legend_fontsize-2)
#     # xlabel("Number of Devs",fontsize=label_fontsize1+1)
#     # # ylabel("Objective",fontsize=label_fontsize1+1)
#     # # yscale("log")
#     # tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
#     # savefig(string(folder,"Scale_obj.pdf"))
#     #
#     # clf()
#     # figure(9,figsize=fig_size)
#     # plot(Numb_devs, Thetas[:,11],linestyle="--",color=colors[1],marker=markers[1], markersize=marker_size, label=string("\$\\theta\$: \$\\kappa\$ =", kaps[11]))
#     # plot(Numb_devs, Thetas[:,15],linestyle="--",color=colors[2],marker=markers[2], markersize=marker_size, label=string("\$\\theta\$: \$\\kappa\$ =", kaps[15]))
#     # plot(Numb_devs, Thetas[:,19],linestyle="--",color=colors[3],marker=markers[3], markersize=marker_size, label=string("\$\\theta\$: \$\\kappa\$ =", kaps[19]))
#     # plot(Numb_devs, Thetas[:,23],linestyle="--",color=colors[4],marker=markers[4], markersize=marker_size, label=string("\$\\theta\$: \$\\kappa\$ =", kaps[23]))
#     # # plot(Numb_devs, Thetas[:,id],linestyle="-",color="k", label="\$\\theta\$")
#     #
#     # legend(loc="best",fontsize=legend_fontsize-2)
#     # xlabel("Number of Devs",fontsize=label_fontsize1+1)
#     # # ylabel("Objective",fontsize=label_fontsize1+1)
#     # # yscale("log")
#     # tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
#     # savefig(string(folder,"Scale_theta.pdf"))
#
#     clf()
#     figure(10,figsize=fig_size)
#     # plot(Numb_devs, Objs_E[:,11],linestyle="--",color=colors[1],marker=markers[1], markersize=marker_size, label=string("\$\\kappa\$ =", kaps[11]))
#     plot(Numb_devs, Objs_E[:,15],linestyle="--",color=colors[2],marker=markers[2], markersize=marker_size, label=string("\$\\kappa\$ =", kaps[15]))
#     plot(Numb_devs, Objs_E[:,19],linestyle="--",color=colors[3],marker=markers[3], markersize=marker_size, label=string("\$\\kappa\$ =", kaps[19]))
#     plot(Numb_devs, Objs_E[:,23],linestyle="--",color=colors[4],marker=markers[4], markersize=marker_size, label=string("\$\\kappa\$ =", kaps[23]))
#
#     legend(loc="best",fontsize=legend_fontsize-2)
#     xlabel("Number of Devs",fontsize=label_fontsize1+1)
#     ylabel("Energy cost",fontsize=label_fontsize1+1)
#     # yscale("log")
#     tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
#     savefig(string(folder,"Scale_obj_E.pdf"))
#
#     clf()
#     figure(11,figsize=fig_size)
#     # plot(Numb_devs, Objs_T[:,11],linestyle="--",color=colors[1],marker=markers[1], markersize=marker_size, label=string("\$\\kappa\$ =", kaps[11]))
#     plot(Numb_devs, Objs_T[:,15],linestyle="--",color=colors[2],marker=markers[2], markersize=marker_size, label=string("\$\\kappa\$ =", kaps[15]))
#     plot(Numb_devs, Objs_T[:,19],linestyle="--",color=colors[3],marker=markers[3], markersize=marker_size, label=string("\$\\kappa\$ =", kaps[19]))
#     plot(Numb_devs, Objs_T[:,23],linestyle="--",color=colors[4],marker=markers[4], markersize=marker_size, label=string("\$\\kappa\$ =", kaps[23]))
#
#     legend(loc="best",fontsize=legend_fontsize-2)
#     xlabel("Number of Devs",fontsize=label_fontsize1+1)
#     ylabel("Time cost",fontsize=label_fontsize1+1)
#     # yscale("log")
#     tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
#     savefig(string(folder,"Scale_obj_T.pdf"))
#
# end

function save_result(Theta1, Obj1, Obj_E, Obj_T, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta)
    h5open(string("result",NumbDevs,".h5"), "w") do file
        # write(file,"kaps", kaps)
        write(file,"Theta1", Theta1)
        write(file,"Obj1", Obj1)
        write(file,"Obj_E", Obj_E)
        write(file,"Obj_T", Obj_T)
        write(file,"T_cmp1", T_cmp1)
        write(file,"Tcmp_N1", Tcmp_N1)
        write(file,"Tcmp_N2", Tcmp_N2)
        write(file,"Tcmp_N3", Tcmp_N3)
        write(file,"E_cmp1", E_cmp1)
        write(file,"T_com1", T_com1)
        write(file,"E_com1", E_com1)
        write(file,"N1", N1)
        write(file,"N2", N2)
        write(file,"N3", N3)
        write(file,"f1", f1)
        write(file,"tau1", tau1)
        write(file,"p1", p1)
        write(file,"d_eta", d_eta)
    end
end

function read_result(filename)
    h5open(filename, "r") do file
        # kaps =read(file,"kaps")
        Theta1 = read(file,"Theta1")
        Obj1  = read(file,"Obj1")
        Obj_E = read(file,"Obj_E")
        Obj_T = read(file,"Obj_T")
        T_cmp1 = read(file,"T_cmp1")
        Tcmp_N1 = read(file,"Tcmp_N1")
        Tcmp_N2 = read(file,"Tcmp_N2")
        Tcmp_N3 = read(file,"Tcmp_N3")
        E_cmp1 = read(file,"E_cmp1")
        T_com1 = read(file,"T_com1")
        E_com1 = read(file,"E_com1")
        N1 = read(file,"N1")
        N2 = read(file,"N2")
        N3 = read(file,"N3")
        f1 = read(file,"f1")
        tau1 = read(file,"tau1")
        p1 = read(file,"p1")
        d_eta = read(file,"d_eta")
        return Theta1, Obj1, Obj_E, Obj_T, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta
    end
end
