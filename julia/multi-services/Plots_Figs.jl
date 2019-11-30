# using Plots
ENV["MPLBACKEND"]="qt5agg" #tkagg, qt5agg, gtkagg
using PyPlot
using PyCall
using Printf
@pyimport matplotlib.patches as patch
# include("plt_FL.jl")

# fig_size = (7.,5.1)
# fig_size = (6.,4.3)
fig_size = (6.,4.)
fig_size1 = (5.5,4.)
label_fontsize = 18-1.5
legend_fontsize = label_fontsize - 3
patterns = ["", ".","/"]

label_fontsize1 = label_fontsize
marker_size=6
l_width=1.2
# stride = convert(Int, max_iters/10) +1
stride = 6

colors=["m","b","coral","g","k","r"]
colors1 = ["dodgerblue", "mediumseagreen", "coral"]
colors2 = ["skyblue", "mediumaquamarine", "sandybrown"]
colors3 = ["dodgerblue", "mediumseagreen", "coral", "crimson", "violet"]
algs = ["PBCD", "Consensus_BCD2", "JP-ADMM", "JP-ADMM_BCD4","IpOpt Solver","Exhaustive Search"]
markers = ["x","o",">","^", "s","."]

folder = string("figs//")

function kappa_finding_sub1(f1)
    global kaps_draw = zeros(3)
    global kaps_draw_idx = zeros(Int32,3)
    # UEs_min = Numb_kaps * ones(NumbDevs, Numb_kaps)
    min_UEs1 = Numb_kaps
    min_UEs2 = 1
    max_UEs  = Numb_kaps

    for n =1:NumbDevs
        UEs_min  = findall(a->abs(a-f_min[n]*1e-9)<1e-3, f1[:,n])
        UEs_max  = findall(a->abs(a-f_max[n]*1e-9)<1e-3, f1[:,n])

        if size(UEs_min)[1] > 0
            min_UEs1 = min(min_UEs1, maximum(UEs_min))
            min_UEs2 = max(min_UEs2, maximum(UEs_min))
        end
        if size(UEs_max)[1] > 0
            max_UEs  = min(max_UEs, minimum(UEs_max))
        end
    end

    kaps_draw[1] = kaps[min_UEs1]
    kaps_draw[2] = kaps[min_UEs2]
    kaps_draw[3] = kaps[max_UEs]
    kaps_draw_idx[1] = min_UEs1
    kaps_draw_idx[2] = min_UEs2
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
        UEs_min  = findall(a->abs(a-Ptx_Min)<1e-4, p1[:,n])
        UEs_max  = findall(a->abs(a-Ptx_Max)<1e-4, p1[:,n])

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

function plot_sub1_T(T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3)
    clf()
    cfig = figure(1,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize-1)
    plot(kaps,T_cmp1 .+ 0.02*maximum(T_cmp1),color=colors[6],linestyle="-",linewidth=l_width+0.3,label="\$T_{cmp}^*\$")
    # plot(kaps,T_cmp + 0.02*maximum(T_cmp1),color="gold",linestyle=":",linewidth=l_width+0.3,label="Solver")
    plot(kaps,Tcmp_N1,color=colors[2],linestyle="--",linewidth=l_width+0.2,label="\$T_{\\mathcal{N}_1}\$")
    plot(kaps,Tcmp_N2,color=colors[4],linestyle="-.",linewidth=l_width+0.2,label="\$T_{\\mathcal{N}_2}\$")
    plot(kaps,Tcmp_N3,color=colors[5],linestyle="-",linewidth=l_width+0.2,label="\$T_{\\mathcal{N}_3}\$")

    r1 = patch.Rectangle([0,0],kaps_draw[1],1.05*maximum(T_cmp1), alpha=0.07,fc="k",ec="blue",linewidth=.7)
    r2 = patch.Rectangle([kaps_draw[1],0],kaps_draw[2] - kaps_draw[1],T_cmp1[kaps_draw_idx[1]]+ 0.02*maximum(T_cmp1), alpha=0.12,fc="k",ec="blue",linewidth=.7)
    r3 = patch.Rectangle([kaps_draw[2],0],kaps_draw[3] - kaps_draw[2],T_cmp1[kaps_draw_idx[2]]+ 0.02*maximum(T_cmp1), alpha=0.16,fc="k",ec="blue",linewidth=.7)
    r4 = patch.Rectangle([kaps_draw[3],0],maximum(kaps)- kaps_draw[3],T_cmp1[kaps_draw_idx[3]]+ 0.02*maximum(T_cmp1), alpha=0.2,fc="k",ec="blue",linewidth=.7)
    ax.add_patch(r1)
    ax.add_patch(r2)
    ax.add_patch(r3)
    ax.add_patch(r4)

    annotate("a", xy=[kaps_draw[1]/2;(T_cmp1[1]+0.07*maximum(T_cmp1))], xycoords="data",size=19)
    annotate("b", xy=[kaps_draw[1] + (kaps_draw[2] - kaps_draw[1])/7; (T_cmp1[kaps_draw_idx[1]]+0.04*maximum(T_cmp1))], xycoords="data",size=19)
    annotate("c", xy=[kaps_draw[2] + (kaps_draw[3] - kaps_draw[2])/6;(T_cmp1[kaps_draw_idx[2]]+0.04*maximum(T_cmp1))], xycoords="data",size=19)
    annotate("d", xy=[kaps_draw[3] + (maximum(kaps)- kaps_draw[3])/450;(T_cmp1[kaps_draw_idx[3]]+0.04*maximum(T_cmp1))], xycoords="data",size=19)

    legend(loc="best",fontsize=legend_fontsize+2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+2)
    xscale("log")
    ylabel("\$T_{cmp}\$ (sec)",fontsize=label_fontsize1+1)
    xlim(1e-3, 1e1)
    ylim(0, 1.15*T_cmp1[1])
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_T.pdf"))
end

function plot_sub1_N(N1, N2, N3)
    clf()
    N3[kaps_draw_idx[end]] = N3[kaps_draw_idx[end]]+1
    cfig = figure(2,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize-1)
    # plot(kaps,T_cmp,color=colors[1],linestyle="-",linewidth=l_width,label="Solver")
    step(kaps,N1,color=colors[4],linestyle="-",linewidth=l_width,label="\$\\mathcal{N}_1\$", where="post", marker=markers[2], markersize=marker_size, markevery=5)
    step(kaps,N2,color=colors[3],linestyle="-",linewidth=l_width,label="\$\\mathcal{N}_2\$", where="pre", marker=markers[3], markersize=marker_size, markevery=5)
    step(kaps,N3,color=colors[2],linestyle="-",linewidth=l_width,label="\$\\mathcal{N}_3\$", where="pre", marker=markers[5], markersize=marker_size, markevery=5)

    r1 = patch.Rectangle([0,0],kaps_draw[1],NumbDevs, alpha=0.07,fc="k",ec="blue",linewidth=.7)
    r2 = patch.Rectangle([kaps_draw[1],0],kaps_draw[2] - kaps_draw[1],NumbDevs, alpha=0.12,fc="k",ec="blue",linewidth=.7)
    r3 = patch.Rectangle([kaps_draw[2],0],kaps_draw[3] - kaps_draw[2],NumbDevs, alpha=0.16,fc="k",ec="blue",linewidth=.7)
    r4 = patch.Rectangle([kaps_draw[3],0],maximum(kaps)- kaps_draw[3],NumbDevs, alpha=0.2,fc="k",ec="blue",linewidth=.7)
    ax.add_patch(r1)
    ax.add_patch(r2)
    ax.add_patch(r3)
    ax.add_patch(r4)

    annotate("a", xy=[kaps_draw[1]/2; NumbDevs/2.], xycoords="data",size=19)
    annotate("b", xy=[kaps_draw[1] + (kaps_draw[2] - kaps_draw[1])/13;NumbDevs/2.], xycoords="data",size=19)
    annotate("c", xy=[kaps_draw[2] + (kaps_draw[3] - kaps_draw[2])/6;NumbDevs/2.], xycoords="data",size=19)
    annotate("d", xy=[kaps_draw[3] + (maximum(kaps)- kaps_draw[3])/450;NumbDevs/2.], xycoords="data",size=19)

    legend(loc=1,fontsize=legend_fontsize)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+2)
    xscale("log")
    xlim(1e-3, 1e1)
    ylabel("Number of elements",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_N.pdf"))
end

function plot_sub1_f(f1)
    kappa_finding_sub1(f1)

    clf()
    cfig = figure(3,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize-1)
    plot(kaps,f_min[1]*ones(Numb_kaps)*1e-9,linestyle=":",color=colors[6])

    if (HETEROGENEOUS == 0) # Homogeneous
        plot(kaps,f_max[1]*ones(Numb_kaps)*1e-9,linestyle="--",color=colors[6])
    end

    for n = 1:5
        if (HETEROGENEOUS > 0)  & (abs(f_max[n]*1e-9 - maximum(f1[:,n])) < 1e-3)
            plot(kaps,f_max[n]*ones(Numb_kaps)*1e-9,linestyle="--",color=colors[n])
            # plot(kaps,f_min[n]*ones(Numb_kaps)*1e-9,linestyle=":",color=colors[n])
        end

        plot(kaps,f1[:,n],color=colors[n],linestyle="-",linewidth=l_width,marker=markers[n], markersize=marker_size-1, markevery=3, label=string("UE ",n))
    end

    r1 = patch.Rectangle([0,0],kaps_draw[1],minimum(f1)+0.1, alpha=0.07,fc="k",ec="blue",linewidth=.7)
    r2 = patch.Rectangle([kaps_draw[1],0],kaps_draw[2] - kaps_draw[1],minimum(f1)+ 0.3, alpha=0.12,fc="k",ec="blue",linewidth=.7)
    r3 = patch.Rectangle([kaps_draw[2],0],kaps_draw[3] - kaps_draw[2],maximum(f1), alpha=0.16,fc="k",ec="blue",linewidth=.7)
    r4 = patch.Rectangle([kaps_draw[3],0],maximum(kaps)- kaps_draw[3],maximum(f1) + 0.15, alpha=0.2,fc="k",ec="blue",linewidth=.7)
    ax.add_patch(r1)
    ax.add_patch(r2)
    ax.add_patch(r3)
    ax.add_patch(r4)

    annotate("a", xy=[kaps_draw[1]/2;(minimum(f1)+0.15)], xycoords="data",size=19)
    annotate("b", xy=[kaps_draw[1] + (kaps_draw[2] - kaps_draw[1])/7;(minimum(f1)+ 0.35)], xycoords="data",size=19)
    annotate("c", xy=[kaps_draw[2] + (kaps_draw[3] - kaps_draw[2])/6;maximum(f1)+ 0.05], xycoords="data",size=19)
    annotate("d", xy=[kaps_draw[3] + (maximum(kaps)- kaps_draw[3])/450;(maximum(f1) + 0.2)], xycoords="data",size=19)

    # axvline(x=kaps_draw[1])
    # axvline(x=kaps_draw[2])
    # axvline(x=kaps_draw[3])

    legend(loc=2,fontsize=legend_fontsize-2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+2)
    xscale("log")
    xlim(1e-3, 1e1)
    ylim(0.2,maximum(f1) + 0.35 )
    ylabel("f (GHz)",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_f.pdf"))
end

function plot_sub2_tau(tau1)
    clf()
    cfig = figure(4,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize+1)

    for n = 1:5
        plot(kaps,tau1[:,n], color=colors[n], linestyle="-",linewidth=l_width, marker=markers[n], markersize=marker_size-1, markevery=3, label=string("UE ",n))
    end

    max_tau = maximum(tau1[1,:])

    # r1 = patch.Rectangle([0,0],kaps_draw2[1], 1.1*max_tau, alpha=0.09,fc="k",ec="blue",linewidth=.7)
    # r2 = patch.Rectangle([kaps_draw2[1],0],kaps_draw2[2] - kaps_draw2[1],maximum(tau1[kaps_draw_idx2[1],:]), alpha=0.14,fc="k",ec="blue",linewidth=.7)
    # r3 = patch.Rectangle([kaps_draw2[2],0],maximum(kaps)- kaps_draw2[2],maximum(tau1[kaps_draw_idx2[2],:]), alpha=0.2,fc="k",ec="blue",linewidth=.7)
    # ax.add_patch(r1)
    # ax.add_patch(r2)
    # ax.add_patch(r3)
    #
    # annotate("a", xy=[kaps_draw2[1]/7;(1.1*max_tau)/2], xycoords="data",size=19)
    # annotate("b", xy=[kaps_draw2[1] + (kaps_draw2[2] - kaps_draw2[1])/50;maximum(tau1[kaps_draw_idx2[1],:])/2], xycoords="data",size=19)
    # annotate("c", xy=[kaps_draw2[2] + (maximum(kaps)- kaps_draw2[2])/10;maximum(tau1[kaps_draw_idx2[2],:])/2.5], xycoords="data",size=19)

    legend(loc="best",fontsize=legend_fontsize-1)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    # yscale("log")
    xscale("log")
    ylim(0, 1.15*max_tau)
    xlim(1e-3, 1e1)
    ylabel("\$\\tau_n\$ (sec)",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub2_Tau.pdf"))
end

function plot_sub2_p(p1)
    kappa_finding_sub2(p1)

    clf()
    cfig = figure(5,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize+1)

    plot(kaps,Ptx_Max*ones(Numb_kaps),linestyle=":",color=colors[6])
    plot(kaps,Ptx_Min*ones(Numb_kaps),linestyle=":",color=colors[6])
    for n = 1:5
        plot(kaps,p1[:,n],color=colors[n],linestyle="-",linewidth=l_width, marker=markers[n], markersize=marker_size-1, markevery=3, label=string("UE ",n))
    end

    # r1 = patch.Rectangle([0,0],kaps_draw2[1],minimum(p1)+0.1, alpha=0.09,fc="k",ec="blue",linewidth=.7)
    # r2 = patch.Rectangle([kaps_draw2[1],0],kaps_draw2[2] - kaps_draw2[1],maximum(p1), alpha=0.14,fc="k",ec="blue",linewidth=.7)
    # r3 = patch.Rectangle([kaps_draw2[2],0],maximum(kaps)- kaps_draw2[2],maximum(p1) + 0.1, alpha=0.2,fc="k",ec="blue",linewidth=.7)
    # ax.add_patch(r1)
    # ax.add_patch(r2)
    # ax.add_patch(r3)
    #
    # annotate("a", xy=[kaps_draw2[1]/7;(minimum(p1)+0.1)/1.5], xycoords="data",size=19)
    # annotate("b", xy=[kaps_draw2[1] + (kaps_draw2[2] - kaps_draw2[1])/50;maximum(p1)/2], xycoords="data",size=19)
    # annotate("c", xy=[kaps_draw2[2] + (maximum(kaps)- kaps_draw2[2])/10;(maximum(p1) + 0.1)/2], xycoords="data",size=19)

    legend(loc=2,fontsize=legend_fontsize-1)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    xscale("log")
    ylim(0.1,maximum(p1) + 0.15)
    xlim(1e-3, 1e1)
    ylabel("p (Watt)",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub2_p.pdf"))
end

function plot_sub3_cvx(Theta1, Obj1, T_cmp1, E_cmp1, T_com1, E_com1)
    clf()
    cfig = figure(6,figsize=fig_size1)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize+2)

    x = collect(1.e-5:0.001:0.99)
    obj   = zeros(size(x)[1])
    glob_cost_iter = zeros(size(x)[1])
    glob_numb_iter = zeros(size(x)[1])
    # id = 36
    id = 32
    # println("Convex for kappa: ",  kaps[id])
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
    println("At kap = ", kaps[id],", Optimal Theta:", Theta1[id])

    plot_FL_MINIST_5users(kaps[id])
    obj = zeros(size(local_iters))
    obj1 = zeros(size(local_iters))

    for i=1:size(local_iters)[1]
        obj[i] = global_iters[i]/Normalized_Global* (E_com1[id] + local_iters[i]/Normalized_Local *E_cmp1[id] + kaps[id] *
        (T_com1[id] + local_iters[i]/Normalized_Local *T_cmp1[id]))
        obj1[i] = 1/(1 - theory_theta2[i])* (E_com1[id] - log(theory_theta2[i]) *E_cmp1[id] + kaps[id] *
        (T_com1[id] - log(theory_theta2[i]) *T_cmp1[id]))
    end
    println("Thoery_obj:",obj1)
    println("True_obj:",obj)

end


function plot_sub3_kappa_theta(Theta, d_eta)
    clf()
    cfig = figure(10,figsize=fig_size1)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize+3.5)
    # plot(Numb_devs, Objs_E[:,11],linestyle="--",color=colors[1],marker=markers[1], markersize=marker_size, label=string("\$\\kappa\$ =", kaps[11]))
    plot(kaps, 1 ./d_eta,linestyle="--",color=colors[3],label="\$\\eta\$")
    plot(kaps, Theta,linestyle="-",color=colors[2],label="\$\\theta^*\$")
    # plot(kaps, Theta1,linestyle="-",color=colors[3],label="Homogeneous:\$\\kappa\$")
    # plot(kaps, 1 ./d_eta,linestyle="-",color=colors[3],label="Homogeneous")

    legend(loc="best",fontsize=legend_fontsize+2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+7)
    ylabel("\$\\theta^*\$ and \$\\eta\$",fontsize=label_fontsize1+4)
    xscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"sub3_kappa_theta.pdf"))

    println("kaps: ", kaps[24], ", theta:", @sprintf("%.3f",Theta[24]), ", eta: ", @sprintf("%.3f",1/d_eta[24]) )
    println("kaps: ", kaps[end], ", theta:", @sprintf("%.3f",Theta[end]), ", eta: ", @sprintf("%.3f",1/d_eta[end]) )
end

function plot_sub3_equation(Theta, d_eta)
    clf()
    x = collect(1.e-6:0.001:0.999)

    cfig = figure(7,figsize=fig_size1)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize+2)
    id1 = 24
    id2 = 32
    plot(x,d_eta[id1]*ones(size(x)),linestyle="-",color="b",label=string("\$\\kappa\$ = ",kaps[id1]))
    plot(x,d_eta[id2]*ones(size(x)),linestyle="-",color="g",label=string("\$\\kappa\$ = ",kaps[id2]))
    # hlines(y=d_eta[20],xmin=0, xmax=Theta[20], linestyle=":",color="k", zorder=1)
    # hlines(y=d_eta[30],xmin=0, xmax=Theta[30], linestyle=":",color="k", zorder=1)
    vlines(x=Theta[id1],ymin=0, ymax=d_eta[id1], linestyle=":",color="k", zorder=2)
    vlines(x=Theta[id2],ymin=0, ymax=d_eta[id2], linestyle=":",color="k", zorder=2)

    plot(x, 1 ./x + log.(x),linestyle="-",color=colors[6], label="\$\\log(e^{1/\\theta} \\theta)\$")
    # for k = 1:Numb_kaps
    #     plot(x,d_eta[k]*ones(size(x)),linestyle=":",color="k")
    # end

    annotate(string("(",@sprintf("%.3f",Theta[id1]),", 1/",@sprintf("%.3f",1/d_eta[id1]),")"), xy=[Theta[id1];1.05*d_eta[id1]], xycoords="data",size=18)
    annotate(string("(",@sprintf("%.3f",Theta[id2]),", 1/",@sprintf("%.3f",1/d_eta[id2]),")"), xy=[0.9*Theta[id2];1.1*d_eta[id2]], xycoords="data",size=18)
    scatter(Theta[id1], d_eta[id1],color="k")
    scatter(Theta[id2], d_eta[id2],color="k")

    legend(loc="best",fontsize=legend_fontsize+6)
    xlim(0, 1.)
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

function plot_convergence(Obj1, Obj2, rs_Obj, r1, r2, Theta1, Theta2, rs_Theta, w1, w2, ws2, rs_w, f1, f2, rs_f, stop1, stop2)
    Stride=8
    clf()
    t = Int(max(stop1[1],stop2[1]))
    figure(2,figsize=fig_size)
    plot(rs_Obj*ones(t),label="Solver",linestyle=":", marker=markers[1], markersize=marker_size,markevery=Stride)
    plot(Obj2[1:stop2[1]],label="miADMM",linestyle="--", marker=markers[2], markersize=marker_size,markevery=Stride)
    plot(Obj1[1:stop1[1]],label="BCD",linestyle="-", marker=markers[3], markersize=marker_size,markevery=4)

    # ylim(6.4,6.8)
    legend(loc="best",fontsize=legend_fontsize-1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Convergence_Obj.pdf"))

    clf()
    figure(3,figsize=fig_size)
    for n=1:NumbDevs
        plot(r1[1,n,1:stop2[1]])
    end
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Convergence_Residual1.pdf"))

    clf()
    figure(4,figsize=fig_size)
    for s=1:Numb_Services
        for n=1:NumbDevs
            plot(r2[1,s,n,1:stop2[1]])
        end
    end
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Convergence_Residual2.pdf"))

    clf()
    figure(5,figsize=fig_size)

    for s=1:Numb_Services
        if(s==1)
            plot(rs_Theta[s]*ones(t),label="Solver",linestyle=":", marker=markers[1], markersize=marker_size,markevery=Stride)
            plot(Theta1[1,s,1:stop1[1]],label="BCD",linestyle="-", marker=markers[3], markersize=marker_size,markevery=4)

        else
            plot(rs_Theta[s]*ones(t),linestyle=":", marker=markers[1], markersize=marker_size,markevery=Stride)
            plot(Theta1[1,s,1:stop1[1]],linestyle="-", marker=markers[3], markersize=marker_size,markevery=4)
        end
        plot(Theta2[1,s,1:stop2[1]],label=string("miADMM:S",s),linestyle="--", marker=markers[2], markersize=marker_size,markevery=Stride)
    end

    legend(loc="best",fontsize=legend_fontsize-1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Convergence_Theta.pdf"))

    clf()
    figure(6,figsize=fig_size)
    for n=1:NumbDevs
        if(n==1)
            plot(rs_w[n]*ones(t),label="Solver",linestyle=":")
            plot(w2[1,n,1:stop2[1]],label="miADMM",linestyle="-")
        else
            plot(rs_w[n]*ones(t),linestyle=":")
            plot(w2[1,n,1:stop2[1]], linestyle="-")
            # for s=1:Numb_Services
            #     plot(ws2[1,s,n,1:stop2[1]],label="miADMM",linestyle="--")
            # end
        end
        # plot(w1[1,n,1:stop1[1]],label="BCD",linestyle="--")
    end
    legend(loc="best",fontsize=legend_fontsize-1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Convergence_w.pdf"))

    clf()
    figure(7,figsize=fig_size)
    for n=1:NumbDevs
        if(n==1)
            for s=1:Numb_Services
                if(s==1)
                    plot(rs_f[s,n]*ones(t),label="Solver",linestyle=":", marker=markers[1], markersize=marker_size,markevery=Stride)
                else
                    plot(rs_f[s,n]*ones(t),linestyle=":", marker=markers[1], markersize=marker_size,markevery=Stride)
                end

                plot(f2[1,s,n,1:stop2[1]],label=string("miADMM:S",s),linestyle="--", marker=markers[2], markersize=marker_size,markevery=Stride)
                # plot(Theta1[1,s,1:stop1[1]],label="BCD",linestyle="--")
            end
        end
    end
    legend(loc="best",fontsize=legend_fontsize-1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Convergence_f.pdf"))


end

# FOR BCD Only
function plot_convergence1(Obj1, Obj2, rs_Obj, r1, r2, Theta1, Theta2, rs_Theta, w1, w2, ws2, rs_w, f1, f2, rs_f, stop1, stop2)
    fig_size = (4.5,5)
    clf()
    t = Int(stop1[1])
    figure(1,figsize=fig_size)
    plot(rs_Obj*ones(t),label="Solver",linestyle=":")
    plot(Obj1[1:stop1[1]],label="BCD",linestyle="--")
    xlabel("Iteration",fontsize=label_fontsize1)
    ylabel("Total Cost",fontsize=label_fontsize1)
    legend(loc="best",fontsize=legend_fontsize-1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Convergence_Obj.pdf"))

    clf()
    figure(4,figsize=fig_size)
    for s=1:Numb_Services
        if(s==1)
            plot(rs_Theta[s]*ones(t),label="Solver",linestyle=":")
        else
            plot(rs_Theta[s]*ones(t),linestyle=":")
        end

        plot(Theta1[1,s,1:stop1[1]],label=string("BCD-Service",s),linestyle="-")
    end
    xlabel("Iteration",fontsize=label_fontsize1)
    ylabel("\$\\theta\$",fontsize=label_fontsize1)
    legend(loc="best",fontsize=legend_fontsize-1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Convergence_Theta.pdf"))

    clf()
    figure(5,figsize=fig_size)
    for n=1:NumbDevs
        if(n==1)
            plot(rs_w[n]*ones(t),label="Solver",linestyle=":")
            plot(w1[1,n,1:stop1[1]],label="BCD",linestyle="-")
        else
            plot(rs_w[n]*ones(t),linestyle=":")
            plot(w1[1,n,1:stop1[1]],linestyle="-")
        end
    end
    xlabel("Iteration",fontsize=label_fontsize1)
    ylabel("\$w\$",fontsize=label_fontsize1)
    legend(loc="best",fontsize=legend_fontsize-1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Convergence_w.pdf"))

    clf()
    figure(6,figsize=fig_size)
    for n=1:NumbDevs
        # for s=1:Numb_Services
        if(n==1)
            plot(rs_f[1,n]*ones(t),label="Solver",linestyle=":")
            plot(f1[1,1,n,1:stop1[1]],label="BCD",linestyle="-")
        else
            plot(rs_f[1,n]*ones(t),linestyle=":")
            plot(f1[1,1,n,1:stop1[1]],linestyle="-")
        end
        # end
    end
    xlabel("Iteration",fontsize=label_fontsize1)
    ylabel("\$f\$",fontsize=label_fontsize1)
    legend(loc="best",fontsize=legend_fontsize-1)
    savefig(string(folder,"Convergence_f.pdf"))
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

end

function plot_comparison(rs_Obj, Obj_E, Obj_T, Heuristic_Obj,Heuristic_Obj_E, Heuristic_Obj_T)
    # println("HERE")
    # fig_size = (8,4)
    clf()
    figure(20,figsize=fig_size)
    alg_labels = ["Heuristic","Optimal"]
    Service_labels = ["Service 1", "Service 2", "Service 3"]
    idx = [1:1:Numb_Services;]
    width = 0.25
    barlist1 = bar(idx.-width/2, Heuristic_Obj_E, width, color=colors1[1], alpha=.8, hatch=patterns[1], label=alg_labels[1])
    barlist2 = bar(idx.+width/2, Obj_E, width, color=colors1[3], alpha=.8, hatch=patterns[2], label=alg_labels[2])

    xticks(idx, Service_labels, fontsize=label_fontsize1)
    ylabel("Energy Consumption", fontsize=label_fontsize1)
    legend(loc="best",fontsize=legend_fontsize-1)
    savefig(string(folder,"E_Comparison.pdf"))
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

    clf()
    figure(21,figsize=fig_size)
    alg_labels = ["Heuristic","Optimal"]
    Service_labels = ["Service 1", "Service 2", "Service 3"]
    idx = [1:1:Numb_Services;]
    width = 0.25
    barlist1 = bar(idx.-width/2, Heuristic_Obj_T, width, color=colors1[1], alpha=.8, hatch=patterns[1], label=alg_labels[1])
    barlist2 = bar(idx.+width/2, Obj_T, width, color=colors1[3], alpha=.8, hatch=patterns[2], label=alg_labels[2])

    xticks(idx, Service_labels, fontsize=label_fontsize1)
    ylabel("Total Time", fontsize=label_fontsize1)
    legend(loc="best",fontsize=legend_fontsize-1)
    savefig(string(folder,"T_Comparison.pdf"))
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

    clf()
    figure(22,figsize=fig_size)
    alg_labels = ["Heuristic","Optimal"]
    Service_labels = ["Service 1", "Service 2", "Service 3"]
    idx = [1:1:Numb_Services;]
    width = 0.25
    barlist1 = bar(idx.-width/2, Heuristic_Obj_E + Heuristic_Obj_T, width, color=colors1[1], alpha=.8, hatch=patterns[1], label=alg_labels[1])
    barlist2 = bar(idx.+width/2, Obj_E + Obj_T, width, color=colors1[3], alpha=.8, hatch=patterns[2], label=alg_labels[2])

    xticks(idx, Service_labels, fontsize=label_fontsize1)
    ylabel("Total cost", fontsize=label_fontsize1)
    legend(loc="best",fontsize=legend_fontsize-1)
    savefig(string(folder,"Service_Cost_Comparison.pdf"))
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

end

function save_result(Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta)
    h5open(string("result",NumbDevs,".h5"), "w") do file
        # write(file,"kaps", kaps)
        write(file,"Theta1", Theta1)
        write(file,"Obj1", Obj1)
        write(file,"Obj_E", Obj_E)
        write(file,"Obj_T", Obj_T)
        write(file,"T_cmp", T_cmp)
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
        T_cmp = read(file,"T_cmp")
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
        return Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta
    end
end


function save_result_iteration(stop1,stop2,stop3)
    h5open(string("result_iter_",NUM_SIM,".h5"), "w") do file
        # write(file,"kaps", kaps)
        write(file,"stop1", stop1)
        write(file,"stop2", stop2)
        write(file,"stop3", stop3)
    end
end

function read_result_iteration()
    h5open(string("result_iter_",NUM_SIM,".h5"), "r") do file
        # write(file,"kaps", kaps)
        stop1 = read(file,"stop1")
        stop2  = read(file,"stop2")
        stop3 = read(file,"stop3")
    end
    return stop1, stop2, stop3
end
