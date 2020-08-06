# using Plots
ENV["MPLBACKEND"]="qt5agg" #tkagg, qt5agg, gtkagg
using PyPlot
using PyCall
using Printf
@pyimport matplotlib.patches as patch
# include("plt_FL.jl")

# fig_size = (7.,5.1)
# fig_size = (6.,4.3)
fig_size = (6.,4.1)
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
    ax.tick_params("both",labelsize =legend_fontsize-1)
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
    ax.tick_params("both",labelsize =legend_fontsize-1)
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
    ax.tick_params("both",labelsize =legend_fontsize-1)
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
    ax.tick_params("both",labelsize =legend_fontsize+1)

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
    ax.tick_params("both",labelsize =legend_fontsize+1)

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

function plot_sub3(rs_eta, Obj_E, Obj_T)
    x = collect(1.e-5:0.001:0.33)
    obj   = zeros(size(x)[1])
    glob_cost_iter = zeros(size(x)[1])
    glob_numb_iter = zeros(size(x)[1])
    K_g = zeros(Numb_Services)
    fig_size = (6.5,4.)
    clf()
    cfig = figure(6,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize =legend_fontsize+1)

    for s=1:Numb_Services
        C_s =  Obj_E[s] + kaps[1]*Obj_T[s]
        for i=1:size(x)[1]
            K_g[s] = 2*rho0*As[s]*(Bs[s]*x[i]^2 +1)/(Cs[s]*x[i] - Ds[s]*x[i]^2)
            obj[i] = K_g[s]* C_s
            glob_cost_iter[i] = C_s
            glob_numb_iter[i] = K_g[s]
        end
        Opt_Obj= 2*rho0*As[s]*(Bs[s]*rs_eta[s]^2 +1) /(Cs[s]*rs_eta[s] - Ds[s]*rs_eta[s]^2)*C_s
        plot(x, obj,linestyle="-",color=colors1[s], label=string("Service ",s))
        # plot(x, glob_cost_iter,linestyle="--",color=colors[2], label=string("\$E_{gl} + \\kappa_",s," * T_{gl}\$"))
        # plot(x, glob_numb_iter,linestyle="--",color=colors[3], label=string("\$ K_1(\\eta_",s,")\$"))
        plot(rs_eta[s], Opt_Obj,color=colors1[s], marker=markers[2], markersize=marker_size+3)
    end

    legend(loc=2,fontsize=legend_fontsize+1)
    xlabel("\$\\eta\$",fontsize=label_fontsize1+1)
    ylabel("SUB1-d Objecive",fontsize=label_fontsize1+1)
    yscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    grid(true)
    savefig(string(folder,"Sub3_obj.pdf"))
end

function plot_numerical_pareto(eta1, T_cmp1, E_cmp1, T_com1, E_com1)
    clf()
    figure(9,figsize=fig_size)

    E_obj   = zeros(Numb_kaps)
    T_obj   = zeros(Numb_kaps)

    for i=1:Numb_kaps
        E_obj[i] = 1/(1 - eta1[i])* (E_com1[i] - log(eta1[i])*E_cmp1[i])
        T_obj[i] = 1/(1 - eta1[i])* (T_com1[i] - log(eta1[i])*T_cmp1[i])
    end
    scatter(E_obj, T_obj)

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("Energy Cost",fontsize=label_fontsize1+1)
    ylabel("Time Cost",fontsize=label_fontsize1+1)
    # yscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"pareto.pdf"))
end

function plot_convergence(Obj1, Obj2, rs_Obj, r1, r2, eta1, eta2, rs_eta, w1, w2, ws2, rs_w, f1, f2, rs_f, stop1, stop2)
    Stride=10
    Stride1=2
    clf()
    t = Int(max(stop1[1],stop2[1]))
    figure(2,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize-1)
    plot(rs_Obj*ones(t),label="Solver", color="k",linestyle=":", marker=markers[1], markersize=marker_size,markevery=Stride)
    plot(Obj2[1:stop2[1]],label="Decentralized",linestyle="--", marker=markers[2], markersize=marker_size,markevery=Stride)
    plot(Obj1[1:3],label="Centralized", color="r",linestyle="-", marker=markers[3], markersize=marker_size,markevery=Stride1)

    # ylim(6.4,6.8)
    legend(loc="best",fontsize=legend_fontsize-1)
    xlabel("Iteration", fontsize=label_fontsize)
    ylabel("Total Cost", fontsize=label_fontsize)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    grid(true)
    savefig(string(folder,"Convergence_Obj.pdf"))

    clf()
    figure(3,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize =legend_fontsize-1)
    for n=1:NumbDevs
        plot(r1[1,n,1:stop2[1]])
    end

    xlabel("Iteration", fontsize=label_fontsize)
    ylabel("Primal Residual (\$r_1\$)", fontsize=label_fontsize)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    grid(true)
    savefig(string(folder,"Convergence_Residual1.pdf"))

    clf()
    figure(4,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize =legend_fontsize-1)
    for s=1:Numb_Services
        for n=1:NumbDevs
            plot(r2[1,s,n,1:stop2[1]])
        end
    end
    # ylim(-0.00001,0.00001)
    xlabel("Iteration", fontsize=label_fontsize)
    ylabel("Primal Residual (\$r_2\$)", fontsize=label_fontsize)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    grid(true)
    savefig(string(folder,"Convergence_Residual2.pdf"))

    clf()
    figure(5,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize =legend_fontsize-1)
    for s=1:Numb_Services
        if(s==1)
            plot(rs_eta[s]*ones(t),label="Solver",linestyle=":", color="k", marker=markers[1], markersize=marker_size,markevery=Stride)
            plot(eta1[1,s,1:3], color="r",label="Centralized",linestyle="-", marker=markers[3], markersize=marker_size,markevery=Stride1)

        else
            plot(rs_eta[s]*ones(t),linestyle=":", color="k", marker=markers[1], markersize=marker_size,markevery=Stride)
            plot(eta1[1,s,1:3], color="r",linestyle="-", marker=markers[3], markersize=marker_size,markevery=Stride1)
        end
        plot(eta2[1,s,1:stop2[1]], color=colors1[s],label=string("Decentralized:S",s),linestyle="--", marker=markers[2], markersize=marker_size,markevery=Stride)
    end
    xlabel("Iteration", fontsize=label_fontsize)
    ylabel("Learning paramter (\$\\eta\$)", fontsize=label_fontsize)
    legend(loc="best",fontsize=legend_fontsize-1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    grid(true)
    savefig(string(folder,"Convergence_eta.pdf"))

    clf()
    figure(30,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize =legend_fontsize-1)
    for n=1:NumbDevs
        if(n==9)
            plot(rs_w[n]*ones(t),label="Solver",linestyle=":", color="k")
            plot(w1[1,n,1:3], color="r",label="Centralized",linestyle="-", marker=markers[3], markersize=marker_size,markevery=Stride1)
            plot(w2[1,n,1:stop2[1]],label=string("Decentralized:UE",n),linestyle="-")
        else
            plot(rs_w[n]*ones(t),linestyle=":", color="k")
            plot(w1[1,n,1:3],color="r",linestyle="-", marker=markers[3], markersize=marker_size,markevery=Stride1)
            plot(w2[1,n,1:stop2[1]], linestyle="-")
            # for s=1:Numb_Services
            #     plot(ws2[1,s,n,1:stop2[1]],label="miADMM",linestyle="--")
            # end
        end
        # plot(w1[1,n,1:stop1[1]],label="BCD",linestyle="--")
    end
    legend(loc="best",fontsize=legend_fontsize-1)
    xlabel("Iteration", fontsize=label_fontsize)
    ylabel("Fraction of Bandwidth (\$w\$)", fontsize=label_fontsize)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    grid(true)
    savefig(string(folder,"Convergence_w.pdf"))

    for n=1:NumbDevs
        if(n==1)
            clf()
            figure(7,figsize=fig_size)
            ax = subplot(1,1,1)
            ax.tick_params("both",labelsize =legend_fontsize-1)
            for s=1:Numb_Services
                if(s==1)
                    plot(rs_f[s,n]*ones(t),label="Solver", color="k",linestyle=":", marker=markers[1], markersize=marker_size,markevery=Stride)
                    plot(f1[1,s,n,1:3], color="r",label=string("Centralized"),linestyle="-", marker=markers[3], markersize=marker_size,markevery=Stride1)
                else
                    plot(rs_f[s,n]*ones(t),linestyle=":", color="k", marker=markers[1], markersize=marker_size,markevery=Stride)
                    plot(f1[1,s,n,1:3], color="r",linestyle="-", marker=markers[3], markersize=marker_size,markevery=Stride1)
                end

                plot(f2[1,s,n,1:stop2[1]], color=colors1[s],label=string("Decentralized:S",s),linestyle="--", marker=markers[2], markersize=marker_size,markevery=Stride)
                # plot(eta1[1,s,1:stop1[1]],label="BCD",linestyle="--")
            end
            xlabel("Iteration", fontsize=label_fontsize)
            ylabel("CPU frequency (\$f\$)", fontsize=label_fontsize)
            legend(loc="best",fontsize=legend_fontsize-1)
            tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
            grid(true)
            savefig(string(folder,"Convergence_f.pdf"))
        end
    end
end

# FOR BCD Only
function plot_convergence1(Obj1, Obj2, rs_Obj, r1, r2, eta1, eta2, rs_eta, w1, w2, ws2, rs_w, f1, f2, rs_f, stop1, stop2)
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
            plot(rs_eta[s]*ones(t),label="Solver",linestyle=":")
        else
            plot(rs_eta[s]*ones(t),linestyle=":")
        end

        plot(eta1[1,s,1:stop1[1]],label=string("BCD-Service",s),linestyle="-")
    end
    xlabel("Iteration",fontsize=label_fontsize1)
    ylabel("\$\\eta\$",fontsize=label_fontsize1)
    legend(loc="best",fontsize=legend_fontsize-1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Convergence_eta.pdf"))

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

function plot_comparison(rs_Obj, Obj_E, Obj_T, Heuristic_Obj,Heuristic_Obj_E, Heuristic_Obj_T, Heuristic_Obj1,Heuristic_Obj_E1, Heuristic_Obj_T1,
    rs_eta, Heuristic_eta, Heuristic_eta1)
    # # println("HERE")
    # # fig_size = (8,4)
    # Theta = zeros(Numb_Services)
    # for s =1:Numb_Services
    #     eta = -B2s[s]/(2*A2s[s])
    #     Theta[s] = A2s[s]*eta^2 + B2s[s]*eta
    #     # Theta[s] = A2s[s]*rs_eta[s]^2 + B2s[s]*rs_eta[s]
    #     println("bound1:", B2s[s]^2 +4*A2s[s])
    #     # println("bound2:", (-B2s[s]+sqrt(B2s[s]^2 +4*A2s[s]) )/(2*A2s[s]) )
    # end
    #
    # println("Big_Theta:",Theta)

    K_g = zeros(Numb_Services)
    for s=1:Numb_Services
        # K_g[s] = A1s[s]/(A2s[s]*rs_eta[s]^2 + B2s[s]*rs_eta[s])
        K_g[s] = 2*rho0*As[s]*(Bs[s]*Heuristic_eta[s]^2 +1)/(Cs[s]*Heuristic_eta[s] - Ds[s]*Heuristic_eta[s]^2)
    end

    clf()
    figure(20,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize =legend_fontsize-1)
    alg_labels = ["Heuristic 1","Heuristic 2","Optimal"]
    Service_labels = ["Service 1", "Service 2", "Service 3"]
    idx = [1:1:Numb_Services;]
    width = 0.25
    barlist1 = bar(idx.-width, K_g.*Heuristic_Obj_E, width, color=colors1, alpha=.5, hatch=patterns[1], label=alg_labels[1])
    barlist2 = bar(idx, K_g.*Heuristic_Obj_E1, width, color=colors1, alpha=.7, hatch=patterns[2], label=alg_labels[2])
    barlist3 = bar(idx.+width, K_g.*Obj_E, width, color=colors1, alpha=.9, hatch=patterns[3], label=alg_labels[3])

    xticks(idx, Service_labels, fontsize=label_fontsize1)
    ylabel("Energy Consumption of UEs", fontsize=label_fontsize1)
    legend(loc="best",fontsize=legend_fontsize-1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"E_Comparison.pdf"))


    K_g = zeros(Numb_Services)
    for s=1:Numb_Services
        # K_g[s] = A1s[s]/(A2s[s]*rs_eta[s]^2 + B2s[s]*rs_eta[s])
        K_g[s] = 2*rho0*As[s]*(Bs[s]*Heuristic_eta1[s]^2 +1)/(Cs[s]*Heuristic_eta1[s] - Ds[s]*Heuristic_eta1[s]^2)
    end

    clf()
    figure(21,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize =legend_fontsize-1)
    # alg_labels = ["Heuristic","Optimal"]
    # Service_labels = ["Service 1", "Service 2", "Service 3"]
    # idx = [1:1:Numb_Services;]
    # width = 0.25
    barlist1 = bar(idx.-width, K_g.*Heuristic_Obj_T, width, color=colors1, alpha=.5, hatch=patterns[1], label=alg_labels[1])
    barlist2 = bar(idx, K_g.*Heuristic_Obj_T1, width, color=colors1, hatch=patterns[2], alpha=.7, label=alg_labels[2])
    barlist3 = bar(idx.+width, K_g.*Obj_T, width, color=colors1, alpha=.9, hatch=patterns[3], label=alg_labels[3])

    ylim(0,3900)
    xticks(idx, Service_labels, fontsize=label_fontsize1)
    ylabel("Total Time", fontsize=label_fontsize1)
    legend(loc="best",fontsize=legend_fontsize-1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"T_Comparison.pdf"))


    clf()
    figure(22,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize =legend_fontsize-1)
    # alg_labels = ["Heuristic","Optimal"]
    # Service_labels = ["Service 1", "Service 2", "Service 3"]
    # idx = [1:1:Numb_Services;]
    # width = 0.25
    K_g = zeros(Numb_Services)
    for s=1:Numb_Services
        # K_g[s] = A1s[s]/(A2s[s]*rs_eta[s]^2 + B2s[s]*rs_eta[s])
        K_g[s] = 2*rho0*As[s]*(Bs[s]*rs_eta[s]^2 +1)/(Cs[s]*rs_eta[s] - Ds[s]*rs_eta[s]^2)
    end

    barlist1 = bar(idx.-width, K_g.*(Heuristic_Obj_E + kaps[1]*Heuristic_Obj_T), width, color=colors1, alpha=.5, hatch=patterns[1], label=alg_labels[1])
    barlist2 = bar(idx, K_g.*(Heuristic_Obj_E1 + kaps[1]*Heuristic_Obj_T1), width, color=colors1, hatch=patterns[2], alpha=.7, label=alg_labels[2])
    barlist3 = bar(idx.+width, K_g.*(Obj_E + kaps[1]*Obj_T), width, color=colors1, alpha=.9, hatch=patterns[3], label=alg_labels[3])

    xticks(idx, Service_labels, fontsize=label_fontsize1)
    ylabel("Total cost", fontsize=label_fontsize1)
    legend(loc="best",fontsize=legend_fontsize-1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Service_Cost_Comparison.pdf"))
end

function plot_data_distance(D_n,dist_list)
    clf()
    fig_size2 = (5.5,4.9)
    figure(23,figsize=fig_size2)
    ax = subplot(4,1,1)
    ax.tick_params("both",labelsize =legend_fontsize-3)
    alg_labels = ["Heuristic","Optimal"]
    # UE_labels = ["UE1", "UE2", "UE3", "UE4", "UE5", "UE6", "UE7", "UE8", "UE9", "UE10"]
    UE_labels = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    idx = [1:1:NumbDevs;]
    width = 0.35
    barlist1 = bar(idx, dist_list, width, color="k", alpha=.8, hatch=patterns[1])

    # xticks(idx, UE_labels, fontsize=label_fontsize1-4)
    ylabel("Distance", fontsize=label_fontsize1-2)
    # xlabel("UE", fontsize=label_fontsize1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.25)
    # savefig(string(folder,"Distance.pdf"))

    ax = subplot(4,1,2)
    # ax.tick_params("both",labelsize =legend_fontsize-1)

    barlist1 = bar(idx, D_n[1,:]*1e-6/8, width, color=colors1[1], alpha=.8, hatch=patterns[1])

    # xticks(idx, UE_labels, fontsize=label_fontsize1-4)
    ylabel("\$D_{1,n}\$", fontsize=label_fontsize1-2)
    # xlabel("UE", fontsize=label_fontsize1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.25)

    ax = subplot(4,1,3)
    # ax.tick_params("both",labelsize =legend_fontsize-1)

    barlist1 = bar(idx, D_n[2,:]*1e-6/8, width, color=colors1[2], alpha=.8, hatch=patterns[1])

    # xticks(idx, UE_labels, fontsize=label_fontsize1-4)
    ylabel("\$D_{2,n}\$", fontsize=label_fontsize1-2)
    # xlabel("UE", fontsize=label_fontsize1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.25)

    ax = subplot(4,1,4)
    # ax.tick_params("both",labelsize =legend_fontsize-1)

    barlist1 = bar(idx, D_n[3,:]*1e-6/8, width, color=colors1[3], alpha=.8, hatch=patterns[1])

    # xticks(idx, UE_labels, fontsize=label_fontsize1-4)
    ylabel("\$D_{3,n}\$", fontsize=label_fontsize1-2)
    xlabel("UE index", fontsize=label_fontsize1-2)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

    savefig(string(folder,"Distance_DataSize.pdf"))

end

function plot_pareto(rs_T_cmp1, rs_E_cmp1, rs_T_com1, rs_E_com1, rs_eta1)
    clf()
    figure(30,figsize=fig_size) #(8.7,5.8)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize =legend_fontsize)

    K_g = zeros(Numb_kaps1,Numb_Services)
    for s=1:Numb_Services
        for k= 1:Numb_kaps1
        # K_g[s] = A1s[s]/(A2s[s]*rs_eta[s]^2 + B2s[s]*rs_eta[s])
            K_g[k,s] = 2*rho0*As[s]*(Bs[s]*rs_eta1[k,s]^2 +1)/(Cs[s]*rs_eta1[k,s] - Ds[s]*rs_eta1[k,s]^2)
        end
    end


    Obj_E, Obj_T=zeros(Numb_kaps1), zeros(Numb_kaps1)

    for k=1:Numb_kaps1
        for s=1:Numb_Services
            Obj_E[k] += K_g[k,s]*(rs_E_com1[k,s] +K_l[s]*rs_E_cmp1[k,s])
            Obj_T[k] += K_g[k,s]*(rs_T_com1[k,s] +K_l[s]*rs_T_cmp1[k,s])
            # println("Obj_E1:",K_g[k,s]*(rs_E_com1[k,s] +K_l[s]*rs_E_cmp1[k,s]))
            # println("Obj_T1:",K_g[k,s]*(rs_T_com1[k,s] +K_l[s]*rs_T_cmp1[k,s]))
        end
    end
    println("Obj_E:",Obj_E)
    println("Obj_T:",Obj_T)

    plot(Obj_E, Obj_T, linestyle=":", marker=markers[2], markersize=marker_size)


    # Obj_E, Obj_T=zeros(Numb_kaps1, Numb_Services), zeros(Numb_kaps1, Numb_Services)
    # for s =1:Numb_Services
    #     Obj_E[:,s] = rs_E_com1[:,s] +K_l[s]*rs_E_cmp1[:,s]
    #     Obj_T[:,s] = rs_T_com1[:,s] +K_l[s]*rs_T_cmp1[:,s]
    # end
    #
    # for s=1:Numb_Services
    # # for s=1:1
    #     plot(K_g[:,s].*Obj_E[:,s],K_g[:,s].*Obj_T[:,s],linestyle=":", color=colors1[s], marker=markers[s], markersize=marker_size)
    # end

    ylabel("Total Time", fontsize=label_fontsize1)
    xlabel("Energy Consumption of UEs", fontsize=label_fontsize1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.25)
    grid(true)
    savefig(string(folder,"Pareto.pdf"))
end

function plot_priority(rs_T_cmp2, rs_E_cmp2, rs_T_com2, rs_E_com2, rs_eta2)
    clf()
    figure(31,figsize=fig_size) #(8.7,5.8)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize =legend_fontsize)
    Service_labels = ["Service 1", "Service 2", "Service 3"]

    K_g = zeros(Numb_Kap_Vals,Numb_Services)
    Obj_T=zeros(Numb_Kap_Vals,Numb_Services)
    Obj_E=zeros(Numb_Kap_Vals,Numb_Services)
    kap_p = zeros(Numb_Kap_Vals)

    for s=1:Numb_Services
        for k=1:Numb_Kap_Vals
            kap_p[k] = k
            K_g[k,s] = 2*rho0*As[s]*(Bs[s]*rs_eta2[k,s]^2 +1)/(Cs[s]*rs_eta2[k,s] - Ds[s]*rs_eta2[k,s]^2)
            Obj_T[k,s] = K_g[k,s]*(rs_T_com2[k,s] + K_l[s]*rs_T_cmp2[k,s])
            Obj_E[k,s] = K_g[k,s]*(rs_E_com2[k,s] +K_l[s]*rs_E_cmp2[k,s])
        end
        plot(kap_p,Obj_T[:,s], linestyle=":", color=colors1[s], marker=markers[s], markersize=marker_size, label = Service_labels[s])
    end

    legend(loc="best",fontsize=legend_fontsize)
    ylabel("Total Time", fontsize=label_fontsize1)
    xlabel("\$\\kappa_3\$", fontsize=label_fontsize1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.25)
    grid(true)
    savefig(string(folder,"Priority_Time.pdf"))

    clf()
    figure(32,figsize=fig_size) #(8.7,5.8)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize =legend_fontsize)

    for s=1:Numb_Services
        plot(kap_p,Obj_E[:,s], linestyle=":", color=colors1[s], marker=markers[s], markersize=marker_size, label = Service_labels[s])
    end

    legend(loc=4,fontsize=legend_fontsize)
    ylabel("Energy Consumption of UEs", fontsize=label_fontsize1)
    xlabel("\$\\kappa_3\$", fontsize=label_fontsize1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.25)
    grid(true)
    savefig(string(folder,"Priority_Energy.pdf"))
end

function save_result(rs_eta, rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, rs_Obj, rs_w, rs_f, Obj1, Obj2, r1, r2,
eta1, eta2, w1, w2, ws2, f1, f2, stop1, stop2, Obj_E, Obj_T, Heuristic_Obj, Heuristic_Obj_E, Heuristic_Obj_T, Heuristic_eta,
Heuristic_Obj1, Heuristic_Obj_E1, Heuristic_Obj_T1, Heuristic_eta1, rs_T_cmp1, rs_E_cmp1, rs_T_com1, rs_E_com1, rs_eta1,
rs_T_cmp2, rs_E_cmp2, rs_T_com2, rs_E_com2, rs_eta2)
    h5open(string("result",NumbDevs,".h5"), "w") do file
        # write(file,"kaps", kaps)
        write(file,"rs_eta", rs_eta)
        write(file,"rs_T_cmp", rs_T_cmp)
        write(file,"rs_E_cmp", rs_E_cmp)
        write(file,"rs_T_com", rs_T_com)
        write(file,"rs_E_com", rs_E_com)
        write(file,"rs_Obj", rs_Obj)
        write(file,"rs_w", rs_w)
        write(file,"rs_f", rs_f)
        write(file,"Obj1", Obj1)
        write(file,"Obj2", Obj2)
        write(file,"r1", r1)
        write(file,"r2", r2)
        write(file,"eta1", eta1)
        write(file,"eta2", eta2)
        write(file,"w1", w1)
        write(file,"w2", w2)
        write(file,"ws2", ws2)
        write(file,"f1", f1)
        write(file,"f2", f2)
        write(file,"stop1", stop1)
        write(file,"stop2", stop2)
        write(file,"Obj_E", Obj_E)
        write(file,"Obj_T", Obj_T)
        write(file,"Heuristic_Obj", Heuristic_Obj)
        write(file,"Heuristic_Obj_E", Heuristic_Obj_E)
        write(file,"Heuristic_Obj_T", Heuristic_Obj_T)
        write(file,"Heuristic_eta", Heuristic_eta)
        write(file,"Heuristic_Obj1", Heuristic_Obj1)
        write(file,"Heuristic_Obj_E1", Heuristic_Obj_E1)
        write(file,"Heuristic_Obj_T1", Heuristic_Obj_T1)
        write(file,"Heuristic_eta1", Heuristic_eta1)
        write(file,"rs_T_cmp1", rs_T_cmp1)
        write(file,"rs_E_cmp1", rs_E_cmp1)
        write(file,"rs_T_com1", rs_T_com1)
        write(file,"rs_E_com1", rs_E_com1)
        write(file,"rs_eta1", rs_eta1)
        write(file,"rs_T_cmp2", rs_T_cmp2)
        write(file,"rs_E_cmp2", rs_E_cmp2)
        write(file,"rs_T_com2", rs_T_com2)
        write(file,"rs_E_com2", rs_E_com2)
        write(file,"rs_eta2", rs_eta2)
    end
end

function read_result(filename)
    h5open(filename, "r") do file
        # kaps =read(file,"kaps")
        rs_eta = read(file,"rs_eta")
        rs_T_cmp  = read(file,"rs_T_cmp")
        rs_E_cmp = read(file,"rs_E_cmp")
        rs_T_com = read(file,"rs_T_com")
        rs_E_com = read(file,"rs_E_com")
        rs_Obj = read(file,"rs_Obj")
        rs_w = read(file,"rs_w")
        rs_f = read(file,"rs_f")
        Obj1 = read(file,"Obj1")
        Obj2 = read(file,"Obj2")
        r1 = read(file,"r1")
        r2 = read(file,"r2")
        eta1 = read(file,"eta1")
        eta2 = read(file,"eta2")
        w1 = read(file,"w1")
        w2 = read(file,"w2")
        ws2 = read(file,"ws2")
        f1 = read(file,"f1")
        f2 = read(file,"f2")
        stop1 = read(file,"stop1")
        stop2 = read(file,"stop2")
        Obj_E = read(file,"Obj_E")
        Obj_T = read(file,"Obj_T")
        Heuristic_Obj = read(file,"Heuristic_Obj")
        Heuristic_Obj_E = read(file,"Heuristic_Obj_E")
        Heuristic_Obj_T = read(file,"Heuristic_Obj_T")
        Heuristic_eta = read(file,"Heuristic_eta")
        Heuristic_Obj1 = read(file,"Heuristic_Obj1")
        Heuristic_Obj_E1 = read(file,"Heuristic_Obj_E1")
        Heuristic_Obj_T1 = read(file,"Heuristic_Obj_T1")
        Heuristic_eta1 = read(file,"Heuristic_eta1")
        rs_T_cmp1  = read(file,"rs_T_cmp1")
        rs_E_cmp1 = read(file,"rs_E_cmp1")
        rs_T_com1 = read(file,"rs_T_com1")
        rs_E_com1 = read(file,"rs_E_com1")
        rs_eta1 = read(file,"rs_eta1")
        rs_T_cmp2  = read(file,"rs_T_cmp2")
        rs_E_cmp2 = read(file,"rs_E_cmp2")
        rs_T_com2 = read(file,"rs_T_com2")
        rs_E_com2 = read(file,"rs_E_com2")
        rs_eta2 = read(file,"rs_eta2")

        return rs_eta, rs_T_cmp, rs_E_cmp, rs_T_com, rs_E_com, rs_Obj, rs_w, rs_f, Obj1, Obj2, r1, r2,
        eta1, eta2, w1, w2, ws2, f1, f2, stop1, stop2, Obj_E, Obj_T, Heuristic_Obj, Heuristic_Obj_E, Heuristic_Obj_T, Heuristic_eta,
        Heuristic_Obj1, Heuristic_Obj_E1, Heuristic_Obj_T1, Heuristic_eta1, rs_T_cmp1, rs_E_cmp1, rs_T_com1, rs_E_com1, rs_eta1,
        rs_T_cmp2, rs_E_cmp2, rs_T_com2, rs_E_com2, rs_eta2
    end
end


function save_result_iteration(stop1,stop2,stop3,stop4,rs_Objs)
    h5open(string("result_iter_",NUM_SIM,".h5"), "w") do file
        # write(file,"kaps", kaps)
        write(file,"stop1", stop1)
        write(file,"stop2", stop2)
        write(file,"stop3", stop3)
        write(file,"stop4", stop4)
        write(file,"rs_Objs", rs_Objs)
    end
end

function read_result_iteration()
    h5open(string("result_iter_",NUM_SIM,".h5"), "r") do file
        # write(file,"kaps", kaps)
        stop1 = read(file,"stop1")
        stop2  = read(file,"stop2")
        stop3 = read(file,"stop3")
        stop4 = read(file,"stop4")
        rs_Objs = read(file,"rs_Objs")
    end
    return stop1, stop2, stop3, stop4, Obj2,Obj3,Obj4
end
