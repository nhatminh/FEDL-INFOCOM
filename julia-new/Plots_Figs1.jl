# using Plots
using PyPlot
using PyCall
@pyimport matplotlib.patches as patch

# fig_size = (7.,5.1)
fig_size = (6.,4.3)
fig_size1 = (5.5,4.)
label_fontsize = 18-1.5
legend_fontsize = label_fontsize - 2
patterns = ["","."]
line_style = ["-","--",":"]
label_fontsize1 = label_fontsize
marker_size=6
l_width=1.2
l_width1=0.5
# stride = convert(Int, max_iters/10) +1
stride = 6

colors=["m","b","coral","g","k","r"]
# algs = ["PBCD", "Consensus_BCD2", "JP-ADMM", "JP-ADMM_BCD4","IpOpt Solver","Exhaustive Search"]
markers = ["x","o",">","^", ".","s"]

folder = string("figs//")

function plot_sub1_T(T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, levels)
    clf()
    cfig = figure(1,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize)

    plot(levels,T_cmp1,color=colors[6],linestyle="-",linewidth=l_width+0.3)
    # plot(levels,T_cmp+ 0.02*maximum(T_cmp1),color="gold",linestyle=":",linewidth=l_width+0.3,label="Solver")
    # # plot(levels,T_cmp1+0.8,color=colors[6],linestyle="-",linewidth=l_width,label="\$T_{cmp}^*\$")
    # plot(levels,Tcmp_N1,color=colors[2],linestyle="--",linewidth=l_width+0.2,label="\$T_{\\mathcal{N}_1}\$")
    # plot(levels,Tcmp_N2,color=colors[4],linestyle="-.",linewidth=l_width+0.2,label="\$T_{\\mathcal{N}_2}\$")
    # plot(levels,Tcmp_N3,color=colors[5],linestyle="-",linewidth=l_width+0.2,label="\$T_{\\mathcal{N}_3}\$")

    # legend(loc="best",fontsize=legend_fontsize+2)
    xlabel("\$L_{cmp}\$",fontsize=label_fontsize1+2)
    xscale("log")
    ylabel("\$T^*_{cmp}\$ (sec)",fontsize=label_fontsize1+2)
    # ylim(-0.3, 1.06*maximum(T_cmp1))
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_T_rs.pdf"))
end

function plot_sub1_N(N1, N2, N3, levels)
    clf()
    cfig = figure(2,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize)
    # plot(kaps,T_cmp,color=colors[1],linestyle="-",linewidth=l_width,label="Solver")
    step(levels,N1,color=colors[4],linestyle="-",linewidth=l_width,label="\$\\mathcal{N}_1\$", where="post", marker=markers[2], markersize=marker_size, markevery=11)
    step(levels,N2,color=colors[3],linestyle="-",linewidth=l_width,label="\$\\mathcal{N}_2\$", where="post", marker=markers[3], markersize=marker_size, markevery=11)
    step(levels,N3,color=colors[2],linestyle="-",linewidth=l_width,label="\$\\mathcal{N}_3\$", where="post", marker=markers[6], markersize=marker_size, markevery=11)

    legend(loc="best",fontsize=legend_fontsize)
    xlabel("\$L_{cmp}\$",fontsize=label_fontsize1+2)
    xscale("log")
    ylabel("Three subsets by Alg.1",fontsize=label_fontsize1+2)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_N_rs.pdf"))
end

function plot_sub1_f(f1,levels)
    clf()
    cfig = figure(3,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize)
    plot(levels,f_min[1]*ones(Numb_D)*1e-9,linestyle=":",color=colors[6])

    for n = 1:NumbDevs
        plot(levels,f1[:,n],linestyle="-",linewidth=l_width1)
    end


    # legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$L_{cmp}\$",fontsize=label_fontsize1+2)
    # ylim(0,maximum(f1) + 0.2 )
    ylabel("f (GHz)",fontsize=label_fontsize1+2)
    xscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_f_rs.pdf"))
end

function plot_sub2_tau(tau1,levels)
    clf()
    cfig = figure(4,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize)

    # for n = 1:5
    #     plot(tau_ratios,tau1[:,n], color=colors[n], linestyle="-",linewidth=l_width,label=string("UE ",n))
    # end
    for n = 1:NumbDevs
        plot(levels,tau1[:,n],linestyle="-",linewidth=l_width1)
    end

    max_tau = maximum(tau1[1,:])

    # legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$L_{com}\$",fontsize=label_fontsize1+2)
    ylabel("\$\\tau_n\$ (sec)",fontsize=label_fontsize1+2)
    xscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub2_Tau_rs.pdf"))
end

function plot_sub2_p(p1,levels)
    clf()
    cfig = figure(5,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize)

    plot(levels,Ptx_Max*ones(Numb_Dis),linestyle=":",color=colors[6])
    plot(levels,Ptx_Min*ones(Numb_Dis),linestyle=":",color=colors[6])
    # for n = 1:5
    #     plot(tau_ratios,p1[:,n],color=colors[n],linestyle="-",linewidth=l_width,label=string("UE ",n))
    # end
    for n = 1:NumbDevs
        plot(levels,p1[:,n],linestyle="-",linewidth=l_width1)
    end

    # legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$L_{com}\$",fontsize=label_fontsize1+2)
    ylabel("p (Watt)",fontsize=label_fontsize1+2)
    xscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub2_p_rs.pdf"))
end

function plot_sub2_Tcom(T_com1,levels)
    clf()
    cfig = figure(5,figsize=fig_size)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize)

    plot(levels,T_com1,linestyle="-",color=colors[6])
    # for n = 1:5
    #     plot(tau_ratios,p1[:,n],color=colors[n],linestyle="-",linewidth=l_width,label=string("UE ",n))
    # end
    # for n = 1:NumbDevs
    #     plot(levels,p1[:,n],linestyle="-",linewidth=l_width1)
    # end

    # legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$L_{com}\$",fontsize=label_fontsize1+2)
    ylabel("\$T^*_{com}\$ (sec)",fontsize=label_fontsize1+2)
    xscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub2_Tcom_rs.pdf"))
end

function plot_sub3_cvx(Theta1, Obj1, T_cmp1, E_cmp1, T_com1, E_com1,levels)
    Numb_Levels = size(levels)[1]

    clf()
    figure(6,figsize=fig_size)
    x = collect(1.e-5:0.001:0.99)
    obj   = zeros(size(x)[1])
    glob_cost_iter = zeros(size(x)[1])
    glob_numb_iter = zeros(size(x)[1])
    id = 40
    println("Convex for kappa: ",  kaps[id])
    for i=1:size(x)[1]
        obj[i] = 1/(1 - x[i])* (E_com1[id] - log(x[i])*E_cmp1[id] + kaps[id] * (T_com1[id] - log(x[i])*T_cmp1[id]))
        glob_cost_iter[i] = E_com1[id] - log(x[i])*E_cmp1[id] + kaps[id] * (T_com1[id] - log(x[i])*T_cmp1[id])
        glob_numb_iter[i] = 1/(1 - x[i])
        # obj[i]   = obj_E[i] + obj_T[i]
    end
    plot(x, obj,linestyle="-",color="k", label=string("Objective: \$\\kappa\$ =", kaps[id]))
    plot(x, glob_cost_iter,linestyle="--",color=colors[2], label=string("\$E_{glob} + \\kappa * T_{glob}\$"))
    plot(x, glob_numb_iter,linestyle="--",color=colors[3], label=string("\$ 1/(1 - \\theta)\$"))
    # println(x)
    plot(Theta1[id], Obj1[id],color="r", marker=markers[2], markersize=marker_size)

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\theta\$",fontsize=label_fontsize1+1)
    # ylabel("Objective",fontsize=label_fontsize1+1)
    yscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub3_obj.pdf"))
    println("Theta: ", minimum(Theta1), " - ", maximum(Theta1))
end

function plot_ratios(T_cmp1, E_cmp1, T_com1, E_com1, T_cmp12, E_cmp12, T_com12, E_com12, levels, levels2)
    clf()
    cfig = figure(11,figsize=fig_size1)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize)

    plot(kaps,T_cmp1[1,:]./(T_cmp1[1,:]+T_com1[1,:]),linestyle=line_style[1],marker=markers[2], markersize=marker_size,markevery=5, label=string("\$L_{cmp}\$=",round(levels[1],digits=2)))
    plot(kaps,T_cmp1[end,:]./(T_cmp1[end,:]+T_com1[end,:]),linestyle=line_style[1],marker=markers[2], markersize=marker_size,markevery=5, label=string("\$L_{cmp}\$=",round(levels[end],digits=2)))

    plot(kaps,T_cmp12[1,:]./(T_cmp12[1,:]+T_com12[1,:]),linestyle=line_style[2],marker=markers[3], markersize=marker_size,markevery=5, label=string("\$L_{com}\$=",round(levels2[1],digits=2)))
    plot(kaps,T_cmp12[end,:]./(T_cmp12[end,:]+T_com12[end,:]),linestyle=line_style[2],marker=markers[3], markersize=marker_size,markevery=5, label=string("\$L_{com}\$=",round(levels2[end],digits=2)))

    legend(loc=1,fontsize=legend_fontsize-1)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+2)
    ylabel("\$T^*_{cmp}/(T^*_{cmp}+ T^*_{com})\$",fontsize=label_fontsize1+1)
    ylim(0, 1.05)
    xscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"sub3_T_ratio_rs.pdf"))

    clf()
    cfig = figure(11,figsize=fig_size1)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize)

    plot(kaps,E_cmp1[1,:]./(E_cmp1[1,:]+E_com1[1,:]),linestyle=line_style[1],marker=markers[2], markersize=marker_size,markevery=5, label=string("\$L_{cmp}\$=",round(levels[1],digits=2)))
    plot(kaps,E_cmp1[end,:]./(E_cmp1[end,:]+E_com1[end,:]),linestyle=line_style[1],marker=markers[2], markersize=marker_size,markevery=5, label=string("\$L_{cmp}\$=",round(levels[end],digits=2)))

    plot(kaps,E_cmp12[1,:]./(E_cmp12[1,:]+E_com12[1,:]),linestyle=line_style[2],marker=markers[3], markersize=marker_size,markevery=5, label=string("\$L_{com}\$=",round(levels2[1],digits=2)))
    plot(kaps,E_cmp12[end,:]./(E_cmp12[end,:]+E_com12[end,:]),linestyle=line_style[2],marker=markers[3], markersize=marker_size,markevery=5, label=string("\$L_{com}\$=",round(levels2[end],digits=2)))

    legend(loc="best",fontsize=legend_fontsize-1)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+2)
    ylabel("\$E^*_{cmp}/(E^*_{cmp}+ E^*_{com})\$",fontsize=label_fontsize1+1)
    ylim(0,1.05)
    xscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"sub3_E_ratio_rs.pdf"))
end

# function plot_sub3_kappa_theta(Theta, d_eta, levels, sub)
#     if(sub==1)
#         round_numb = 2
#         lbl_lv = "\$L_{cmp}\$"
#     else
#         round_numb = 2
#         lbl_lv = "\$L_{com}\$"
#     end
#
#     Numb_Levels = size(levels)[1]
#
#     clf()
#     cfig = figure(10,figsize=fig_size1)
#     ax = subplot(1,1,1)
#     ax.tick_params("both",labelsize=legend_fontsize)
#
#     # plot(Numb_devs, Objs_E[:,11],linestyle="--",color=colors[1],marker=markers[1], markersize=marker_size, label=string("\$\\kappa\$ =", kaps[11]))
#     plot(kaps, 1 ./d_eta[1,:],linestyle=line_style[1],marker=markers[2], markersize=marker_size,markevery=5,label=string("\$\\eta\$, ",lbl_lv,"=",round(levels[1],digits=round_numb)))
#     plot(kaps, Theta[1,:],linestyle=line_style[1],marker=markers[3], markersize=marker_size,markevery=5,label=string("\$\\theta\$, ",lbl_lv,"=",round(levels[1],digits=round_numb)))
#     plot(kaps, 1 ./d_eta[end,:],linestyle=line_style[2],marker=markers[2], markersize=marker_size,markevery=5,label=string("\$\\eta\$, ",lbl_lv,"=",round(levels[end],digits=round_numb)))
#     plot(kaps, Theta[end,:],linestyle=line_style[2],marker=markers[3], markersize=marker_size,markevery=5,label=string("\$\\theta\$, ",lbl_lv,"=",round(levels[end],digits=round_numb)))
#
#     legend(loc="best",fontsize=legend_fontsize-1)
#     xlabel("\$\\kappa\$",fontsize=label_fontsize1+2)
#     ylabel("\$\\theta^*\$ and \$\\eta\$",fontsize=label_fontsize1+1)
#     xscale("log")
#     tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
#     savefig(string(folder,"sub3_kappa_theta_rs",sub,".pdf"))
# end

function plot_sub3_kappa_theta_eta(Theta, theta, eta, levels, sub)
    if(sub==1)
        round_numb = 2
        lbl_lv = "\$L_{cmp}\$"
    else
        round_numb = 2
        lbl_lv = "\$L_{com}\$"
    end

    Numb_Levels = size(levels)[1]

    clf()
    cfig = figure(10,figsize=fig_size1)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize)

    # plot(kaps, Theta,linestyle="--",color=colors[1],label="\$\\Theta^*\$")
    # plot(kaps, theta,linestyle="--",color=colors[2],label="\$\\theta^*\$")
    # plot(kaps, eta,linestyle="--",color=colors[3],label="\$\\eta^*\$")

    plot(kaps, eta[1,:],linestyle=line_style[1],marker=markers[2], markersize=marker_size,markevery=5,label=string("\$\\eta^*\$, ",lbl_lv,"=",round(levels[1],digits=round_numb)))
    plot(kaps, theta[1,:],linestyle=line_style[1],marker=markers[3], markersize=marker_size,markevery=5,label=string("\$\\theta^*\$, ",lbl_lv,"=",round(levels[1],digits=round_numb)))
    plot(kaps, eta[end,:],linestyle=line_style[2],marker=markers[2], markersize=marker_size,markevery=5,label=string("\$\\eta^*\$, ",lbl_lv,"=",round(levels[end],digits=round_numb)))
    plot(kaps, theta[end,:],linestyle=line_style[2],marker=markers[3], markersize=marker_size,markevery=5,label=string("\$\\theta^*\$, ",lbl_lv,"=",round(levels[end],digits=round_numb)))

    legend(loc="best",fontsize=legend_fontsize-1)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+2)
    ylabel("\$\\theta^*\$ and \$\\eta\$",fontsize=label_fontsize1+1)
    xscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"sub3_kappa_theta_rs",sub,".pdf"))
end

function plot_sub3_equation(d_eta,levels)
    Numb_Levels = size(levels)[1]

    clf()
    x = collect(1.e-6:0.001:0.999)

    figure(7,figsize=fig_size)
    plot(x, 1 ./x + log.(x),linestyle="-",color=colors[6])
    for k = 1:Numb_kaps
        plot(x,d_eta[k]*ones(size(x)),linestyle=":",color="k")
    end

    legend(loc="best",fontsize=legend_fontsize-1)
    xlim(0, 0.5)
    ylim(0.98,maximum(d_eta)+0.1*maximum(d_eta))
    xlabel("\$\\theta\$",fontsize=label_fontsize1+1)
    ylabel("\$\\log(e^{1/\\theta} \\theta)\$",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub3_eq.pdf"))
end

function plot_numerical_pareto(Theta1, theta, T_cmp1, E_cmp1, T_com1, E_com1, levels, sub)
    # println("Theta1",Theta1)
    # println("T_cmp1",T_cmp1)
    # println("E_cmp1",E_cmp1)
    # println("T_com1",T_com1)
    # println("E_com1",E_com1)
    # println("levels",levels)

    if(sub==1)
        round_numb = 2
        lbl_lv = "\$L_{cmp}\$"
    else
        round_numb = 2
        lbl_lv = "\$L_{com}\$"
    end

    Numb_Levels = size(levels)[1]
    # if sub == 1
    #     # idx_levels = [1, 9]
    #     idx_levels = 1:3
    # else
    #     # idx_levels = [1, 9]
    #     idx_levels = 1:3
    # end
    idx_levels = [1,3,5]
    clf()
    cfig = figure(30,figsize=fig_size1)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize)

    x = 1

    E_obj   = zeros(Numb_kaps,3)
    T_obj   = zeros(Numb_kaps,3)

    for idx in idx_levels
        for i=1:Numb_kaps
            E_obj[i,x] = 1/(Theta1[idx,i])* (E_com1[idx,i] + 1/gamma*(log(C) - log(theta[idx,i]))*E_cmp1[idx,i])
            T_obj[i,x] = 1/(Theta1[idx,i])* (T_com1[idx,i] + 1/gamma*(log(C) - log(theta[idx,i]))*T_cmp1[idx,i])
        end

        x+=1
    end
    # plot(E_obj[:,1], T_obj[:,1],linestyle=line_style[1])
    # plot(E_obj[11,1], T_obj[11,1],linestyle="--",color="k", marker=markers[2], markersize=marker_size)
    # plot(E_obj[13,1], T_obj[13,1],linestyle="--",color="k", marker=markers[2], markersize=marker_size)
    # plot(E_obj[22,1], T_obj[22,1],linestyle="--",color="k", marker=markers[2], markersize=marker_size)
    #
    # annotate(string("\$\\kappa\$=",kaps[11]), xy=[E_obj[11,1]+0.5;T_obj[11,1]], xycoords="data",size=19)
    # annotate(string("\$\\kappa\$=",kaps[13]), xy=[E_obj[13,1]+0.5;T_obj[13,1]], xycoords="data",size=19)
    # annotate(string("\$\\kappa\$=",kaps[22]), xy=[E_obj[22,1]-2.5;T_obj[22,1]+0.5], xycoords="data",size=19)

    plot(E_obj[:,1], T_obj[:,1],linestyle=line_style[1],label=string(lbl_lv,"=",round(levels[1],digits=round_numb)))
    plot(E_obj[:,2], T_obj[:,2],linestyle="-.",label=string(lbl_lv,"=",round(levels[3],digits=round_numb)))
    plot(E_obj[:,3], T_obj[:,3],linestyle=line_style[3],label=string(lbl_lv,"=",round(levels[5],digits=round_numb)))

    plot(E_obj[11,:], T_obj[11,:],linestyle="--",color="k", marker=markers[2], markersize=marker_size)
    plot(E_obj[13,:], T_obj[13,:],linestyle="--",color="k", marker=markers[2], markersize=marker_size)
    plot(E_obj[22,:], T_obj[22,:],linestyle="--",color="k", marker=markers[2], markersize=marker_size)

    annotate(string("\$\\kappa\$=",kaps[11]), xy=[E_obj[11,3]+0.5;T_obj[11,3]], xycoords="data",size=19)
    annotate(string("\$\\kappa\$=",kaps[13]), xy=[E_obj[13,3]+0.5;T_obj[13,3]], xycoords="data",size=19)
    annotate(string("\$\\kappa\$=",kaps[22]), xy=[E_obj[22,3]-2.5;T_obj[22,3]+0.5], xycoords="data",size=19)

    legend(loc="best",fontsize=legend_fontsize)
    xlabel("Energy Cost",fontsize=label_fontsize1+1)
    ylabel("Time Cost",fontsize=label_fontsize1+1)

    if(sub==1)
        xlim(0,1250)
        ylim(50,700)
    else
        xlim(0,1250)
        ylim(50,550)
    end

    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"pareto_rs",sub,".pdf"))
end

function plot_total_cost(Obj, levels, sub)
    if(sub==1)
        round_numb = 2
        lbl_lv = "\$L_{cmp}\$"
    else
        round_numb = 2
        lbl_lv = "\$L_{com}\$"
    end

    Numb_Levels = size(levels)[1]
    clf()
    cfig = figure(32,figsize=fig_size1)
    ax = subplot(1,1,1)
    ax.tick_params("both",labelsize=legend_fontsize)

    plot(kaps, Obj[1,:],linestyle=line_style[1],label=string(lbl_lv,"=",round(levels[1],digits=round_numb)))
    plot(kaps, Obj[3,:],linestyle=line_style[2],label=string(lbl_lv,"=",round(levels[3],digits=round_numb)))
    plot(kaps, Obj[end,:],linestyle=line_style[3],label=string(lbl_lv,"=",round(levels[end],digits=round_numb)))

    legend(loc="best",fontsize=legend_fontsize)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+2)
    ylabel("FEDL's obj",fontsize=label_fontsize1+1)
    xscale("log")
    ylim(0,1000)
    xlim(1e-1,5e0)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"total1_rs",sub,".pdf"))

    println("kappa:", [kaps[11], kaps[13], kaps[22]])
    println("Cost scale up:", [round(Obj[3,11]/Obj[1,11],digits=3), round(Obj[end,11]/Obj[1,11],digits=3)])
    println("Cost scale up:", [round(Obj[3,13]/Obj[1,13],digits=3), round(Obj[end,13]/Obj[1,13],digits=3)])
    println("Cost scale up:", [round(Obj[3,22]/Obj[1,22],digits=3), round(Obj[end,22]/Obj[1,22],digits=3)])

    # clf()
    # figure(10,figsize=fig_size)
    # plot(levels, Obj[:,1],linestyle="-",label=string("\$\\kappa\$=",kaps[1]))
    # plot(levels, Obj[:,10],linestyle="-",label=string("\$\\kappa\$=",kaps[10]))
    # plot(levels, Obj[:,end],linestyle="-",label=string("\$\\kappa\$=",kaps[end]))
    #
    # legend(loc="best",fontsize=legend_fontsize-2)
    # xlabel(string("level",sub),fontsize=label_fontsize1+1)
    # ylabel("Total Cost",fontsize=label_fontsize1+1)
    # # yscale("log")
    # tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    # savefig(string(folder,"total2_rs",sub,".pdf"))
end

function save_result(filename,Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta, levels1, levels2,
    Theta, theta, eta)
    h5open(filename, "w") do file
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
        write(file,"levels1", levels1)
        write(file,"levels2", levels2)
        write(file,"Theta", Theta)
        write(file,"theta", theta)
        write(file,"eta", eta)
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
        levels1 = read(file,"levels1")
        levels2 = read(file,"levels2")
        Theta = read(file,"Theta")
        theta = read(file,"theta")
        eta = read(file,"eta")
        return Theta1, Obj1, Obj_E, Obj_T, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta, levels1, levels2,
            Theta, theta, eta
    end
end
