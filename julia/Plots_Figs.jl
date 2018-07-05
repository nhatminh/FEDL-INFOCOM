using PyPlot

# fig_size = (7.,5.1)
fig_size = (6.,4.3)
label_fontsize = 18-1.5
legend_fontsize = label_fontsize - 4
patterns = ["","."]

label_fontsize1 = label_fontsize
marker_size=6
l_width=1.2
# stride = convert(Int, max_iters/10) +1
stride = 6

colors=["m","b","coral","g","k","r"]
algs = ["PBCD", "Consensus_BCD2", "JP-ADMM", "JP-ADMM_BCD4","IpOpt Solver","Exhaustive Search"]
markers = ["x","o",">","^", "."]

folder = string("figs//")

function plot_sub1_T(kaps, T_cmp, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3)
    clf()
    figure(1,figsize=fig_size)
    # plot(kaps,T_cmp,color=colors[1],linestyle="-",linewidth=l_width,label="Solver")
    plot(kaps,T_cmp1,color=colors[6],linestyle="-",linewidth=l_width,label="T_cmp_opt")
    plot(kaps,Tcmp_N1,color=colors[2],linestyle="-.",linewidth=l_width,label="Tcmp_N1")
    plot(kaps,Tcmp_N2,color=colors[4],linestyle="-.",linewidth=l_width,label="Tcmp_N2")
    plot(kaps,Tcmp_N3,color=colors[5],linestyle="-.",linewidth=l_width,label="Tcmp_N3")

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    xscale("log")
    ylabel("\$T_{cmp}\$ (sec)",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_T.png"))
end

function plot_sub1_N(kaps, N1, N2, N3)
    clf()
    figure(1,figsize=fig_size)
    # plot(kaps,T_cmp,color=colors[1],linestyle="-",linewidth=l_width,label="Solver")
    plot(kaps,N1,color=colors[2],linestyle="--",linewidth=l_width,label="N1")
    plot(kaps,N2,color=colors[4],linestyle="--",linewidth=l_width,label="N2")
    plot(kaps,N3,color=colors[5],linestyle="--",linewidth=l_width,label="N3")

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    xscale("log")
    ylabel("Numb of Devs",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_N.png"))
end

function plot_sub1_f(kaps, f, f1)
    clf()
    figure(2,figsize=fig_size)
    plot(kaps,f_min[1]*ones(size(kaps))*1e-9,linestyle=":",color=colors[6])

    if (HETEROGENEOUS == 0) # Homogeneous
        plot(kaps,f_max[1]*ones(size(kaps))*1e-9,linestyle="--",color=colors[6])
    end

    for n = 1:NumbDevs
        if (HETEROGENEOUS > 0)
            plot(kaps,f_max[n]*ones(size(kaps))*1e-9,linestyle="--",color=colors[n])
            # plot(kaps,f_min[n]*ones(size(kaps))*1e-9,linestyle=":",color=colors[n])
        end

        plot(kaps,f1[:,n],color=colors[n],linestyle="-",linewidth=l_width,label=string("Dev ",n))
    end
    # plot(kaps,T_comps1,marker=markers[2],color=colors[2], markersize=marker_size,linestyle=":",linewidth=l_width,label="Tcmp1")

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    xscale("log")
    ylabel("f (GHz)",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_f.png"))
end

function plot_sub2_tau(kaps, tau, tau1)
    clf()
    figure(3,figsize=fig_size)
    for n = 1:NumbDevs
        plot(kaps,tau1[:,n], color=colors[n], linestyle="-",linewidth=l_width,label=string("Dev ",n))
    end
    # plot(kaps,T_comps1,marker=markers[2],color=colors[2], markersize=marker_size,linestyle=":",linewidth=l_width,label="Tcmp1")

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    xscale("log")
    ylabel("\$\\tau_n\$ (sec)",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub2_Tau.png"))
end

function plot_sub2_p(kaps, p, p1)
    clf()
    figure(4,figsize=fig_size)
    plot(kaps,Ptx_Max*ones(size(kaps)),linestyle=":",color=colors[6])
    plot(kaps,Ptx_Min*ones(size(kaps)),linestyle=":",color=colors[6])
    for n = 1:NumbDevs
        plot(kaps,p1[:,n],color=colors[n],linestyle="-",linewidth=l_width,label=string("Dev ",n))
    end
    # plot(kaps,T_comps1,marker=markers[2],color=colors[2], markersize=marker_size,linestyle=":",linewidth=l_width,label="Tcmp1")

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    xscale("log")
    ylabel("p (Watt)",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub2_p.png"))
end

function plot_sub3_cvx(kaps, Theta, Theta1, Obj, Obj1, T_cmp1, E_cmp1, T_com1, E_com1)
    clf()
    figure(5,figsize=fig_size)
    # plot(Theta,Obj,color=colors[1],linestyle="-",linewidth=l_width,label="Solver")
    # plot(Theta1,Obj1,color=colors[2],linestyle=":",linewidth=l_width,label="Closed-Form")
    #
    # legend(loc="best",fontsize=legend_fontsize-2)
    # xlabel("\$\\Theta\$",fontsize=label_fontsize1+1)
    # ylabel("Objective",fontsize=label_fontsize1+1)
    # tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    # savefig(string(folder,"Sub3_obj.png"))

    x = collect(1.e-5:0.001:0.99)
    obj   = zeros(size(x)[1])
    glob_cost_iter = zeros(size(x)[1])
    glob_numb_iter = zeros(size(x)[1])
    id = 19
    println("Convex for kappa: ",  kaps[id])
    for i=1:size(x)[1]
        obj[i] = 1/(1 - x[i])* (E_com1[id] - log(x[i])*E_cmp1[id] + kaps[id] * (T_com1[id] - log(x[i])*T_cmp1[id]))
        glob_cost_iter[i] = E_com1[id] - log(x[i])*E_cmp1[id] + kaps[id] * (T_com1[id] - log(x[i])*T_cmp1[id])
        glob_numb_iter[i] = 1/(1 - x[i])
        # obj[i]   = obj_E[i] + obj_T[i]
    end
    plot(x, obj,linestyle="-",color="k", label=string("Objective: \$\\kappa\$ =", kaps[id]))
    plot(x, glob_cost_iter,linestyle="--",color=colors[2], label=string("\$E_{glob} + \\kappa * T_{glob}\$"))
    plot(x, glob_numb_iter,linestyle="--",color=colors[3], label=string("\$ 1/(1 - \\Theta)\$"))
    # println(x)
    plot(Theta1[id], Obj1[id],color="r", marker=markers[2], markersize=marker_size)

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\Theta\$",fontsize=label_fontsize1+1)
    # ylabel("Objective",fontsize=label_fontsize1+1)
    yscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub3_obj.png"))
    println("Theta: ", minimum(Theta1), " - ", maximum(Theta1))
end

function plot_sub3_equation(kaps, d_eta)
    Numb_kaps = size(kaps)[1]
    clf()
    x = collect(1.e-6:0.001:0.999)

    # if ((maximum(d_eta)>200) &  (maximum(d_eta)< 500))
    #     x[1] = 2.5e-3
    # elseif (maximum(d_eta)> 500)
    #     x[1] = 1.e-3
    # end

    figure(7,figsize=fig_size)
    plot(x, 1./x + log.(x),linestyle="-",color=colors[6])
    for k = 1:Numb_kaps
        plot(x,d_eta[k]*ones(size(x)),linestyle=":",color="k")
    end

    legend(loc="best",fontsize=legend_fontsize-2)
    ylim(0.98,maximum(d_eta)+0.1*maximum(d_eta))
    xlabel("\$\\Theta\$",fontsize=label_fontsize1+1)
    ylabel("\$\\log(e^{1/\\Theta} \\Theta)\$",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub3_eq.png"))
end
