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

function plot_sub1(kaps, T_cmp, T_cmp1)
    clf()
    figure(1,figsize=fig_size)
    plot(kaps,T_cmp,color=colors[1],linestyle="-",linewidth=l_width,label="Solver")
    plot(kaps,T_cmp1,color=colors[2],linestyle=":",linewidth=l_width,label="Closed-Form")

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    xscale("log")
    ylabel("\$T_{cmp}\$",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1.png"))
end


function plot_sub2_tau(kaps, tau, tau1)
    clf()
    figure(2,figsize=fig_size)
    for n = 1:NumbDevs
        plot(kaps,tau1[:,n], color=colors[n], linestyle="-",linewidth=l_width,label=string("Dev ",n))
    end
    # plot(kaps,T_comps1,marker=markers[2],color=colors[2], markersize=marker_size,linestyle=":",linewidth=l_width,label="Tcmp1")

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    xscale("log")
    ylabel("\$\\tau_n\$",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub2a.png"))
end

function plot_sub2_p(kaps, p, p1)
    clf()
    figure(3,figsize=fig_size)
    plot(kaps,Ptx_Max*ones(size(kaps)),linestyle=":",color=colors[6])
    plot(kaps,Ptx_Min*ones(size(kaps)),linestyle=":",color=colors[6])
    for n = 1:NumbDevs
        plot(kaps,p1[:,n],color=colors[n],linestyle="-",linewidth=l_width,label=string("Dev ",n))
    end
    # plot(kaps,T_comps1,marker=markers[2],color=colors[2], markersize=marker_size,linestyle=":",linewidth=l_width,label="Tcmp1")

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    xscale("log")
    ylabel("p",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub2b.png"))
end

function plot_sub3_cvx(kaps, Theta, Theta1, Obj, Obj1)
    clf()
    figure(4,figsize=fig_size)
    plot(Theta,Obj,color=colors[1],linestyle="-",linewidth=l_width,label="Solver")
    plot(Theta1,Obj1,color=colors[2],linestyle=":",linewidth=l_width,label="Closed-Form")

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\Theta\$",fontsize=label_fontsize1+1)
    ylabel("Objective",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub3_obj.png"))
end

function plot_sub3_equation(kaps, d_eta)
    Numb_kaps = size(kaps)[1]
    clf()
    x = collect(1e-6:0.0001:0.9999)
    # println(x)
    figure(5,figsize=fig_size)
    plot(x,log.(e.^(1./x).*x),linestyle="-",color=colors[6])
    # println(x)
    println(log(e^(1/x[1])x[1]))
    for k = 1:Numb_kaps
        plot(x,d_eta[k]*ones(size(x)),linestyle=":",color="k")
    end

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\Theta\$",fontsize=label_fontsize1+1)
    ylabel("\$\\log(e^{1/\\Theta} \\Theta)\$",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub3_eq.png"))
end
