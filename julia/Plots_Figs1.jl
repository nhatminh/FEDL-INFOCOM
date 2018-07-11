# using Plots
using PyPlot
using PyCall
@pyimport matplotlib.patches as patch

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
markers = ["x","o",">","^", ".","s"]

folder = string("figs//")

function plot_sub1_T(T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, levels)
    clf()
    cfig = figure(1,figsize=fig_size)
    ax = cfig[:add_subplot](1,1,1)
    # plot(kaps,T_cmp,color=colors[1],linestyle="-",linewidth=l_width,label="Solver")
    plot(levels,T_cmp1+0.5,color=colors[6],linestyle="-",linewidth=l_width,label="\$T_{cmp}^*\$")
    plot(levels,Tcmp_N1,color=colors[2],linestyle="-.",linewidth=l_width,label="\$T_{\\mathcal{N}_1}\$")
    plot(levels,Tcmp_N2,color=colors[4],linestyle="-.",linewidth=l_width,label="\$T_{\\mathcal{N}_2}\$")
    plot(levels,Tcmp_N3,color=colors[5],linestyle="-.",linewidth=l_width,label="\$T_{\\mathcal{N}_3}\$")

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("Homogeneous_Sub1_ratio",fontsize=label_fontsize1+1)
    ylabel("\$T_{cmp}\$ (sec)",fontsize=label_fontsize1+1)
    ylim(0, T_cmp1[1]+1.2)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_T_rs.pdf"))
end

function plot_sub1_N(N1, N2, N3, levels)
    clf()
    cfig = figure(2,figsize=fig_size)
    ax = cfig[:add_subplot](1,1,1)
    # plot(kaps,T_cmp,color=colors[1],linestyle="-",linewidth=l_width,label="Solver")
    step(levels,N1,color=colors[4],linestyle="-",linewidth=l_width,label="\$\\mathcal{N}_1\$", where="post", marker=markers[2], markersize=marker_size, markevery=11)
    step(levels,N2,color=colors[3],linestyle="-",linewidth=l_width,label="\$\\mathcal{N}_2\$", where="post", marker=markers[3], markersize=marker_size, markevery=11)
    step(levels,N3,color=colors[2],linestyle="-",linewidth=l_width,label="\$\\mathcal{N}_3\$", where="post", marker=markers[6], markersize=marker_size, markevery=11)

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("Homogeneous_Sub1_ratio",fontsize=label_fontsize1+1)
    ylabel("Three subsets by Alg.1",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_N_rs.pdf"))
end

function plot_sub1_f(f1,levels)
    clf()
    cfig = figure(3,figsize=fig_size)
    ax = cfig[:add_subplot](1,1,1)
    plot(levels,f_min[1]*ones(Numb_D)*1e-9,linestyle=":",color=colors[6])

    if (HETEROGENEOUS != 1 ) # Homogeneous
        plot(levels,f_max[1]*ones(Numb_D)*1e-9,linestyle="--",color=colors[6])
    end

    # for n = 1:5
    #     if (HETEROGENEOUS > 0)  & (abs(f_max[n]*1e-9 - maximum(f1[:,n])) < 1e-3)
    #         plot(D_ratios1,f_max[n]*ones(Numb_D)*1e-9,linestyle="--",color=colors[n])
    #         # plot(kaps,f_min[n]*ones(size(kaps))*1e-9,linestyle=":",color=colors[n])
    #     end
    #
    #     plot(D_ratios1,f1[:,n],color=colors[n],linestyle="-",linewidth=l_width,label=string("UE ",n))
    # end

    for n = 1:NumbDevs
        plot(levels,f1[:,n],linestyle="-",linewidth=l_width)
    end


    # legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("Homogeneous_Sub1_ratio",fontsize=label_fontsize1+1)
    # ylim(0,maximum(f1) + 0.2 )
    ylabel("f (GHz)",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub1_f_rs.pdf"))
end

function plot_sub2_tau(tau1,levels)
    clf()
    cfig = figure(4,figsize=fig_size)
    ax = cfig[:add_subplot](1,1,1)

    # for n = 1:5
    #     plot(tau_ratios,tau1[:,n], color=colors[n], linestyle="-",linewidth=l_width,label=string("UE ",n))
    # end
    for n = 1:NumbDevs
        plot(levels,tau1[:,n],linestyle="-",linewidth=l_width)
    end

    max_tau = maximum(tau1[1,:])

    # legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("Homogeneous_Sub2_ratio",fontsize=label_fontsize1+1)
    ylabel("\$\\tau_n\$ (sec)",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub2_Tau_rs.pdf"))
end

function plot_sub2_p(p1,levels)
    clf()
    cfig = figure(5,figsize=fig_size)
    ax = cfig[:add_subplot](1,1,1)

    plot(levels,Ptx_Max*ones(Numb_Dis),linestyle=":",color=colors[6])
    plot(levels,Ptx_Min*ones(Numb_Dis),linestyle=":",color=colors[6])
    # for n = 1:5
    #     plot(tau_ratios,p1[:,n],color=colors[n],linestyle="-",linewidth=l_width,label=string("UE ",n))
    # end
    for n = 1:NumbDevs
        plot(levels,p1[:,n],linestyle="-",linewidth=l_width)
    end

    # legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("Homogeneous_Sub2_ratio",fontsize=label_fontsize1+1)
    ylabel("p (Watt)",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub2_p_rs.pdf"))
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
    plot(x, glob_numb_iter,linestyle="--",color=colors[3], label=string("\$ 1/(1 - \\Theta)\$"))
    # println(x)
    plot(Theta1[id], Obj1[id],color="r", marker=markers[2], markersize=marker_size)

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\Theta\$",fontsize=label_fontsize1+1)
    # ylabel("Objective",fontsize=label_fontsize1+1)
    yscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub3_obj.pdf"))
    println("Theta: ", minimum(Theta1), " - ", maximum(Theta1))
end

function plot_sub3_kappa_theta(Theta, d_eta,levels, sub)
    Numb_Levels = size(levels)[1]

    clf()
    figure(10,figsize=fig_size)
    # plot(Numb_devs, Objs_E[:,11],linestyle="--",color=colors[1],marker=markers[1], markersize=marker_size, label=string("\$\\kappa\$ =", kaps[11]))
    plot(kaps, 1./d_eta[1,:],linestyle="-",label=string("\$\\eta\$, level",sub,"=",round(levels[1],2)))
    plot(kaps, Theta[1,:],linestyle="-",label=string("\$\\Theta\$, level",sub,"=",round(levels[1],2)))
    plot(kaps, 1./d_eta[4,:],linestyle="-",label=string("\$\\eta\$, level",sub,"=",round(levels[4],2)))
    plot(kaps, Theta[4,:],linestyle="-",label=string("\$\\Theta\$, level",sub,"=",round(levels[4],2)))
    plot(kaps, 1./d_eta[end,:],linestyle="-",label=string("\$\\eta\$, level",sub,"=",round(levels[end],2)))
    plot(kaps, Theta[end,:],linestyle="-",label=string("\$\\Theta\$, level",sub,"=",round(levels[end],2)))

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    # ylabel("\$\\Theta\$",fontsize=label_fontsize1+1)
    xscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"sub3_kappa_theta_rs",sub,".pdf"))
end

function plot_sub3_equation(d_eta,levels)
    Numb_Levels = size(levels)[1]

    clf()
    x = collect(1.e-6:0.001:0.999)

    figure(7,figsize=fig_size)
    plot(x, 1./x + log.(x),linestyle="-",color=colors[6])
    for k = 1:Numb_kaps
        plot(x,d_eta[k]*ones(size(x)),linestyle=":",color="k")
    end

    legend(loc="best",fontsize=legend_fontsize-2)
    xlim(0, 0.5)
    ylim(0.98,maximum(d_eta)+0.1*maximum(d_eta))
    xlabel("\$\\Theta\$",fontsize=label_fontsize1+1)
    ylabel("\$\\log(e^{1/\\Theta} \\Theta)\$",fontsize=label_fontsize1+1)
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"Sub3_eq.pdf"))
end

function plot_numerical_pareto(Theta1, T_cmp1, E_cmp1, T_com1, E_com1, levels, sub)
    Numb_Levels = size(levels)[1]
    if sub == 1
        # idx_levels = [1, 9]
        idx_levels = 1:9
    else
        # idx_levels = [1, 9]
        idx_levels = 1:8
    end
    clf()
    figure(9,figsize=fig_size)
    for idx in idx_levels
        E_obj   = zeros(Numb_kaps)
        T_obj   = zeros(Numb_kaps)

        for i=1:Numb_kaps
            E_obj[i] = 1/(1 - Theta1[idx,i])* (E_com1[idx,i] - log(Theta1[idx,i])*E_cmp1[idx,i])
            T_obj[i] = 1/(1 - Theta1[idx,i])* (T_com1[idx,i] - log(Theta1[idx,i])*T_cmp1[idx,i])
        end
        scatter(E_obj, T_obj,label=string("level",sub,"=",round(levels[idx],2)))
    end

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("E_obj",fontsize=label_fontsize1+1)
    ylabel("T_obj",fontsize=label_fontsize1+1)
    # yscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"pareto_rs",sub,".pdf"))
end

function plot_total_cost(Obj, levels, sub)
    Numb_Levels = size(levels)[1]
    clf()
    figure(9,figsize=fig_size)
    plot(kaps, Obj[1,:],linestyle="-",label=string("level",sub,"=",round(levels[1],2)))
    plot(kaps, Obj[4,:],linestyle="-",label=string("level",sub,"=",round(levels[4],2)))
    plot(kaps, Obj[end,:],linestyle="-",label=string("level",sub,"=",round(levels[end],2)))

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel("\$\\kappa\$",fontsize=label_fontsize1+1)
    ylabel("Total Cost",fontsize=label_fontsize1+1)
    xscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"total1_rs",sub,".pdf"))

    clf()
    figure(10,figsize=fig_size)
    plot(levels, Obj[:,1],linestyle="-",label=string("\$\\kappa\$=",kaps[1]))
    plot(levels, Obj[:,10],linestyle="-",label=string("\$\\kappa\$=",kaps[10]))
    plot(levels, Obj[:,end],linestyle="-",label=string("\$\\kappa\$=",kaps[end]))

    legend(loc="best",fontsize=legend_fontsize-2)
    xlabel(string("level",sub),fontsize=label_fontsize1+1)
    ylabel("Total Cost",fontsize=label_fontsize1+1)
    # yscale("log")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"total2_rs",sub,".pdf"))
end

function save_result(filename,Theta1, Obj1, Obj_E, Obj_T, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta, levels)
    h5open(filename, "w") do file
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
        write(file,"levels", levels)
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
        levels = read(file,"levels")
        return Theta1, Obj1, Obj_E, Obj_T, T_cmp1, Tcmp_N1, Tcmp_N2, Tcmp_N3, E_cmp1, T_com1, E_com1, N1, N2, N3, f1, tau1, p1, d_eta, levels
    end
end
