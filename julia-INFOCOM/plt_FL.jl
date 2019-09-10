using PyPlot
using PyCall
@pyimport matplotlib.patches as patch
fig_size = (6.1,4.3)
folder = string("figs//")
label_fontsize = 18-1.5
legend_fontsize = label_fontsize - 2
patterns = ["","."]

label_fontsize1 = label_fontsize
marker_size=6
l_width=1.2
# stride = convert(Int, max_iters/10) +1
stride = 6

colors=["m","b","coral","g","k","r"]
markers = ["x","o",">","^", "s","."]

local_iters = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
global_iters = [227, 121, 86, 69, 59, 51, 34, 30, 27, 23, 21, 20, 18]
Normalized_Local = 50
Normalized_Global = 40
theory_theta2 = [ 1/(exp(x/Normalized_Local)) for x in local_iters ]
println("Theta: ",theory_theta2)
println("Theory numb of global iters: ",[ 1/(1-x) for x in theory_theta2 ])
println("Normlized true numb of global iters: ",[ x/Normalized_Global for x in global_iters ])


local_iters1 = [5, 10, 15, 20, 25, 30, 35, 40]
global_iters_FedProx = [152, 48, 20, 12, 8, 6, 5, 5]
global_iters_FedAvg = [141, 39, 20, 12, 8, 6, 5, 4]

function plot_FL_MINIST_comparison()
    clf()
    cfig = figure(2,figsize=fig_size)
    ax = cfig[:add_subplot](1,1,1)
    ax[:tick_params]("both",labelsize=legend_fontsize-1)

    step(local_iters1,global_iters_FedProx, where="post",label="FedProx",linestyle="--", marker=markers[2], markersize=marker_size, markevery=1)
    step(local_iters1,global_iters_FedAvg, where="post",label="FedAvg",linestyle="--", marker=markers[3], markersize=marker_size, markevery=1)
    legend(loc="best",fontsize=legend_fontsize)
    xlabel("Number of Local Iterations",fontsize=label_fontsize-1)
    ylabel("Number of Global Iterations",fontsize=label_fontsize-1)
    # title("Relationship between the number of local and global iterations")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"local_global_iters_comparison.pdf"))
    show()
end

function plot_FL_MINIST_5users(opt_theta)
    println("Theta = ", opt_theta, ", Opt Numb local iters: ", round(log(1/opt_theta)*Normalized_Local/10)*10)

    clf()
    cfig = figure(2,figsize=fig_size)

    step(local_iters,global_iters, where="post")
    xlabel("Numb of Local Iterations(\$x_{loc}\$)")
    ylabel("Numb of Global Iterations (\$K(\\Theta)\$)")
    title("Relationship between the number of local and global iterations")
    tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    savefig(string(folder,"local_global_iters.pdf"))
    show()

    clf()
    cfig = figure(2,figsize=fig_size)
    # max_local = maximum(local_iters)
    # local_iters1 = collect(10:10:200)
    # theory_theta1 = [ 1/(exp(x/190)) for x in local_iters1 ]


    step(local_iters,theory_theta2, where="post")
    xlabel("Numb of Local Iterations (\$x_{loc}\$)")
    ylabel("Local Error (\$\\Theta=\\frac{1}{e^{x_{loc}/50}} \$)")
    title("Relationship between the number of local iters and \$\\Theta\$")
    savefig(string(folder,"local_theta.pdf"))
    show()

    clf()
    cfig = figure(2,figsize=fig_size)
    step(global_iters, theory_theta2, where="post")
    xlabel("Numb of Global Iterations (\$K(\\Theta)\$)")
    ylabel("Local Error (\$\\Theta\$)")
    title("Relationship between the number of global iters and \$\\Theta\$")
    savefig(string(folder,"global_theta.pdf"))
    show()
end

# opt_theta = 0.405
# plot_FL_MINIST_5users(opt_theta)
plot_FL_MINIST_comparison()
