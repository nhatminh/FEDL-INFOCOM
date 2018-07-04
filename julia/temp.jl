# using Roots
# fx(x) = log(x)/x - 0.1
# # ## bracketing
# # fzero(fx, 8, 9)		          # 8.613169456441398
# # fzero(fx, -10, 0)		      # -0.8155534188089606
# # fzeros(fx, -10, 10)            # -0.815553, 1.42961  and 8.61317
#
# ## use a derivative free method
# fzero(fx, 1)			          # 1.4296118247255558
#
# # ## use a different order
# # fzero(sin, 3, order=16)		  # 3.141592653589793
#
# println(log(e))
#
# UEs = Dict{Float64,Float64}()
# UEs[1] = 2
# UEs[2] = 1
# sorted_UEs_array = sort(collect(UEs), by=x->x[2])
# println(size(sorted_UEs_array)[1])
# # for k,v in size(sorted_UEs_array)
# #     print(k)
# # end
#
# a = 1:3
# println(typeof(a))


# using Distributions
# println(rand(Uniform(1, 2), 5))

# using PyCall
# pygui(:qt5)

using PyPlot
clf()
figure(1)
x = linspace(0,2*pi,1000);
 y = sin.(5*x + 4*cos.(2*x))
plot(x, y, color="red", linewidth=2.0, linestyle="--")
show()
savefig("plot.png")

# using Plots
# pyplot() # Switch to using the PyPlot.jl backend
# plot(rand(5,5),linewidth=2,title="My Plot") # The same plotting command works

x = collect(0:.1:10)
println(x.^2)
