using Roots
fx(x) = log(x)/x - 0.1
# ## bracketing
# fzero(fx, 8, 9)		          # 8.613169456441398
# fzero(fx, -10, 0)		      # -0.8155534188089606
# fzeros(fx, -10, 10)            # -0.815553, 1.42961  and 8.61317

## use a derivative free method
fzero(fx, 1)			          # 1.4296118247255558

# ## use a different order
# fzero(sin, 3, order=16)		  # 3.141592653589793
