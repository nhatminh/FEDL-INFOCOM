# OnDevAI

This repository is a simulation code for the paper

Nguyen H. Tran, Wei Bao, Albert Zomaya, Minh N. H. Nguyen, Choong Seon Hong 
"Federated Learning over Wireless Networks: Optimization Model Design and Analysis", IEEE INFOCOM 2019 (Accepted).

We use Julia version 0.6.

All of Figures will be generated in `figs` folder.

Main.jl will includes 5 files
- "Setting.jl"
- "TrafficGen.jl"
- "Solving1.jl"
- "Utils.jl" 
- "Plots_Figs.jl"/"Plots_Figs1.jl"

```
Several Setting config parameters as 
if(NUMERICAL_RS)
    # Result for Section V: NUMERICAL RESULTS
    include("Plots_Figs1.jl")
else
    # Result for Section IV: Closed-Form solution
    include("Plots_Figs.jl")
end


READ_RESULT = true		# Read results from result files
NUMERICAL_RS = true     # Result for Section V: NUMERICAL RESULTS in the paper (50 devs)
# NUMERICAL_RS = false  # Result for Section IV: Closed-Form solution in the paper (5 devs)
```

We apply the Federated Learning code from
https://github.com/shaoxiongji/federated-learning
