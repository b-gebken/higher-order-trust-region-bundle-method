# Script which applies SuperPolyak from https://github.com/COR-OPT/SuperPolyak.jl to the function (6.1).
# (Install via Pkg.add(url="https://github.com/COR-OPT/SuperPolyak.jl"))

using LinearAlgebra
using SuperPolyak
using MAT
using Dates

# Flag for saving the result of SuperPolyak
save_flag = true

include("max_root.jl")

n = 100
f = max_root
gradf = grad_max_root
x0 = 0.1*collect(1:n)/n

# Run SuperPolyak
result = SuperPolyak.superpolyak(f,gradf,x0)

# Run SuperPolyak another 20 times (without saving the result) to get the average runtime
reps = 20
total_time = 0
for i in 1:reps
    global total_time += @elapsed SuperPolyak.superpolyak(f,gradf,x0)
end
println("Avg. time = ",total_time/reps)

# Print results
x_min = zeros(n)
f_min = 0

num_iter = size(result.loss_history,1)-1
println("||x - x^*|| = ",norm(result.solution - x_min,2))
println("#iter       = ",num_iter)
println("#subgrad    = ",sum(result.oracle_calls))

# Export result to mat-file
if(save_flag)
    file = matopen("SuperPolyak-result-" * Dates.format(Dates.now(), dateformat"dd-mm-yy_HH-MM-SS") * ".mat", "w")
    write(file, "oracle_calls_SupPol", result.oracle_calls)
    write(file, "loss_history_SupPol", result.loss_history)
    write(file, "solution_SupPol", result.solution)
    close(file)
end