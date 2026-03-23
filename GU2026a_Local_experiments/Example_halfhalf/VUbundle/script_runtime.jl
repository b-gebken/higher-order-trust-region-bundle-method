# Script which applies the VUbundle algorithm from https://github.com/GillesBareilles/NonSmoothSolvers.jl to the half-and-half problem and measures the runtime.

using NonSmoothSolvers
using NonSmoothProblems

n = 8
pb = Halfhalf()
x = zeros(8) .+ 20.08

optparams = OptimizerParams(show_trace = false)

# Run VU bundle method 20 times to get the average runtime
o = VUbundle()

reps = 20
total_time = 0
for i in 1:reps
    global total_time += @elapsed xfinal_vu, tr = optimize!(pb, o, x; optparams);
end

println("Avg. time = ",total_time/reps)
