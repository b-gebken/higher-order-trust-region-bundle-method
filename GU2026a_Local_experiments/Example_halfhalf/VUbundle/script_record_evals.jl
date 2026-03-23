# Script which applies the VUbundle algorithm from https://github.com/GillesBareilles/NonSmoothSolvers.jl to the half-and-half problem and records all points at which the oracle was called.
# The points are recorded via global variables in the modified implementation Halfhalf_mod.jl of the half-and-half function.
# (For measuring the runtime, the unmodified implementation of half-and-half from NonSmoothProblems is used in the script "script_runtime".)

using NonSmoothSolvers
using NonSmoothProblems
using MAT
using Dates
using DataStructures

# Flag for saving the result of VUbundle
save_flag = true

include("Halfhalf_mod.jl")

n = 8

# Global variables for recording all oracle calls
global F_eval_x_arr = zeros(n,0)
global F_eval_fx_arr = zeros(1,0)
global gradF_eval_x_arr = zeros(n,0)
global gradF_eval_fx_arr = zeros(1,0)

pb = Halfhalf_mod()
x = zeros(8) .+ 20.08

# Run VU bundle method
o = VUbundle()
xfinal_vu, tr = optimize!(pb, o, x);

# Process output
num_iter = length(tr)
f_arr_VU = zeros(num_iter)
cumul_nullsteps_arr_VU = zeros(num_iter)
for i in 1:num_iter
    f_arr_VU[i] = tr[i].Fx;

    if(i >= 2)
        cumul_nullsteps_arr_VU[i] = cumul_nullsteps_arr_VU[i-1] + tr[i].additionalinfo.nnullsteps
    end
end

# Export result to mat-file
if(save_flag)
    file = matopen("VUbundle-result-" * Dates.format(Dates.now(), dateformat"dd-mm-yy_HH-MM-SS") * ".mat", "w")
    write(file, "F_eval_x_arr", F_eval_x_arr)
    write(file, "F_eval_fx_arr", F_eval_fx_arr)
    write(file, "gradF_eval_x_arr", gradF_eval_x_arr)
    write(file, "gradF_eval_fx_arr", gradF_eval_fx_arr)
    write(file, "f_arr_VU", f_arr_VU)
    write(file, "cumul_nullsteps_arr_VU", cumul_nullsteps_arr_VU)
    close(file)
end