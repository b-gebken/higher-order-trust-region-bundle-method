using NonSmoothProblems
using LinearAlgebra
using MAT

import NonSmoothProblems.F
import NonSmoothProblems.∂F_elt

struct Halfhalf_mod{Tf} <: NonSmoothPb
    n::Int64
    A::Diagonal{Tf}
    B::Diagonal{Tf}
    Halfhalf_mod(;Tf = Float64) = new{Tf}(
        8,
        Diagonal(Tf[isodd(i) ? 1 : 0 for i in 1:8]),
        Diagonal([Tf(1)/Tf(i)^2 for i in 1:8])
    )
end
projectodd(x) = @view x[1:2:7]

f(::Halfhalf_mod, y) = norm(y)
g(::Halfhalf_mod, x) = projectodd(x)

function F(pb::Halfhalf_mod, x) 

    global F_eval_x_arr, F_eval_fx_arr
    
    Fx = f(pb, g(pb, x)) + dot(x, pb.B*x)

    F_eval_x_arr = cat(F_eval_x_arr,x;dims = 2)
    F_eval_fx_arr = cat(F_eval_fx_arr,Fx;dims = 2)

    return Fx
end

function ∂F_elt(pb::Halfhalf_mod, x)
    
    global gradF_eval_x_arr, gradF_eval_fx_arr

    gradF_eval_x_arr = cat(gradF_eval_x_arr,x;dims = 2)
    gradF_eval_fx_arr = cat(gradF_eval_fx_arr,f(pb, g(pb, x)) + dot(x, pb.B*x);dims = 2)

    Ax = pb.A * x
    if norm(Ax) == 0
        return 2*pb.B*x
    else
        return Ax ./ norm(Ax) + 2*pb.B*x
    end
end

is_differentiable(::Halfhalf_mod, x) = (norm(projectodd(x)) != 0)