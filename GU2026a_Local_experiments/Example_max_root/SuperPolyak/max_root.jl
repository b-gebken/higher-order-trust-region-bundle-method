function max_root(x::Vector{Float64})
    a = 1/4
    return maximum(sqrt.(abs.(x) .+ a) .- sqrt(a))
end

function grad_max_root(x::Vector{Float64})

    a = 1/4
    (~,I) = findmax(sqrt.(abs.(x) .+ a) .- sqrt(a)) 

    n = size(x,1)
    grad = zeros(n);
    
    s = sign(x[I])
    if(s == 0)
        s = 1
    end
    grad[I] = 1/2*(abs(x[I]) + a)^(-1/2) * s

    return grad
end