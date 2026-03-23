% An implementation of the k-bundle Newton method for strongly convex
% functions from [Lewis, Wylie (2019), "A simple Newton method for local
% nonsmooth optimization" (https://arxiv.org/abs/1907.11742)], Alg. 2.2. 

function [x_hat_arr,delta_arr,diam_arr,best_x_arr,best_f_arr] = k_bundle_newton_method(S,problem_data,algo_options)

eps_bar = algo_options.eps_bar;
delta_bar = algo_options.delta_bar;
max_iter = algo_options.max_iter;

[n,k] = size(S);
delta_arr = zeros(1,max_iter);
diam_arr = zeros(1,max_iter);
x_hat_arr = zeros(n,max_iter);
best_f_arr = zeros(1,max_iter);
best_x_arr = zeros(n,max_iter);

% Evaluate objective function and derivatives for initial bundle
f_S = zeros(1,k);
grad_S = zeros(n,k);
hess_S = cell(1,k);
for i = 1:k
    f_S(i) = problem_data.oracle{1}(S(:,i));
    grad_S(:,i) = problem_data.oracle{2}(S(:,i));
    hess_S{i} = problem_data.oracle{3}(S(:,i));
end

% Initialize the best x_hat found so far
best_f = Inf;
best_x = [];

optns = optimoptions("quadprog","Display","none","OptimalityTolerance",10^-15,"StepTolerance",10^-15);
for i = 1:max_iter

    % Compute delta and lambda
    qp_H = grad_S'*grad_S;
    qp_f = zeros(k,1);

    lambda = quadprog(qp_H,qp_f,[],[],ones(1,k),1,zeros(1,k),ones(1,k),[],optns);
    delta = norm(grad_S*lambda,2);
    delta_arr(i) = delta;

    % Stopping criterion
    diam_S = max(pdist(S));
    diam_arr(i) = diam_S;
    if(diam_S < eps_bar && delta < delta_bar )
        x_hat_arr = x_hat_arr(:,1:i-1);
        delta_arr = delta_arr(1:i);
        diam_arr = diam_arr(1:i);
        best_x_arr = best_x_arr(:,1:i-1);
        best_f_arr = best_f_arr(1:i-1);
        break
    end

    % Compute x_hat
    qp_H = zeros(n);
    for s = 1:k
        qp_H = qp_H + lambda(s) * hess_S{s}; % 1/2 already in quadprog
    end

    tmp = zeros(n,1);
    for s = 1:k
        tmp = tmp - lambda(s) * (S(:,s)' * hess_S{s})';
    end
    qp_f = grad_S * lambda + tmp;

    qp_Aeq = (grad_S(:,2:k) - grad_S(:,1))';

    qp_beq = zeros(k-1,1);
    for s = 2:k
        qp_beq(s-1) = -f_S(s) + f_S(1) + grad_S(:,s)'*S(:,s) - grad_S(:,1)'*S(:,1);
    end

    % Solve the quadratic problem either via quadprog or directly (cf.
    % https://en.wikipedia.org/wiki/Quadratic_programming#Equality_constraints)

    % x_hat = quadprog(qp_H,qp_f,[],[],qp_Aeq,qp_beq,[],[],[],optns);
    
    M = [qp_H, qp_Aeq'; qp_Aeq, zeros(k-1)];
    rhs = [-qp_f; qp_beq];
    sol = M^-1 * rhs;
    x_hat = sol(1:n);

    x_hat_arr(:,i) = x_hat;

    % (Skipping nonsmoothness check)

    % Choose s
    grad_x_hat = problem_data.oracle{2}(x_hat);
    Theta_s_arr = zeros(1,k);
    for s = 1:k
        bundle_grad_tmp = grad_S;
        bundle_grad_tmp(:,s) = grad_x_hat;

        qp_H = bundle_grad_tmp'*bundle_grad_tmp;
        qp_f = zeros(k,1);

        lambda = quadprog(qp_H,qp_f,[],[],ones(1,k),1,zeros(1,k),ones(1,k),[],optns);
        Theta_s_arr(s) = norm(bundle_grad_tmp*lambda,2);
    end

    [~,s_min] = min(Theta_s_arr);

    S(:,s_min) = x_hat;
    f_S(s_min) = problem_data.oracle{1}(x_hat);
    grad_S(:,s_min) = grad_x_hat;
    hess_S{s_min} = problem_data.oracle{3}(x_hat);

    if(f_S(s_min) < best_f)
        best_f = f_S(s_min);
        best_x = x_hat;
    end

    best_f_arr(i) = best_f;
    best_x_arr(:,i) = best_x;
end

end

