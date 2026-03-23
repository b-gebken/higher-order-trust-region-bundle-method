% Function for solving the subproblem (3.5) in [GU2026a] via linprog.

function [z_bar,theta,mu] = solve_subproblem_linprog(W,W_f,W_grad,x,eps,optns)

    % Define input for linprog
    [n,N_sample] = size(W);
    lp_fun = [zeros(n,1); 1];
    lp_A = [W_grad',-ones(N_sample,1)];
    lp_b = zeros(N_sample,1);
    for i = 1:N_sample
        lp_b(i) = -W_f(i) + W_grad(:,i)'*W(:,i);
    end

    lp_lb = [x - eps; -Inf];
    lp_ub = [x + eps; Inf];

    % Run the solver
    [lp_sol,~,~,~,lambda] = linprog(lp_fun,lp_A,lp_b,[],[],lp_lb,lp_ub,optns);

    % Check if solution was found
    if(isempty(lp_sol) || lambda.lower(n+1) > 0 || lambda.upper(n+1) > 0)
        fprintf(2,'No solution via linprog!\n')
        z_bar = Inf(n,1);
        theta = Inf;
        mu = Inf;
    else
        z_bar = lp_sol(1:n);
        theta = max(lp_A*lp_sol - lp_b + lp_sol(n+1));
        mu = max([lambda.upper;lambda.lower],[],1);
    end
    
end