% Function for solving the subproblem (3.5) in [GU2026a] via fmincon. (Less optimized
% than the method using IPOPT.) 

function [z_bar,theta,mu] = solve_subproblem_fmincon(W,W_f,W_grad,W_hess,x,eps,optns) 

n = size(W_grad,1);

% Define function handles for fmincon
fun = @(in) get_fun(in,n);
nonlcon = @(in) get_nonlcon(in,n,W,W_f,W_grad,W_hess,x,eps);
laghess = @(in,lambda) get_laghess(in,lambda,n,W_hess);

% Set options for fmincon
optns.SpecifyObjectiveGradient = true;
optns.SpecifyConstraintGradient = true;
if(strcmp(optns.Algorithm,'interior-point'))
    optns.HessianFcn = laghess;
end

% Run the solver
[sol,~,exitflag,~,lambda] = fmincon(fun,[x;max(nonlcon([x;0])) + 1],[],[],[],[],[(x - eps);-Inf],[(x + eps);Inf],nonlcon,optns);

% Check the exitflag
if(exitflag < 1)
    fprintf(2,['Warning: fmincon exitflag = ',num2str(exitflag),'.\n'])
end

% Prepare output
mu = lambda.ineqnonlin(end);
z_bar = sol(1:n);
k = numel(W_hess);
nonlcon_sol = nonlcon(sol);
theta = max(nonlcon_sol(1:k) + sol(n+1));

end

% Auxiliary functions
function [y,grady] = get_fun(in,n)
    y = in(n+1);
    grady = [zeros(n,1);1];
end

function [c,ceq,gradc,gradceq] = get_nonlcon(in,n,W,W_f,W_grad,W_hess,x,eps)
    ceq = 0;
    gradceq = zeros(n+1,1);
    
    k = numel(W_hess);
    c = zeros(k+1,1);
    gradc = zeros(n+1,k+1);
    for i = 1:k
        c(i) = W_f(i) + W_grad(:,i)'*(in(1:n) - W(:,i)) + 1/2*(in(1:n) - W(:,i))'*W_hess{i}*(in(1:n) - W(:,i)) - in(n+1);

        gradc(:,i) = [W_grad(:,i) + W_hess{i}*(in(1:n) - W(:,i)); -1];
    end
    c(k+1) = norm(in(1:n) - x,2)^2 - eps^2;
    gradc(:,k+1) = [2*(in(1:n) - x);0];
end

function laghess = get_laghess(in,lambda,n,W_hess)
    
    laghess = zeros(n+1);
    k = numel(W_hess);
    for i = 1:k
        laghess = laghess + lambda.ineqnonlin(i)*[W_hess{i},zeros(n,1);zeros(1,n+1)];
    end
    laghess = laghess + lambda.ineqnonlin(k+1)*[2*eye(n),zeros(n,1); zeros(1,n),0];
end