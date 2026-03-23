%{ 
An implementation of Alg. 4.2 from [GU2026a]. Models of order q > 2 are
currently not supported. 

Inputs:
    x1: Initial point
    eps1: Initial trust-region radius
    problem_data: A struct with fields:
        n: Number of variables
        oracle: A cell array in which the i-th cell contains a function
            handle to evaluate the (i-1)-th derivative. For example,
            oracle{1}(x) is the objective value and oracle{2}(x) is the
            gradient at a point x.
    algo_options: A struct with fields:
        q: Order of Taylor expansion to be used
        p: Order of growth of the objective function around the minimum
        sigma: Approximation parameter for Alg. 4.1
        disp_flag: Flag for displaying information of iterations
        memory_max_size: Maximum size of memory of past iterates for the
            bundle initialization in Alg. 4.1  (cf. Remark 4.2(a)).
        local_options: A struct with fields:
            kappa: Rate parameter in (4.2)
            eps_thr: Lower bound for eps (cf. Section 6)
            j_thr: For detecting unsuccessful runs of Alg. 4.2 (cf.
                Subsection 5.1) 
            act_thr: Tolerance for activity of trust-region constraint
            init_N_sample: Number of randomly sampled points for the
                initialization of Alg. 4.1. The first point is always the
                center of the trust-region, so init_N_sample = 1 means that
                no additional points are randomly sampled. (Always set to 1
                in our experiments.)   
            norm_flag: Which p-norm to use for the trust-region. (Not to be
                confused with the growth order p.)
            sp_solver: Which solver to use for subproblem (3.5). Currently
                supports "IPOPT" and "fmincon" (for q = 2 and norm_flag =
                2) and "linprog" (for q = 1 and norm_flag = inf). ("IPOPT"
                uses a slightly modified version of mexIPOPT, see the
                comments in "solve_subproblem_IPOPT.m" for details.)  
            sp_solver_optns: Can be used to set solver specific options

Output:
    result_local_method: A struct with fields:
        x_arr: Sequence (x^j)_j in Alg. 4.2
        f_arr: Sequence (f(x^j))_j
        eps_arr: The trust-region for each j (as in (4.1)) 
        act_arr: The relative activity of the trust-region constraint
        mu_arr: The multiplier of the trust-region constraint in subproblem
            (3.5) for each j
        numsample_arr: The number of oracle calls in Alg. 4.1 for each j
        best_x: Point with the smallest objective value that was
            encountered. (Recall that Alg. 4.2 is not a descent method.)
        best_f_val: Objective value at best_x
        eval_counter: Number of oracle calls for each order
        eval_counter_arr: Number of oracle calls for each order that was
            required for each j. (The last row is eval_counter.)
        success_flag: Returns the success of the run (cf. Subsection 5.1)
%} 

function result_local_method = local_method(x1,eps1,problem_data,algo_options)

% Read problem specification
n = problem_data.n;
oracle = problem_data.oracle;

% Read algorithm options
q = algo_options.q;
p = algo_options.p;
sigma = algo_options.sigma;
disp_flag = algo_options.disp_flag;
memory_max_size = algo_options.memory_max_size;

local_options = algo_options.local_options;
kappa = local_options.kappa;
eps_thr = local_options.eps_thr;
j_thr = local_options.j_thr;
act_thr = local_options.act_thr;

norm_flag = local_options.norm_flag;

% Tic for measuring runtime
if(disp_flag >= 1)
    start_tic = tic;
end

% Initialize variables
j_max = floor(1 + (log(log(eps_thr/eps1)/log(kappa) + 1))/log((q+sigma)/p)); % (solving eps_j = eps_thr for j)
x_arr = zeros(n,j_max+1);
f_arr = zeros(1,j_max+1);
eps_arr = zeros(1,j_max);
act_arr = zeros(1,j_max);
mu_arr = zeros(1,j_max);
numsample_arr = zeros(1,j_max);
eval_counter = zeros(1,q+1);
eval_counter_arr = zeros(j_max+1,q+1);
success_flag = true;

% Initialize the memory
memory.sample_pts = [];
memory.oracle_vals = {};
memory.max_size = memory_max_size;

x_arr(:,1) = x1;
f_x = oracle{1}(x1); eval_counter(1) = eval_counter(1) + 1;
f_x1 = f_x;
f_arr(1) = f_x;

%%%%%%%%%%%%%%%%%%
%%%%% Step 1 %%%%%
%%%%%%%%%%%%%%%%%%

for j = 1:j_max
    
    %%%%%%%%%%%%%%%%%%
    %%%%% Step 2 %%%%%
    %%%%%%%%%%%%%%%%%%

    eps_arr(j) = eps1 * kappa.^(((q + sigma)/p).^(j-1) - 1);

    if(disp_flag >= 1)
        disp(['            ----- Iteration j = ',num2str(j),'/',num2str(j_max),' -------------------'])
        disp(['                eps_j = ',num2str(eps_arr(j))]);
    end

    %%%%%%%%%%%%%%%%%%
    %%%%% Step 3 %%%%%
    %%%%%%%%%%%%%%%%%%

    if(disp_flag >= 1)
        disp('                Applying Alg. 4.1...')
    end

    % Execute Alg. 4.1 for generating the set W
    [z_bar,f_z_bar,mu,eval_counter,memory] = generate_W(x_arr(:,j),eps_arr(j),f_x,Inf,0,memory,eval_counter,problem_data,algo_options,local_options);

    % Save number of oracle calls up to this iteration
    eval_counter_arr(j+1,:) = eval_counter;

    % Save number of samples required
    numsample_arr(j) = eval_counter_arr(j+1,2) - eval_counter_arr(j,2);

    % Save multiplier of trust-region constraint
    mu_arr(j) = mu;

    % Save activity of trust-region constraint
    act_arr(j) = norm(z_bar - x_arr(:,j),norm_flag)/eps_arr(j);

    if(disp_flag >= 1)
        disp(['                activity   = ',num2str(act_arr(j))]);
        disp(['                act_thresh = ',num2str(act_thr)]);
        disp(['                j     = ',num2str(j)])
        disp(['                j_thr = ',num2str(j_thr)])
    end

    %%%%%%%%%%%%%%%%%%
    %%%%% Step 4 %%%%%
    %%%%%%%%%%%%%%%%%%

    % Update iterate
    x_arr(:,j+1) = z_bar;
    f_arr(j+1) = f_z_bar;
    f_x = f_z_bar;

    % Check for nonsuccess (cf. Subsection 5.1)
    if(j >= j_thr && act_arr(j) > act_thr)
        x_arr = x_arr(:,1:j+1);
        f_arr = f_arr(1:j+1);
        act_arr = act_arr(1:j);
        mu_arr = mu_arr(1:j);
        numsample_arr = numsample_arr(1:j);
        success_flag = false;

        break
    end

end

% Find point with smallest objective value. (Recall that Alg. 4.2 is not a
% descent method.) 
[best_f_val,I_best_f_val] = min(f_arr);
best_x = x_arr(:,I_best_f_val);

if(disp_flag >= 1)
    if(success_flag)
        disp('        Local method successful!')
    else
        disp('        Local method failed!')
    end

    disp(['            f(x1)    = ',num2str(f_x1)]);
    disp(['            best_f   = ',num2str(best_f_val)]);
    disp(['            f(z_bar) = ',num2str(f_z_bar)]);

    disp(' ')
    disp(['            Runtime = ',num2str(toc(start_tic)),'s']);
    disp('            Required evaluations: ');
    disp(['                f    - ',num2str(eval_counter(1))])
    disp(['                grad - ',num2str(eval_counter(2))])
    if(q > 1)
        disp(['                hess - ',num2str(eval_counter(3))])
    end

end

result_local_method.x_arr = x_arr;
result_local_method.f_arr = f_arr;
result_local_method.eps_arr = eps_arr;
result_local_method.act_arr = act_arr;
result_local_method.mu_arr = mu_arr;
result_local_method.numsample_arr = numsample_arr;
result_local_method.best_x = best_x;
result_local_method.best_f_val = best_f_val;
result_local_method.eval_counter = eval_counter;
result_local_method.eval_counter_arr = eval_counter_arr;
result_local_method.success_flag = success_flag;

end

