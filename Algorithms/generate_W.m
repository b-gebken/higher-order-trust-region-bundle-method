% An implementation of Alg. 4.1 from [GU2026a], including the modification
% from [GU2026b] of adding the parameter c. Models of order q > 2 are
% currently not supported.  

function [z_bar,f_z_bar,mu,eval_counter,memory,numsample] = generate_W(x,eps,f_x,c,reusing_eps_tolerance,memory,eval_counter,problem_data,algo_options,sp_options)

% Read problem specification
n = problem_data.n;
oracle = problem_data.oracle;

% Read algorithm options
q = algo_options.q;
sigma = algo_options.sigma;
disp_flag = algo_options.disp_flag;

init_N_sample = sp_options.init_N_sample;
norm_flag = sp_options.norm_flag;
sp_solver = sp_options.sp_solver;
sp_solver_optns = sp_options.sp_solver_optns;

% For computing numsample in the end
eval_counter_old = eval_counter;

% Tolerance for reusing points
if(eps < 2*reusing_eps_tolerance)
    fprintf(2,['Warning: eps in magnitude of reuse tolerance (eps = ',num2str(eps),', reuse tolerance = ',num2str(reusing_eps_tolerance),').\n'])
end

% Reuse memorized elements
if(size(memory.sample_pts,2) > 0)
    inds = vecnorm(memory.sample_pts - x,norm_flag,1) <= eps + reusing_eps_tolerance;

    W = memory.sample_pts(:,inds);
    oracle_W = memory.oracle_vals(:,inds);

    N_reused = size(W,2);

    if(disp_flag >= 2)
        disp(['                    Reusing ',num2str(N_reused),' sample points. (Memory size = ',num2str(size(memory.sample_pts,2)),'.)'])
    end
else
    W = [];
    oracle_W = {};

    N_reused = 0;
end

% Sample new elements in the trust-region
if(disp_flag >= 2)
    disp(['                    Initially sampling ',num2str(init_N_sample),' point(s).'])
end

W(:,end+1:end+init_N_sample) = x + eps*sample_hypersphere(n,init_N_sample);
oracle_W(:,end+1:end+init_N_sample) = cell(q+1,init_N_sample);
for k = 1:init_N_sample
    if(k == 1)
        oracle_W{1,N_reused + k} = f_x;
    else
        oracle_W{1,N_reused + k} = oracle{1}(W(:,N_reused + k)); eval_counter(1) = eval_counter(1) + 1;
    end

    for order = 2:q+1
        oracle_W{order,N_reused + k} = oracle{order}(W(:,N_reused + k)); eval_counter(order) = eval_counter(order) + 1;
    end
end
N_W = N_reused + init_N_sample;

%%%%%%%%%%%%%%%%%%
%%%%% Step 1 %%%%%
%%%%%%%%%%%%%%%%%%

while(1)

    %%%%%%%%%%%%%%%%%%
    %%%%% Step 2 %%%%%
    %%%%%%%%%%%%%%%%%%

    % Solve subproblem (3.5)
    if(strcmp(sp_solver,'fmincon'))
        [z_bar,theta,mu] = solve_subproblem_fmincon(W,cell2mat(oracle_W(1,:)),cell2mat(oracle_W(2,:)),oracle_W(3,:),x,eps,sp_solver_optns);
    elseif(strcmp(sp_solver,'IPOPT'))
        [z_bar,theta,mu] = solve_subproblem_IPOPT(W,cell2mat(oracle_W(1,:)),cell2mat(oracle_W(2,:)),oracle_W(3,:),x,eps,sp_solver_optns);
    elseif(strcmp(sp_solver,'linprog'))
        [z_bar,theta,mu] = solve_subproblem_linprog(W,cell2mat(oracle_W(1,:)),cell2mat(oracle_W(2,:)),x,eps,sp_solver_optns);
    end

    if(norm(z_bar - x,norm_flag) > 1.5*eps)
        fprintf(2,['Warning: eps-constraint is heavily violated! eps = ',num2str(eps), ', ||z-x|| = ',num2str(norm(z_bar - x,norm_flag)),', frac = ',num2str(norm(z_bar - x,norm_flag)/eps),'\n'])
    end

    f_z_bar = oracle{1}(z_bar); eval_counter(1) = eval_counter(1) + 1;

    %%%%%%%%%%%%%%%%%%
    %%%%% Step 3 %%%%%
    %%%%%%%%%%%%%%%%%%
    
    if(f_z_bar - theta <= min(c,eps^(q+sigma)))

        %%%%%%%%%%%%%%%%%%
        %%%%% Step 4 %%%%%
        %%%%%%%%%%%%%%%%%%

        break
    else

        %%%%%%%%%%%%%%%%%%
        %%%%% Step 6 %%%%%
        %%%%%%%%%%%%%%%%%%

        % Sample new information at z_bar
        N_W = N_W + 1;
        W = [W,z_bar];
        oracle_W{1,N_W} = f_z_bar;
        for order = 2:q+1
            oracle_W{order,N_W} = oracle{order}(z_bar); eval_counter(order) = eval_counter(order) + 1;
        end
    end
end

if(disp_flag >= 2)
    disp('                ...Alg. 4.1 finished!')
    disp(['                ',num2str(N_W - init_N_sample - N_reused),' additional sample points required.']);
    disp(['                f(z_bar) = ',sprintf('%9.2e', f_z_bar)]);
    disp(['                theta    = ',sprintf('%9.2e', theta)]);
    disp(['                f(z_bar) - theta  = ',sprintf('%9.2e', f_z_bar - theta)]);
    disp(['                eps^(q+sigma)     = ',sprintf('%9.2e', eps^(q+sigma))]);
    disp(' ');
end

% Store newly sampled jet elements
if(memory.max_size > 0)
    memory = add_to_memory(W(:,N_reused + 1 : end),...
        oracle_W(:,N_reused + 1 : end),...
        memory);
end

numsample = eval_counter(2) - eval_counter_old(2);

end

