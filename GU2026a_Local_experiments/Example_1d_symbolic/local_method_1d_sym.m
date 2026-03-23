% A modified version of local_method.m using symbolic arrays. Only intended
% for use in Example 6.1!

function result_local_method = local_method_1d_sym(x1,eps1,problem_data,algo_options)

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

% Initialize arrays
j_max = double( floor(1 + (log(log(eps_thr/eps1)/log(kappa) + 1))/log((q+sigma)/p)) ); % (solving eps_j = eps_thr for j)
x_arr = vpa(zeros(n,j_max+1));
f_arr = vpa(zeros(1,j_max+1));
eps_arr = vpa(zeros(1,j_max));
act_arr = vpa(zeros(1,j_max));
mu_arr = vpa(zeros(1,j_max));
numsample_arr = NaN(1,j_max);
eval_counter = zeros(1,double(q)+1);
eval_counter_arr = zeros(j_max+1,double(q)+1);
success_flag = true;
syms z;

% Initialize the memory
memory.sample_pts = [];
memory.oracle_vals = {};
memory.max_size = memory_max_size;

% Tic for measuring runtime
if(disp_flag >= 1)
    start_tic = tic;
end

x_arr(:,1) = x1;
f_x = oracle{1}(x1); eval_counter(1) = eval_counter(1) + 1;
f_x1 = f_x;
f_arr(1) = f_x;

%%%%%%%%%%%%%%%%%%
%%%%% Step 1 %%%%%
%%%%%%%%%%%%%%%%%%

for j = 1:j_max

    if(disp_flag >= 1)
        disp(['            ----- Iteration j = ',num2str(j),'/',num2str(j_max),' -------------------'])
        disp(['                eps_j = ',num2str(double(eps_arr(j)))]);
    end

    %%%%%%%%%%%%%%%%%%
    %%%%% Step 2 %%%%%
    %%%%%%%%%%%%%%%%%%

    eps_arr(j) = eps1 * kappa.^(((q + sigma)/p).^(j-1) - 1);

    %%%%%%%%%%%%%%%%%%
    %%%%% Step 3 %%%%%
    %%%%%%%%%%%%%%%%%%

    if(disp_flag >= 1)
        disp('                Applying Alg. 4.1...')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Alg. 4.1, Initialization %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    init_N_sample = local_options.init_N_sample;

    N_reused = 0;

    % Sample new elements of the eps-jet
    if(disp_flag >= 1)
        disp(['                    Initially sampling ',num2str(init_N_sample),' point(s).'])
    end

    W = x_arr(:,j);
    oracle_W = cell(q+1,1);
    oracle_W{1} = f_x;
    for order = 2:q+1
        oracle_W{order,1} = oracle{order}(W(:,1)); eval_counter(order) = eval_counter(order) + 1;
    end
    N_W = 1;

    while(1)

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Alg. 4.1, Step 2 %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%

        % Solve subproblem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if(N_W == 1)
            
            poly1 = vpa(0);
            for i = double(0:q)
                poly1 = poly1 + 1/factorial(vpa(i)) * oracle_W{i+1} * (z -  W(:,1))^i;
            end

            % Critical points
            r = vpa(root(diff(poly1,z)));
            r = r(imag(r) == 0);

            z_list = r;

            % Filter our infeasible
            z_list = z_list(z_list >= x_arr(:,j) - eps_arr(j) & z_list <= x_arr(:,j) + eps_arr(j));
            
            % Boundary points
            z_list = [z_list,x_arr(:,j) - eps_arr(j),x_arr(:,j) + eps_arr(j)];
            
            
            % Find minimum
            theta_list = vpa(zeros(1,numel(z_list)));
            for i = 1:numel(z_list)
                theta_list(i) = subs(poly1,z,z_list(i));
            end
            
            [~,min_i] = min(theta_list);
            z_bar = z_list(min_i);
            theta = theta_list(min_i);

        elseif(N_W == 2)

            poly1 = vpa(0);
            for i = double(0:q)
                poly1 = poly1 + 1/factorial(vpa(i)) * oracle_W{i+1,1} * (z - W(:,1))^i;
            end

            poly2 = vpa(0);
            for i = double(0:q)
                poly2 = poly2 + 1/factorial(vpa(i)) * oracle_W{i+1,2} * (z - W(:,2))^i;
            end

            % Critical points
            r1 = vpa(root(diff(poly1,z)));
            r1 = r1(imag(r1) == 0);

            r2 = vpa(root(diff(poly2,z)));
            r2 = r2(imag(r2) == 0);

            z_list = [r1,r2];

            % Nonsmooth points
            r_nonsm = vpa(root(poly2 - poly1,z));
            r_nonsm = r_nonsm(imag(r_nonsm) == 0);

            z_list = [z_list,r_nonsm'];
            
            % Filter our infeasible
            z_list = z_list(z_list >= x_arr(:,j) - eps_arr(j) & z_list <= x_arr(:,j) + eps_arr(j));

            % Boundary points
            z_list = [z_list,x_arr(:,j) - eps_arr(j),x_arr(:,j) + eps_arr(j)];
            

            % Find minimum
            theta_list = vpa(zeros(1,numel(z_list)));
            for i = 1:numel(z_list)
                poly1_z = subs(poly1,z,z_list(i));
                poly2_z = subs(poly2,z,z_list(i));
                theta_list(i) = max(poly1_z,poly2_z);
            end
            
            [~,min_i] = min(theta_list);
            z_bar = z_list(min_i);
            theta = theta_list(min_i);

        else
            input('N_W > 2 not supported.')
        end

        mu = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

        if(norm(z_bar - x_arr(:,j),norm_flag) > 1.5*eps_arr(j))
            fprintf(2,['Warning: eps-constraint is heavily violated! eps = ',num2str(eps_arr(j)), ', ||z-x|| = ',num2str(norm(z_bar - x_arr(:,j),norm_flag)),', frac = ',num2str(norm(z_bar - x_arr(:,j),norm_flag)/eps_arr(j)),'\n'])
        end

        f_z_bar = oracle{1}(z_bar); eval_counter(1) = eval_counter(1) + 1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Alg. 4.1, Step 3 %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % [f_z_bar - theta, eps_arr(j)^(q+sigma)]
        if(f_z_bar - theta <= eps_arr(j)^(q+sigma) )

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Alg. 4.1, Step 4 %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            break
        else

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Alg. 4.1, Step 6 %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            % Sample new information at z_bar
            N_W = N_W + 1;
            W = [W,z_bar];
            oracle_W{1,N_W} = f_z_bar;
            for order = 2:q+1
                oracle_W{order,N_W} = oracle{order}(z_bar); eval_counter(order) = eval_counter(order) + 1;
            end
        end
    end

    if(disp_flag >= 1)
        disp('                ...Alg. 4.1 finished!')
        disp(['                ',num2str(N_W - init_N_sample - N_reused),' additional sample points required.']);
        disp(['                f(z_bar) = ',sprintf('%9.2e', f_z_bar)]);
        disp(['                theta    = ',sprintf('%9.2e', theta)]);
        disp(['                f(z_bar) - theta  = ',sprintf('%9.2e', f_z_bar - theta)]);
        disp(['                eps^(q+sigma)     = ',sprintf('%9.2e', eps_arr(j)^(q+sigma))]);
        disp(' ');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% End of Alg. 4.1 %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%

    % Save number of oracle calls up to this iteration
    eval_counter_arr(j+1,:) = eval_counter;

    % Save number of samples required
    numsample_arr(j) = eval_counter_arr(j+1,2) - eval_counter_arr(j,2);

    % Save multiplier of trust-region constraint
    mu_arr(j) = mu;

    % Save activity of trust-region constraint
    act_arr(j) = norm(z_bar - x_arr(:,j),norm_flag)/eps_arr(j);

    if(disp_flag >= 1)
        disp(['                activity   = ',num2str(double(act_arr(j)))]);
        disp(['                act_thresh = ',num2str(double(act_thr))]);
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
    if(j > j_thr && act_arr(j) > act_thr)
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

    disp(['            f(x1)    = ',num2str(double(f_x1))]);
    disp(['            best_f   = ',num2str(double(best_f_val))]);
    disp(['            f(z_bar) = ',num2str(double(f_z_bar))]);

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

