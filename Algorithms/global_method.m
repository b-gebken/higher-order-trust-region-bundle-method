%{ 
An implementation of Alg. 2 from [GU2026b]. Models of order q > 2 are
currently not supported. 

Inputs:
    problem_data: A struct with fields:
        n: Number of variables
        oracle: A cell array in which the i-th cell contains a function
            handle to evaluate the (i-1)-th derivative. For example,
            oracle{1}(x) is the objective value and oracle{2}(x) is the
            gradient at a point x.
        x0: Initial point
    algo_options: A struct with fields:
        q: Order of Taylor expansion to be used
        p: Order of growth of the objective function around the minimum
        disp_flag: Flag for displaying information of iterations
        memory_max_size: Maximum size of memory of past iterates for the
            bundle initialization in Alg. 4.1  (cf. [GU2026a], Remark 4.2(a)).
        local_flag: Flag for whether the local superlinear method from
            [GU2026a] should be applied after every outer j-iteration
        global_options: A struct with fields:
            delta_arr: Sequence of trust-region radii
            tau_arr: Sequence of threshold values deciding when to change
                the trust-region radius (Step 4 in [GU2026b])
            c: Threshold value for Alg. 1 in [GU2026a]
            i_max: Maximum number of inner i-iterations
            init_N_sample: Number of randomly sampled points for the
                initialization of Alg. 1. The first point is always the
                center of the trust-region, so init_N_sample = 1 means that
                no additional points are randomly sampled. (Always set to 1
                in our experiments.) 
            norm_flag: Which p-norm to use for the trust-region. (Not to be
                confused with the growth order p.)
            sp_solver: Which solver to use for trust-region subproblem.
                Currently supports "IPOPT" and "fmincon" (for q = 2 and
                norm_flag = 2) and "linprog" (for q = 1 and norm_flag =
                inf). ("IPOPT" uses a slightly modified version of
                mexIPOPT, see the comments in "solve_subproblem_IPOPT.m"
                for details.)   
            sp_solver_optns: Can be used to set solver specific options
        local_options: See "local_method.m"

Output:
    result_global_method: A struct with fields:
        local_success_flag: Whether the trust-region constraint in the
            local method was active from iteration j_thr onwards (see
            [GU2026a], Subsection 5.1)
        best_f_val: Smallest objective value found in global and local
            phases
        best_x: Point corresponding to best_f_val
        x_cell: Each cell contains the sequence generated for a fixed
            trust-region radius. In paper notation: 
            x^{j,i} = x_cell{j}(:,i+1)
        eval_counter: Number of oracle calls for each order
        numsample_cell: The number of oracle calls in Alg. 1, indexed as in
            x_cell
%} 

function [result_global_phase,result_local_phase] = global_method(problem_data,algo_options)

% Read problem specification
n = problem_data.n;
x0 = problem_data.x0;
oracle = problem_data.oracle;

% Read algorithm options
q = algo_options.q;
p = algo_options.p;
disp_flag = algo_options.disp_flag;
memory_max_size = algo_options.memory_max_size;
local_flag = algo_options.local_flag;

global_options = algo_options.global_options;
delta_arr = global_options.delta_arr;
tau_arr = global_options.tau_arr;
c = global_options.c;
norm_flag = global_options.norm_flag;
i_max = global_options.i_max;

% Tic for measuring runtime
if(disp_flag >= 1)
    start_tic = tic;
end

% Initialize variables
j_max = numel(delta_arr);
x_cell = cell(1,j_max);
numsample_cell = cell(1,j_max);
numsample_arr = zeros(1,i_max);
eval_counter = zeros(1,q+1);
result_local_phase(j_max) = struct('x_arr',[],'f_arr',[],'eps_arr',[],'act_arr',[],'mu_arr',[],'numsample_arr',[],'best_x',[],'best_f_val',[],'eval_counter',[],'eval_counter_arr',[],'success_flag',[]);

% Initialize the memory
memory.sample_pts = [];
memory.oracle_vals = {};
memory.max_size = memory_max_size;

x_arr = zeros(n,i_max+1);
x_arr(:,1) = x0;
f_x = oracle{1}(x0); eval_counter(1) = eval_counter(1) + 1;

% Tolerance for determining whether memorized points lie in the current
% trust region.
reusing_eps_tolerance = 10^-7;

%%%%%%%%%%%%%%%%%%
%%%%% Step 1 %%%%%
%%%%%%%%%%%%%%%%%%

% Outer j-loop
for j = 1:j_max

    if(disp_flag >= 2)
        disp(['----- Outer iteration j = ',num2str(j),' -------------------']);
        disp(['    delta_j = ',num2str(delta_arr(j)),', tau_j = ',num2str(tau_arr(j))]);
    end

    %%%%%%%%%%%%%%%%%%
    %%%%% Step 2 %%%%%
    %%%%%%%%%%%%%%%%%%

    % Inner i-loop
    for i = 1:i_max
        
        if(disp_flag >= 2)
            disp(['    ----- Inner iteration j = ',num2str(j),', i = ',num2str(i),' -----']);
            disp(['        delta_j = ',num2str(delta_arr(j)),', tau_j = ',num2str(tau_arr(j))]);
            disp(['        Applying Alg. 4.1...'])
        end

        %%%%%%%%%%%%%%%%%%
        %%%%% Step 3 %%%%%
        %%%%%%%%%%%%%%%%%%

        [z_bar,f_z_bar,mu,eval_counter,memory,numsample] = generate_W(x_arr(:,i),delta_arr(j),f_x,c,reusing_eps_tolerance,memory,eval_counter,problem_data,algo_options,global_options);

        % Save number of sample points
        numsample_arr(i) = numsample;

        if(disp_flag >= 2)
            disp(['        f(x) = ',sprintf('%9.2e', f_x)]);
            disp(['        f(x) - f(z) = ',sprintf('%9.2e', f_x - f_z_bar )]);
            disp(['        (f(x) - f(z))/delta_j^p = ',sprintf('%9.2e', (f_x - f_z_bar)/(delta_arr(j)^p) )]);
            disp(['        tau_j                 = ',sprintf('%9.2e', tau_arr(j))]);
            disp(['        ||z - x||/delta_j = ',sprintf('%9.2e', norm(z_bar - x_arr(:,i),norm_flag)/delta_arr(j))]);
        end

        %%%%%%%%%%%%%%%%
        %%%% Step 4 %%%%
        %%%%%%%%%%%%%%%%

        if((f_x - f_z_bar)/(delta_arr(j)^p) < tau_arr(j))
            
            %%%%%%%%%%%%%%%%
            %%%% Step 5 %%%%
            %%%%%%%%%%%%%%%%

            if(disp_flag >= 2)
                disp('        Decision: Break i-loop.')
                disp(' ');
            end

            x_arr = x_arr(:,1:i);
            numsample_arr = numsample_arr(1:i);
            break
        else
            %%%%%%%%%%%%%%%%
            %%%% Step 7 %%%%
            %%%%%%%%%%%%%%%%

            if(disp_flag >= 2)
                disp('        Decision: Sufficient decrease.')
                disp(' ');
            end

            x_arr(:,i+1) = z_bar;
            f_x = f_z_bar;
        end 
    end

    % Save information of inner iteration to cell arrays
    x_cell{j} = x_arr;
    numsample_cell{j} = numsample_arr;

    %%%%%%%%%%%%%%%%%
    %%%% Step 10 %%%%
    %%%%%%%%%%%%%%%%%

    if(local_flag)
        if(disp_flag >= 2)
            disp('        Attempting local method...')
        end
        result_local_method = local_method(x_arr(:,i),delta_arr(j),problem_data,algo_options);

        result_local_phase(j) = result_local_method;
        eval_counter = eval_counter + result_local_phase(j).eval_counter;
        local_success_flag = result_local_phase(j).success_flag;

        if(local_success_flag)
            result_local_phase = result_local_phase(1:j);
            x_cell = x_cell(1:j);
            numsample_cell = numsample_cell(1:j);
            break
        end
    else
        if(disp_flag >= 2)
            disp('        Skipping local method.')
        end
        local_success_flag = false;
    end

    %%%%%%%%%%%%%%%%%
    %%%% Step 11 %%%%
    %%%%%%%%%%%%%%%%%

    x_arr = zeros(n,i_max+1);
    x_arr(:,1) = x_cell{j}(:,end);

    if(i == i_max)
        fprintf(2,'Warning: Maximal number of iterations reached before stopping criterion is satisfied.\n')
    end
end

% Prepare output
result_global_phase.local_success_flag = local_success_flag;
if(local_success_flag)
    result_global_phase.best_f_val = result_local_phase(end).best_f_val;
    result_global_phase.best_x = result_local_phase(end).best_x;
else
    result_global_phase.best_f_val = f_x;
    result_global_phase.best_x = x_cell{j}(:,end);
end

result_global_phase.x_cell = x_cell;
result_global_phase.eval_counter = eval_counter;
result_global_phase.numsample_cell = numsample_cell;

if(disp_flag >= 1)
    disp(' ')
    if(local_success_flag)
        disp('Algorithm finished in local phase');
        disp(['    Final objective value =',sprintf('%9.2e', result_local_phase(end).best_f_val)]);
    else
        disp('Algorithm finished outside of local phase');
        disp(['    Final objective value =',sprintf('%9.2e', result_global_phase.best_f_val)]);
    end
    disp(['    Runtime               = ',num2str(toc(start_tic)),'s']);
    disp('    Required evaluations: ');
    disp(['        f       - ',num2str(eval_counter(1))])
    disp(['        subgrad - ',num2str(eval_counter(2))])
    if(q > 1)
        disp(['        subhess - ',num2str(eval_counter(3))])
    end
end

end

