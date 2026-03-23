% Script that performs the experiment in Example 6.5 and creates the plots
% shown in the article.

% Flag for saving the plots
save_figure_flag = 0;

%% Define the problem

addpath('max_root');
n = 100;

problem_data.n = n;
problem_data.oracle = {@(x) max_root(x); @(x) grad_max_root(x); @(x) hess_max_root(x)};

%% Set the parameters for the method
clear algo_options local_options global_options

% General parameters
algo_options.q = 1;
algo_options.p = 1;
algo_options.sigma = 0.5;
algo_options.disp_flag = 0;
algo_options.memory_max_size = 10000;

% Parameters for local phase
local_options.kappa = 0.75;
local_options.eps_thr = 10^-7;
local_options.j_thr = Inf; % Disables termination due to active trust-region constraint
local_options.act_thr = 0.95;
local_options.init_N_sample = 1;

local_options.norm_flag = Inf;
local_options.sp_solver = 'linprog';
local_options.sp_solver_optns = optimoptions("linprog","Display","none","OptimalityTolerance",10^-10);

algo_options.local_options = local_options;

%% Run the algorithm
addpath('../../Algorithms')

eps1 = 0.5;
x1 = 0.1 * (1:n)'/n;

reps = 1;
total_time = 0;
disp(['Running Alg. 4.2 ',num2str(reps),' times for the average...'])
for i = 1:reps
    start_tic = tic;
    result_local_method = local_method(x1,eps1,problem_data,algo_options);
    total_time = total_time + toc(start_tic);
end
runtime_alg4_2 = total_time/reps;
disp(['Average runtime = ',num2str(runtime_alg4_2)]);

x_arr = result_local_method.x_arr;
f_arr = result_local_method.f_arr;
act_arr = result_local_method.act_arr;
mu_arr = result_local_method.mu_arr;
numsample_arr = result_local_method.numsample_arr;
best_x = result_local_method.best_x;
best_f_val = result_local_method.best_f_val;
eval_counter = result_local_method.eval_counter;
eval_counter_arr = result_local_method.eval_counter_arr;
success_flag = result_local_method.success_flag;

%% Plots 

x_min = zeros(n,1);
f_min = 0;

load('SuperPolyak/SuperPolyak-result-16-12-25_12-48-59')

% Use vpa to handle small objective value from SuperPolyak
f_SupPol = double(max_root(vpa(solution_SupPol,32)));
loss_history_SupPol(end) = f_SupPol;

disp(['f(Alg. 4.2)    = ',num2str(result_local_method.best_f_val)])
disp(['f(SuperPolyak) = ',num2str(f_SupPol)])
disp(['VUbundle total gradient calls = ',num2str(sum(oracle_calls_SupPol))])
disp(['Alg. 4.2 total oracle calls   = ',num2str(eval_counter(1))])
disp(['Alg. 4.2 maximum oracle calls in Alg 4.1 = ',num2str(max(numsample_arr))])
disp(['Runtime Alg. 4.2 = ',num2str(runtime_alg4_2)]);

lw = 1.5;
ms = 13;

kappa = local_options.kappa;
eps_thr = local_options.eps_thr;
act_thr = local_options.act_thr;
sigma = algo_options.sigma;
q = algo_options.q;
p = algo_options.p;
size_x_arr = size(x_arr,2);
j_max = floor(1 + (log(log(eps_thr/eps1)/log(kappa) + 1))/log((q+sigma)/p)) + 1;
eps_arr = eps1 * kappa.^(((q + sigma)/p).^((1:j_max)-1) - 1);

fig1 = figure; % ----------------------------------------------------------

plot(1:size_x_arr,log10(vecnorm(x_arr - x_min,2,1)),'k.-','LineWidth',lw,'MarkerSize',ms);
hold on
plot(1:j_max,log10(eps_arr),'r.:','LineWidth',lw,'MarkerSize',ms);

% dummy plots for legend
h1 = plot(-10,-10,'k.-','LineWidth',lw,'MarkerSize',15);
h2 = plot(-10,-10,'r.:','LineWidth',lw,'MarkerSize',15);

legend([h1,h2],{'$\| x^j - x^* \|$','$\varepsilon_j$'},'Interpreter','latex','FontSize',18,'Location','sw');

xlabel('$j$','Interpreter','latex');
grid on

xlim([1,size_x_arr]);
ylim([-15.0,2.5])
set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

% Log tick labeling (y)
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)

if(save_figure_flag)
    export_fig 'plot1' '-png' '-r500' '-transparent'
    close(fig1)
end

fig2 = figure; % ----------------------------------------------------------

cumul_oracle_calls_SupPol = zeros(1,numel(oracle_calls_SupPol));
for j = 1:numel(oracle_calls_SupPol)
    cumul_oracle_calls_SupPol(j) = sum(oracle_calls_SupPol(1:j));
end

plot(eval_counter_arr(:,1),log10(f_arr - f_min),'k.-','LineWidth',lw,'MarkerSize',15);
hold on
plot(cumul_oracle_calls_SupPol,log10(loss_history_SupPol - f_min),'b.--','LineWidth',lw,'MarkerSize',15);

% dummy plots for legend
h1 = plot(-100,-100,'k.-','LineWidth',lw,'MarkerSize',15);
h2 = plot(-100,-100,'b.--','LineWidth',lw,'MarkerSize',15);

legend([h1,h2],{'Alg.\ 4.2','\verb+SuperPolyak+'},'Interpreter','latex','FontSize',18,'Location','ne');

xlabel('Oracle calls $l$','Interpreter','latex');
grid on

xlim([0,eval_counter_arr(end,1)])
ylim([-17,1])
set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

% Log tick labeling (y)
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)

if(save_figure_flag)
    export_fig 'plot2' '-png' '-r500' '-transparent'
    close(fig2)
end