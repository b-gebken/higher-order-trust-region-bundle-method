% Script that performs the experiment in Example 6.3 and creates the plots
% shown in the article.

% Flag for saving the plots
save_figure_flag = 0;

%% Define the problem

addpath('LW2019_eigval');
n = 50;
m = 25;

rng("default")

A_tensor = zeros(n+1,m,m);
for i = 1:n+1
    tmp = 2*rand(m) - 1;
    A_tensor(i,:,:) = (1/2)*(tmp + tmp');
end

problem_data.n = n;
problem_data.oracle = {@(x) lw2019_eigval_f(x,A_tensor); @(x) lw2019_eigval_grad(x,A_tensor); @(x) lw2019_eigval_hess(x,A_tensor)};

%% Set the parameters for the method
clear algo_options local_options global_options

% General parameters
algo_options.q = 2;
algo_options.p = 2;
algo_options.sigma = 0.5;
algo_options.disp_flag = 0;
algo_options.memory_max_size = 10000;

% Parameters for local phase
local_options.kappa = 0.75;
local_options.eps_thr = 10^-3;
local_options.j_thr = Inf; % Disables termination due to active trust-region constraint
local_options.act_thr = 0.95;
local_options.init_N_sample = 1;

local_options.norm_flag = 2;
local_options.sp_solver = 'IPOPT';
local_options.sp_solver_optns.tol = 10^-10;

% Code for using fmincon if desired (less accurate and slower than IPOPT in
% our experience)
% local_options.norm_flag = 2;
% local_options.sp_solver = 'fmincon';
% local_options.sp_solver_optns = optimoptions("fmincon","Display","none","Algorithm","interior-point","MaxFunctionEvaluations",10^6);

algo_options.local_options = local_options;

%% Run the algorithm
addpath('../../Algorithms')

eps1 = 0.5;
x1 = zeros(n,1);

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

load('hanso_sol_10-Mar-2026_13_52_48')

disp(['f(Alg. 4.2) - f(HANSO) = ',num2str(result_local_method.best_f_val - f_hanso)])
disp(['Runtime Alg. 4.2 = ',num2str(runtime_alg4_2)]);
disp(['Runtime HANSO    = ',num2str(runtime_hanso)]);
disp(['HANSO total oracle calls    = ',num2str(cumul_oracle_calls_hanso(end))])
disp(['Alg. 4.2 total oracle calls = ',num2str(eval_counter(1))])
disp(['Alg. 4.2 maximum oracle calls in Alg 4.1 = ',num2str(max(numsample_arr))])

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

plot(1:size_x_arr,log10(vecnorm(x_arr - x_hanso,2,1)),'k.-','LineWidth',lw,'MarkerSize',ms);
hold on
plot(1:j_max,log10(eps_arr),'r.:','LineWidth',lw,'MarkerSize',ms);

% dummy plots for legend
h1 = plot(-10,-10,'k.-','LineWidth',lw,'MarkerSize',15);
h2 = plot(-10,-10,'r.:','LineWidth',lw,'MarkerSize',15);

legend([h1,h2],{'$\| x^j - \tilde{x}^* \|$','$\varepsilon_j$'},'Interpreter','latex','FontSize',18,'Location','sw');

xlabel('$j$','Interpreter','latex');
grid on

xlim([1,size_x_arr]);
ylim([-5.5,0.5])
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

fevalrec_arr_hanso = [fevalrec_hanso{:}];

disp(['Oracle calls for HANSO to reach accuracy of Alg. 4.2 = ',num2str(find(fevalrec_arr_hanso < result_local_method.best_f_val,1,'first'))])

plot(eval_counter_arr(:,1),log10(f_arr - f_hanso),'k.-','LineWidth',lw,'MarkerSize',ms);
hold on
plot(1:numel(fevalrec_arr_hanso),log10(fevalrec_arr_hanso - f_hanso),'b--','LineWidth',lw,'MarkerSize',ms);

% dummy plots for legend
h1 = plot(-100,-100,'k.-','LineWidth',lw,'MarkerSize',15);
h2 = plot(-100,-100,'b--','LineWidth',lw,'MarkerSize',15);

legend([h1,h2],{'Alg.\ 4.2','\verb+HANSO+'},'Interpreter','latex','FontSize',18,'Location','ne');

xlabel('Oracle calls $l$','Interpreter','latex');
grid on

xlim([0,570]);
ylim([-7.0,1.5])
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