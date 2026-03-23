% Script that performs the experiment in Example 6.2 and creates the plots
% shown in the article.

% Flag for saving the plots
save_figure_flag = 0;

%% Define the problem

addpath('LW2019_85');
n = 50;
k = 40;

rng("default");

lambda = rand(k,1); lambda = lambda/sum(lambda);
tmp = 2*rand(n,k)-1; g_arr = tmp - tmp*lambda;

H_cell = cell(1,k);
for i = 1:k
    tmp = 2*rand(n)-1;
    H_cell{i} = tmp * tmp';
end

c_arr = rand(1,k);

problem_data.n = n;
problem_data.oracle = {@(x) lw2019_85_f(x,g_arr,H_cell,c_arr); @(x) lw2019_85_grad(x,g_arr,H_cell,c_arr); @(x) lw2019_85_hess(x,g_arr,H_cell,c_arr)};

%% Set the parameters for the method
clear algo_options local_options global_options

% General parameters
algo_options.q = 2;
algo_options.p = 2;
algo_options.sigma = 0.5;
algo_options.disp_flag = 2;
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

eps1 = 10;
x1 = ones(n,1);
result_local_method = local_method(x1,eps1,problem_data,algo_options);

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
ylim([-8.5,1.5])
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

plot(1:size_x_arr-1,numsample_arr,'k.-','LineWidth',lw,'MarkerSize',ms);
hold on

% dummy plots for legend
h1 = plot(-10,-10,'k.-','LineWidth',lw,'MarkerSize',15);

legend(h1,{'Oracle calls'},'Interpreter','latex','FontSize',18,'Location','nw');

xlabel('$j$','Interpreter','latex');
grid on

xlim([1,size_x_arr-1]);
ylim([0,18]);
set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

if(save_figure_flag)
    export_fig 'plot2' '-png' '-r500' '-transparent'
    close(fig2)
end

fig3 = figure; % ----------------------------------------------------------

plot(eval_counter_arr(:,1),log10(vecnorm(x_arr - x_min,2,1)),'k.-','LineWidth',lw,'MarkerSize',ms);
hold on

% dummy plots for legend
h1 = plot(-10,-10,'k.-','LineWidth',lw,'MarkerSize',15);

legend(h1,{'$\| x^{j(l)} - x^* \|$'},'Interpreter','latex','FontSize',18,'Location','ne');

xlabel('Oracle calls $l$','Interpreter','latex');
grid on

xlim([0,eval_counter_arr(end,1)]);
ylim([-8.5,1.5])
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
    export_fig 'plot3' '-png' '-r500' '-transparent'
    close(fig3)
end