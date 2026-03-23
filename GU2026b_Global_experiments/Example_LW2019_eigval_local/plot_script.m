% Script that performs the experiment in Example 5.4 and creates the plot
% shown in the article. 

save_figure_flag = 0;

%% Define the problem

% LW2019, (8.7)
addpath('../Example_LW2019_eigval_global/LW2019_eigval/');
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
problem_data.x0 = ones(n,1);

%% Set the parameters for the method
clear algo_options local_options global_options

% General parameters
algo_options.q = 2;
algo_options.p = 2;
algo_options.sigma = 0.5;
algo_options.disp_flag = 2;
algo_options.memory_max_size = 100;
algo_options.local_flag = 1;

% Parameters for global phase
global_options.delta_arr = 10.^(0:-1:-4);
global_options.tau_arr = 10^(-5)*ones(1,numel(global_options.delta_arr));
global_options.c = 0.1;
global_options.i_max = 10000;

global_options.init_N_sample = 1;

global_options.norm_flag = 2;
global_options.sp_solver = 'IPOPT';
global_options.sp_solver_optns.tol = 10^-10;

algo_options.global_options = global_options;

% Parameters for local phase
local_options.kappa = 0.75;
local_options.eps_thr = 10^-3;
local_options.j_thr = 2;
local_options.act_thr = 0.95;
local_options.init_N_sample = 1;

local_options.norm_flag = 2;
local_options.sp_solver = 'IPOPT';
local_options.sp_solver_optns.tol = 10^-10;

algo_options.local_options = local_options;

%% Run the algorithm
addpath('../../Algorithms')
[result_global_phase,result_local_phase] = global_method(problem_data,algo_options);

x_cell = result_global_phase.x_cell;
eval_counter = result_global_phase.eval_counter;
numsample_cell = result_global_phase.numsample_cell;
best_f_val = result_global_phase.best_f_val;
best_x = result_global_phase.best_x;
local_success_flag = result_global_phase.local_success_flag;

%% Plots

load('hanso_sol_10-Mar-2026_13_52_48')
x_min = x_hanso;

f_x_min = problem_data.oracle{1}(x_min);

x_arr_all = [x_cell{:}];
f_x_arr_all = zeros(1,size(x_arr_all,2));
for i = 1:size(x_arr_all,2)
    f_x_arr_all(i) = problem_data.oracle{1}(x_arr_all(:,i));
end

j_max = numel(x_cell);
ind_cell = cell(1,j_max);
ind_cell{1} = 1:size(x_cell{1},2);
for j = 2:j_max
    ind_cell{j} = (ind_cell{j-1}(end) + 1):(ind_cell{j-1}(end) + size(x_cell{j},2));
end

j_subinds = cellfun(@(in) in(:,end),ind_cell);

lw = 1.5;
ms = 15;

figure
h1 = plot(log10(vecnorm(x_arr_all - x_min,2,1)),'k.-','MarkerSize',ms,'LineWidth',lw);
hold on
h4 = plot(j_subinds,log10(vecnorm(x_arr_all(:,j_subinds) - x_min,2,1)),'ko','MarkerSize',9,'LineWidth',lw);
for j = 1:j_max
    h2 = plot([ind_cell{j}(1)-0.5,ind_cell{j}(end)+0.5],log10(algo_options.global_options.delta_arr(j))*[1,1],'r-','LineWidth',lw);

    if(j < j_max)
        xline(ind_cell{j}(end)+0.5,'k--','LineWidth',lw)
    end
end

if(local_success_flag)
    xline(ind_cell{end}(end)+0.5,'k--','LineWidth',lw)
    h3 = plot(size(x_arr_all,2) + (1:size(result_local_phase.x_arr,2)), log10(vecnorm(result_local_phase.x_arr - x_min,2,1)),'ks-','MarkerSize',8,'LineWidth',lw);
    plot(size(x_arr_all,2) + (1:size(result_local_phase.x_arr,2)-1), log10(result_local_phase.eps_arr), 'r:','LineWidth',lw)
end

legend([h1,h4,h2],{'$\| x^{j,i} - x^* \|$', '$\| x^j - x^* \|$', '$\Delta_j$'},'Interpreter','latex','FontSize',15,'Location','sw');

if(local_success_flag)
    legend([h1,h4,h2,h3],{'$\| x^{j,i} - x^* \|$', '$\| x^j - x^* \|$', '$\Delta_j$', 'Local method'},'Interpreter','latex','FontSize',23,'Location','sw');
end

grid on
axis square

set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

if(local_success_flag)
    xlim([0.5,size(x_arr_all,2) + size(result_local_phase.x_arr,2) + 0.5])
else
    xlim([0.5,ind_cell{end}(end)+0.5])
end


xticks(cellfun(@(in) in(1),ind_cell));
label_cell = cell(1,j_max);
for j = 1:j_max
    label_cell{j} = ['(',num2str(j),',0)'];
end
xticklabels(label_cell);

if(local_success_flag)
    xticks([xticks,size(x_arr_all,2)+1])
    xticklabels([label_cell, {['(',num2str(j_max+1),',0)']}])
end

xtickangle(60);
ax = gca;
ax.XAxis.FontSize = 13;

ylim([-5.5,1.1])

% Log tick labeling (y)
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)


%%
if(save_figure_flag)
    export_fig 'plot' '-png' '-r500' '-transparent'
end