% A script for visualizing the example in Remark 5.1 using Alg. 2. (Due to
% the fast decrease of the involved sequences, only few iterations can be
% executed.) 

save_figure_flag = 0;

%% Define the problem

n = 1;

problem_data.n = n;
problem_data.oracle = {@(x) x^2; @(x) 2*x; @(x) 2};
problem_data.x0 = 0.5;

% Paste this into generate_W for exact subproblem solves
% ---------------------------
% if(abs(x) <= eps)
%     z_bar = 0;
% else
%     z_bar = x - sign(x)*eps;
% end
% theta = z_bar^2;
% mu = 0;
% ---------------------------

%% Set the parameters for the method
clear algo_options local_options global_options

% General parameters
algo_options.q = 2;
algo_options.p = 2;
algo_options.sigma = 0.5;
algo_options.disp_flag = 2;
algo_options.memory_max_size = 0;
algo_options.local_flag = 0;

% Parameters for global phase
e = 2.0;
a = @(j) (1/2).^(j.^e);
j_max = 5;
global_options.delta_arr = a(1:(j_max-1)) - a(2:j_max);
global_options.tau_arr = 2*( 1./( (1/2).^((1:j_max).^e - ((1:j_max)+1).^e) - 1 ) ).^2;
global_options.c = 0.1;
global_options.i_max = 10000;

global_options.init_N_sample = 1;

global_options.norm_flag = 2;
global_options.sp_solver = 'fmincon';
global_options.sp_solver_optns = optimoptions("fmincon","Display","none","Algorithm","sqp","MaxFunctionEvaluations",10^6);

algo_options.global_options = global_options;


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

x_min = zeros(n,1);
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
    h3 = plot(size(x_arr_all,2) + (1:size(result_local_phase.x_arr,2)), log10(vecnorm(result_local_phase.x_arr - x_min,2,1)),'ko-','MarkerSize',6,'LineWidth',lw);
    plot(size(x_arr_all,2) + (1:size(result_local_phase.x_arr,2)-1), log10(result_local_phase.eps_arr), 'r:','LineWidth',lw)
end

legend([h1,h4,h2],{'$\| x^{j,i} - x^* \|$', '$\| x^j - x^* \|$', '$\Delta_j$'},'Interpreter','latex','FontSize',15,'Location','sw');

if(local_success_flag)
    legend([h1,h4,h2,h3],{'$\| x^{j,i} - x^* \|$', '$\| x^j - x^* \|$', '$\Delta_j$', 'Local method'},'Interpreter','latex','FontSize',15,'Location','sw');
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