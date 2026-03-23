% A script that applies the local method ([GU2026a], Alg. 4.2) to a simple
% example. 

%% Define the problem
% In this case half-and-half from Section 5.5 in [Lewis, Overton (2008),
% "Nonsmooth optimization via BFGS"].

n = 8;
A = diag([1,0,1,0,1,0,1,0]);
B = diag(1./(1:8).^2);

f = @(x) sqrt(x'*A*x) + x'*B*x;
grad_f = @(x) (x'*A*x)^(-1/2)*A*x + 2*B*x;
hess_f = @(x) (-(x'*A*x)^(-3/2)*A*x)*(A*x)' + (x'*A*x)^(-1/2)*A + 2*B;

problem_data.n = n;
problem_data.oracle = {f; grad_f; hess_f};

%% Set the parameters for the method
clear algo_options local_options global_options

% General parameters
algo_options.q = 2;
algo_options.p = 2;
algo_options.sigma = 0.5;
algo_options.disp_flag = 1;
algo_options.memory_max_size = 10000;

% Parameters for local phase
local_options.kappa = 0.75;
local_options.eps_thr = 10^-3;
local_options.j_thr = Inf;
local_options.act_thr = 0.95;
local_options.init_N_sample = 1;

local_options.norm_flag = 2;
local_options.sp_solver = 'IPOPT';
local_options.sp_solver_optns.tol = 10^-10;

% local_options.norm_flag = 2;
% local_options.sp_solver = 'fmincon';
% local_options.sp_solver_optns = optimoptions("fmincon","Display","none","Algorithm","interior-point","MaxFunctionEvaluations",10^6);

algo_options.local_options = local_options;

%% Run the algorithm

eps1 = 30;
x1 = 20.08 * ones(n,1);

result_local_method = local_method(x1,eps1,problem_data,algo_options);

x_arr = result_local_method.x_arr;
f_arr = result_local_method.f_arr;
eps_arr = result_local_method.eps_arr;
act_arr = result_local_method.act_arr;
mu_arr = result_local_method.mu_arr;
numsample_arr = result_local_method.numsample_arr;
best_x = result_local_method.best_x;
best_f_val = result_local_method.best_f_val;
eval_counter = result_local_method.eval_counter;
eval_counter_arr = result_local_method.eval_counter_arr;
success_flag = result_local_method.success_flag;


%% Plots

x_min = zeros(8,1);
f_min = 0;

lw = 1.5;
ms = 13;

size_x_arr = size(x_arr,2);
j_max = numel(eps_arr);

figure

subplot(1,2,1) % ----------------------------------------------------------
plot(1:size_x_arr,log10(vecnorm(x_arr - x_min,2,1)),'k.-','LineWidth',lw,'MarkerSize',ms);
hold on
h = plot(1:j_max,log10(eps_arr),'r.:','LineWidth',lw,'MarkerSize',ms);

legend(h,{'$\varepsilon_j$'},'Interpreter','latex','FontSize',18,'Location','sw');

title('$\| x^j - x^* \|$','Interpreter','latex');

xlabel('$j$','Interpreter','latex');
grid on

xlim([1,size_x_arr]);
ylim([-11.5,2.5])
set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

% Log tick labeling (y)
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)

subplot(1,2,2) % ----------------------------------------------------------

plot(eval_counter_arr(:,1),log10(f_arr - f_min),'k.-','LineWidth',lw,'MarkerSize',ms);

title('$f(x^{j(l)}) - f(x^*)$','Interpreter','latex')

xlabel('Oracle calls $l$','Interpreter','latex');
grid on

ylim([-11.5,3.5])
set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

% Log tick labeling (y)
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)