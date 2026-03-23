% Script that performs the experiment in Example 6.1 and creates the plots
% shown in the article. It uses a modified version of local_method.m (Alg.
% 4.2) that works with symbolic arrays.

% Flag for saving the plots
save_figure_flag = 0;

%% Define the problem

% Accuracy for vpa
digits(2000)

n = 1;

% Maximum model order q
max_q = 5;

oracle_cell = cell(1,max_q+1);
oracle_cell{1} = @(x) vpa((abs(x) + vpa(1/4))^(1/2) - vpa(1/2));
for m = 1:max_q
    oracle_cell{m+1} = @(x) vpa((prod(1:-2:(1-2*(m-1))))/(2^m) * (abs(x) + vpa(1/4))^((1-2*m)/2) * sign(x)^m);
end

problem_data.n = n;

%% Set the parameters for the method
clear algo_options local_options global_options

% General parameters
algo_options.p = vpa(1);
algo_options.sigma = vpa(1/2);
algo_options.disp_flag = 2;
algo_options.memory_max_size = 0;

% Parameters for local phase
local_options.kappa = vpa(0.75);
local_options.eps_thr = vpa(10)^-500;
local_options.j_thr = Inf; % Disables termination due to active trust-region constraint
local_options.act_thr = vpa(0.95);
local_options.init_N_sample = 1;

local_options.norm_flag = 2;

algo_options.local_options = local_options;

%% Run the algorithm for q = 1:max_q

eps1 = vpa(0.5);
x1 = vpa(0.1);

result_cell = cell(1,max_q);

for q = 1:max_q
    
    problem_data.oracle = oracle_cell(1:q+1);
    algo_options.q = vpa(q);

    result_cell{q} = local_method_1d_sym(x1,eps1,problem_data,algo_options);

end

%% Plots

x_min = zeros(n,1);
f_min = 0;

kappa = local_options.kappa;
eps_thr = local_options.eps_thr;
sigma = algo_options.sigma;
p = algo_options.p;

lw = 1.5;
ms = 8;

marker_cell = {'d','*','o','x','s'};
handle_arr = zeros(1,6);

figure
hold on
for q = 1:max_q
    x_arr = result_cell{q}.x_arr;
    
    j_max = floor(1 + (log(log(eps_thr/eps1)/log(kappa) + 1))/log((q+sigma)/p)) + 1;
    eps_arr = eps1 * kappa.^(((q + sigma)/p).^((1:j_max)-1) - 1);

    handle_arr(q) = plot(1:size(x_arr,2),double(log10(abs(x_arr - x_min))),'k-','LineWidth',lw,'MarkerSize',ms,'Marker',marker_cell{q});
    plot(1:numel(eps_arr),log10(eps_arr),'r:','LineWidth',lw,'MarkerSize',ms,'Marker',marker_cell{q})
end

% dummy plot for legend
handle_arr(6) = plot(1:numel(eps_arr),log10(eps_arr),'r:','LineWidth',lw);

legend(handle_arr,{'$q = 1$','$q = 2$','$q = 3$','$q = 4$','$q = 5$','$\varepsilon_j$'},'Position',[0.6192    0.1763    0.1903    0.4269],'Interpreter','latex','FontSize',18)

xlabel('$j$','Interpreter','latex');
grid on
box on

xlim([1,22]);
set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

% Log tick labeling (y)
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)

% Check |W| = 2
for q = 1:max_q
    if(all(result_cell{q}.numsample_arr == 2))
        disp(['|W| = 2 for q = ',num2str(q)]);
    end
end

if(save_figure_flag)
    export_fig 'plot' '-png' '-r500' '-transparent'
end

