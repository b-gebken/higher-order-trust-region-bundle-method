% A script for visualizing the example in Remark 5.1 symbolically

n = 1;
f = @(x) x^2;
x0 = sym(0.5);

j_max = 20;

j_arr = 1:j_max;
delta_arr = sym( (1/2).^(j_arr.^2) - (1/2).^((j_arr + 1).^2) );
tau_arr = sym(2 * (1/2).^(2*(j_arr + 1).^2) ./ delta_arr.^2 );

%%

i_max = 5;

x_cell = cell(1,j_max);
x_arr = sym(zeros(1,i_max));
x_arr(1) = x0;

quot_cell = cell(1,j_max);
quot_arr = sym(zeros(1,i_max));

ji_cell = cell(1,j_max);
ji_arr = zeros(1,i_max);

for j = 1:j_max
    for i = 1:i_max
        ji_arr(i) = j;

        if(abs(x_arr(i)) > delta_arr(j))
            z_bar =  x_arr(i) - sign(x_arr(i))*delta_arr(j);
        else
            z_bar = sym( 0 );
        end

        quot_arr(i) = (f(x_arr(i)) - f(z_bar)) / delta_arr(j)^2;
        if( (f(x_arr(i)) - f(z_bar)) / delta_arr(j)^2  <= tau_arr(j))
            x_arr = x_arr(1:i);
            quot_arr = quot_arr(1:i);
            ji_arr = ji_arr(1:i);
            break
        else
            x_arr(i+1) = z_bar;
        end
    end

    x_cell{j} = x_arr(1:i);
    x_arr = sym(zeros(1,i_max));
    x_arr(1) = x_cell{j}(end);

    quot_cell{j} = quot_arr;
    quot_arr = sym(zeros(1,i_max));

    ji_cell{j} = ji_arr;
    ji_arr = zeros(1,i_max);
end

%% 

problem_data.oracle = {@(x) x^2};

x_min = zeros(n,1);
f_x_min = problem_data.oracle{1}(x_min);

x_arr_all = double([x_cell{:}]);
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
    h2 = plot([ind_cell{j}(1)-0.5,ind_cell{j}(end)+0.5],log10(delta_arr(j))*[1,1],'r-','LineWidth',lw);

    if(j < j_max)
        xline(ind_cell{j}(end)+0.5,'k--','LineWidth',lw)
    end
end

legend([h1,h4,h2],{'$\| x^{j,i} - x^* \|$', '$\| x^j - x^* \|$', '$\Delta_j$'},'Interpreter','latex','FontSize',15,'Location','sw');

grid on
axis square

set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

xlim([0.5,ind_cell{end}(end)+0.5])

xticks(cellfun(@(in) in(1),ind_cell));
label_cell = cell(1,j_max);
for j = 1:j_max
    label_cell{j} = ['(',num2str(j),',0)'];
end
xticklabels(label_cell);

xtickangle(60);
ax = gca;
ax.XAxis.FontSize = 13;

% ylim([-7,1.5])

% Log tick labeling (y)
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)

%%

x_arr_all = [x_cell{:}];

lw = 1.5;
ms = 15;

figure

quot_arr_all = [quot_cell{:}];
ji_arr_all = [ji_cell{:}];

h1 = plot(double(vpa(log10(quot_arr_all))),'k.-','MarkerSize',ms,'LineWidth',lw);
hold on
h2 = plot(double(vpa(log10(tau_arr(ji_arr_all)))),'r.-','MarkerSize',ms,'LineWidth',lw);

legend([h2,h1],{'$\tau_j$', '$(f(x^{j,i}) - f(\bar{z}^{j,i})) / \Delta_j^2$'},'Interpreter','latex','FontSize',15,'Location','sw');

grid on

set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

title('Stopping criterion')

xlim([0.5,ind_cell{end}(end)+0.5])

xticks(cellfun(@(in) in(1),ind_cell));
label_cell = cell(1,j_max);
for j = 1:j_max
    label_cell{j} = ['(',num2str(j),',0)'];
end
xticklabels(label_cell);

xtickangle(60);
ax = gca;
ax.XAxis.FontSize = 13;

% ylim([-7,1.5])

% Log tick labeling (y)
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)