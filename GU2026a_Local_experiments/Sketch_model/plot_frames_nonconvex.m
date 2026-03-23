% Script that creates Fig. 1

% Flag for saving the plots
save_figure_flag = 0;

%% Define the function and compute its derivatives

x_sym = sym('x',[1,1]);

f1_sym = -(x_sym+0.5)^2 +0.25*abs(x_sym)^(3/2) + 0.5;  
f2_sym = x_sym^2 + 0.5*abs(x_sym)^(3/2) - 0.25;
f3_sym = -1/(abs(x_sym) + 0.25) + 2;
f_sym_tmp = (1/2) * (f1_sym + f2_sym + abs(f1_sym - f2_sym));
f_sym = (1/2) * (f_sym_tmp + f3_sym + abs(f_sym_tmp - f3_sym));

W = [-1.2,-0.9,-0.3,0.75,1.25];
num_W = size(W,2);

q_max = 5;

disp('Computing derivatives...')
tic;
deriv_sym_arr = sym(zeros(q_max+1,1));
deriv_sym_arr(1) = subs(f_sym,conj(x_sym),x_sym);
for i = 1:q_max
    deriv_sym_arr(i+1) = subs(diff(deriv_sym_arr(i),x_sym,1),conj(x_sym),x_sym);
end
deriv_mat = double(vpa(subs(deriv_sym_arr,x_sym,W),10));
t = toc;
disp(['... done in ',num2str(t),'s.'])

%% Plots
f = matlabFunction(f_sym,'Vars',{x_sym(:)});

plot_N = 5000;
plot_min = -1.3;
plot_max = 1.3;
plot_arr = linspace(plot_min,plot_max,plot_N);

lw = 2.4;
ms = 27;
text_size = 25;

% Create plots for q in {1,2,5}
for plot_i = 1:3

    fig = figure;
    hold on;

    switch plot_i

        case 1
            h0 = plot(plot_arr,f(plot_arr),'k-','LineWidth',lw);
            plot(W,f(W),'k.','MarkerSize',ms);

            q = 1;
            model_mat = NaN(num_W,plot_N);
            set(gca,'ColorOrderIndex',1)
            for i = 1:num_W
                model_mat(i,:) = eval_Taylor(plot_arr,deriv_mat(1:q+1,i),W(:,i));
                set(gca,'ColorOrderIndex',i)
                tmp_h = plot(plot_arr,model_mat(i,:),'--','LineWidth',lw);
                h2 = plot(W(i),f(W(i)),'.','MarkerSize',ms,'Color',tmp_h.Color);
            end

            text(0,2.25,['$q = ',num2str(q),'$'],'FontSize',text_size,'Interpreter','latex','VerticalAlignment','middle','HorizontalAlignment','center');
            h3 = plot(plot_arr,max(model_mat,[],1),'r-','LineWidth',lw);

            %dummy plot for legend
            h1 = plot(-100,-100,'k--','LineWidth',lw);

            legend([h0,h1,h3],{'$f$','$T^q f$','Model'},'Location','sw','FontSize',18,'Interpreter','latex');

            [max_error_1,I_1] = max(abs(max(model_mat,[],1) - f(plot_arr)));
            disp(['max |R(z)| = ',num2str(max_error_1),' at z = ',num2str(plot_arr(I_1))])

        case 2
            plot(plot_arr,f(plot_arr),'k-','LineWidth',lw);
            plot(W,f(W),'k.','MarkerSize',ms);

            q = 2;
            model_mat = NaN(num_W,plot_N);
            set(gca,'ColorOrderIndex',1)
            for i = 1:num_W
                model_mat(i,:) = eval_Taylor(plot_arr,deriv_mat(1:q+1,i),W(:,i));
                set(gca,'ColorOrderIndex',i)
                tmp_h = plot(plot_arr,model_mat(i,:),'--','LineWidth',lw);
                plot(W(i),f(W(i)),'.','MarkerSize',ms,'Color',tmp_h.Color);
            end

            text(0,2.25,['$q = ',num2str(q),'$'],'FontSize',text_size,'Interpreter','latex','VerticalAlignment','middle','HorizontalAlignment','center');
            h = plot(plot_arr,max(model_mat,[],1),'r-','LineWidth',lw);

            [max_error_2,I_2] = max(abs(max(model_mat,[],1) - f(plot_arr)));
            disp(['max |R(z)| = ',num2str(max_error_2),' at z = ',num2str(plot_arr(I_2))])

        case 3
            plot(plot_arr,f(plot_arr),'k-','LineWidth',lw);
            plot(W,f(W),'k.','MarkerSize',ms);

            q = 5;
            model_mat = NaN(num_W,plot_N);
            set(gca,'ColorOrderIndex',1)
            for i = 1:num_W
                model_mat(i,:) = eval_Taylor(plot_arr,deriv_mat(1:q+1,i),W(:,i));
                set(gca,'ColorOrderIndex',i)
                tmp_h = plot(plot_arr,model_mat(i,:),'--','LineWidth',lw);
                plot(W(i),f(W(i)),'.','MarkerSize',ms,'Color',tmp_h.Color);
            end

            text(0,2.25,['$q = ',num2str(q),'$'],'FontSize',text_size,'Interpreter','latex','VerticalAlignment','middle','HorizontalAlignment','center');
            h = plot(plot_arr,max(model_mat,[],1),'r-','LineWidth',lw);

            [max_error_3,I_3] = max(abs(max(model_mat,[],1) - f(plot_arr)));
            disp(['max |R(z)| = ',num2str(max_error_3),' at z = ',num2str(plot_arr(I_3))])
    end

    set(gca,'linewidth',1.2)
    set(gca,'fontsize',18)
    grid on
    box on
    xlim([plot_min,plot_max])
    ylim([-0.4,2.6])

    xlabel('$x$','Interpreter','latex')

    if(save_figure_flag)
        fname = ['plot_',num2str(plot_i)];
        export_fig(fname,'-png','-r500','-transparent','-painters')
        close(fig)
    end

end
