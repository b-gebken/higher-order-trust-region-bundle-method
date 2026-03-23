% Script that creates the plots shown in Fig. 1.

save_plot_flag = 0;

%%
x_sym = sym('x'); assume(x_sym,'real');

f = @(x,p) x.^(p+1) .* sin(1./x) + 1/p * abs(x).^p;

x_lower = -0.08;
x_upper = 0.08;

x_arr = linspace(x_lower,x_upper,10000);
tmp = [-1./(2*pi*(2:10000)),1./(2*pi*(2:10000))];
x_arr = sort([x_arr,tmp]);

lw = 1.5;

for p = [1,2,4]
    figure

    plot(x_arr,f(x_arr,p),'k-','LineWidth',lw);

    grid on
    xlim([x_lower,x_upper]);
    xticks(-0.1:0.05:0.1);

    set(gca,'linewidth',1.1)
    set(gca,'fontsize',15)

    ylim_vals = ylim;
    text(-0.045,0.9*ylim_vals(2),['$p = ',num2str(p),'$'],'Interpreter','latex','HorizontalAlignment','center','FontSize',24);

    axis square

    if(p == 4)
        % Subfigure -----------------------------------------------------------

        ax = axes('position',[.47 .6 .3 .3]);

        min_zoom = 0.052050644954723;
        max_zoom = 0.055616954831905;
        zoom_arr = linspace(min_zoom,max_zoom,1000);

        plot(zoom_arr,f(zoom_arr,p),'k-','LineWidth',lw);
        hold on

        xlim([min_zoom,max_zoom])
        axis square
        ax.ClippingStyle = "rectangle";

        xticks([0.053,0.055])
        yticks([])

        set(gca,'linewidth',1.1)
        set(gca,'fontsize',10)
    end

    if(save_plot_flag)
        export_fig(['plot_p_',num2str(p)],'-png','-r500','-transparent')
    end
end