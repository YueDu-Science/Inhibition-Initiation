function f = plot_individual_DTST_simulation_mle(P_rec, RTsim, D)
% Initialization
% color
%cl = [0 0 1; 1 0 0; 0.8516 0.6445 0.1250; 0.9297 0.5078 0.9297; 0.5 0.5 0]; 
cl = [0 0 1; 1 0 0; 0.9297 0.5078 0.9297; 0.5 0.5 0]; 
%c = [176,224,230; 1 0 0; 0.9297 0.5078 0.9297; 0.5 0.5 0]; 

mks = 3; lw = 1.4;
cols = [0 0 0; 0 100 255; 255 0 0; 100 0 255]/256;


% model;
Xfit = P_rec.Xfit;
st = P_rec.st;

%head(DATA,3)  % check if data look right
xplot= D.xplot;
% deal with individual subject's data
figure_individual_1 = figure('name','Individual_1');
figure_individual_2 = figure('name','Individual_2');    
for s = 1:40
    data = [];
    data = RTsim(s,:)';
    if s <= 20
        figure(figure_individual_1);

        subplot(4,5,s)
        set(gcf,'color','w');
        hold on
        sub_tt = s;
        set(gca,'TickDir','out');
        set(gca,'fontsize',10)
        axis([0 0.52 0.05 0.8])
        set(gca,'xTick',[0:0.1:0.5],'xTickLabel',[0:100:500],'FontSize',10, 'FontWeight','normal','FontName','Arial');
        set(gca,'yTick',[0:0.1:0.8],'yTickLabel',[0:100:800],'FontSize',10, 'FontWeight','normal','FontName','Arial');
        title(sub_tt,'FontSize',12, 'FontWeight','normal');
        xlabel('DT(ms)','FontSize',12, 'FontWeight','normal');
        ylabel('Resp Time(ms)','FontSize',12, 'FontWeight','normal');
        
%         ind = find(data.t_choice ~= -1);
%         plot(data.t_prep(ind),data.t_choice(ind),'ro','markersize',mks,'LineWidth',lw,'MarkerFaceColor','r');
        ind = find(~isnan(data));
        plot(st(ind), data(ind),'ro','markersize',mks,'LineWidth',lw,'MarkerFaceColor','w');
        
%         ind = find(data.block_key == min(data.block_key));
%         plot(data.t_prep(ind),data.button(ind),'bo','markersize',mks,'LineWidth',lw,'MarkerFaceColor','b');
        
        % smooth t_choice
%         t_prep = data.t_prep(ind);
%         t_choice = smooth(t_prep,data.t_choice(ind),0.4,'rloess');
%         [B,I] = sort(t_prep);
%         plot(t_prep(I),t_choice(I),'k-','markersize',mks,'LineWidth',lw,'MarkerFaceColor','k');
        
        %find out first index when Resp Time enter the target zone
        %t_prep_sort = t_prep(I); t_choice_sort = t_choice(I);

        %         for i = 2:(numel(t_choice_sort)-1)
%             if t_choice_sort(i-1) > 0.53 && t_choice_sort(i) <= 0.53 && t_choice_sort(i+1) <= 0.53
%                 cut_off = t_prep_sort(i);
%                 break
%             end
%         end
%         cut_off = mu_nogo(s);
%         plot([cut_off, cut_off],[0,1],'k--','LineWidth',lw - 0.4)
        
        ind = find(isnan(data));
        data(ind) = 0.12;
        rng('default'); ss = rng; rng(ss);
        tmp1 = 0 + normrnd(0,0.01,numel(ind),1);
        scatter(st(ind),data(ind)+tmp1,20,'k','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
        
       
        st_sort = sort(st);
        fit_nogo_rt_1 = (linspace(0.05,Xfit(s,1),10) - Xfit(s,1))*Xfit(s,5) + Xfit(s,2);
        fit_nogo_rt_2 = (linspace(Xfit(s,1),0.5,10) - Xfit(s,1))*Xfit(s,7) + Xfit(s,2);

        plot(linspace(0.05,Xfit(s,1),10), fit_nogo_rt_1, 'k-', 'linewidth',2)
        plot(linspace(Xfit(s,1),0.5,10), fit_nogo_rt_2, 'k-', 'linewidth',2);
        plot([Xfit(s,1) Xfit(s,1)], [0.2 0.4],'b-','linewidth',2);
    
%         if s == 1
%             legend([f1, f2],{'Original' ,'ExGau-ST'},'Location',[0.05, 0.05, 0.05 , 0.05],'NumColumns',1,...
%                  'fontsize',12,'textcolor','k');
%                 legend('boxoff');
%         end
    elseif s > 20
        figure(figure_individual_2);

        subplot(4,5,s-20)
        set(gcf,'color','w');
        hold on
        sub_tt = s;
        set(gca,'TickDir','out');
        set(gca,'fontsize',10)
        axis([0 0.52 0.05 0.8])
        set(gca,'xTick',[0:0.1:0.5],'xTickLabel',[0:100:500],'FontSize',10, 'FontWeight','normal','FontName','Arial');
        set(gca,'yTick',[0.1:0.1:0.8],'yTickLabel',[100:100:800],'FontSize',10, 'FontWeight','normal','FontName','Arial');
        title(sub_tt,'FontSize',12, 'FontWeight','normal');
        xlabel('DT(ms)','FontSize',12, 'FontWeight','normal');
        ylabel('Resp Time(ms)','FontSize',12, 'FontWeight','normal');
        
        
        ind = find(~isnan(data));
        plot(st(ind),data(ind),'ro','markersize',mks,'LineWidth',lw,'MarkerFaceColor','w');
        
%         ind = find(data.block_key == min(data.block_key));
%         plot(data.t_prep(ind),data.button(ind),'bo','markersize',mks,'LineWidth',lw,'MarkerFaceColor','b');
        % smooth t_choice
%         t_prep = data.t_prep(ind);
%         t_choice = smooth(t_prep,data.t_choice(ind),0.4,'rloess');
%         [B,I] = sort(t_prep);
%         plot(t_prep(I),t_choice(I),'k-','markersize',mks,'LineWidth',lw,'MarkerFaceColor','k');
        
        %find out first index when Resp Time enter the target zone
        %t_prep_sort = t_prep(I); t_choice_sort = t_choice(I);
%         for i = 2:(numel(t_choice_sort)-1)
%             if t_choice_sort(i-1) > 0.53 && t_choice_sort(i) <= 0.53 && t_choice_sort(i+1) <= 0.53
%                 cut_off = t_prep_sort(i);
%                 break
%             end
%         end
%         cut_off = mu_nogo(s);
%         plot([cut_off, cut_off],[0,1],'k--','LineWidth',lw - 0.4)
        
        ind = find(isnan(data));
        data(ind) = 0.12;
        rng('default'); ss = rng; rng(ss);
        tmp1 = 0 + normrnd(0,0.01,numel(ind),1);
        scatter(st(ind),data(ind)+tmp1,20,'k','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
        
%         % plot off-diagonal lines
%         x = 0:0.1:0.8;
%         intercept = 0;
%         for jj = 1:20
%             y = intercept - x;
%             plot(x,y,'b.-')
%             intercept = intercept + 0.05;
%         end
       
        st_sort = sort(st);
        fit_nogo_rt_1 = (linspace(0.05,Xfit(s,1),10) - Xfit(s,1))*Xfit(s,5) + Xfit(s,2);
        fit_nogo_rt_2 = (linspace(Xfit(s,1),0.5,10) - Xfit(s,1))*Xfit(s,7) + Xfit(s,2);

        plot(linspace(0.05,Xfit(s,1),10), fit_nogo_rt_1, 'k-', 'linewidth',2)
        plot(linspace(Xfit(s,1),0.5,10), fit_nogo_rt_2, 'k-', 'linewidth',2);
        plot([Xfit(s,1) Xfit(s,1)], [0.35 0.6],'b-','linewidth',2);
        
%         if s == 21
%             legend([f1, f2],{'Original' ,'ExGau-ST'},'Location',[0.05, 0.05, 0.05 , 0.05],'NumColumns',1,...
%                  'fontsize',12,'textcolor','k');
%                 legend('boxoff');
%         end
    end
end
%%
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 8],...
%      'PaperUnits', 'Inches', 'PaperSize', [8, 8])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal\analysis\Figure');
% print(figure_individual_1,'PrepTime_individual_1_buttonclick', '-depsc','-r600');
% 
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 8],...
%      'PaperUnits', 'Inches', 'PaperSize', [8, 8])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal\analysis\Figure');
% print(figure_individual_2,'PrepTime_individual_2_buttonclick', '-depsc','-r600');
% cd ..




% %%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% this code plot Phit for individuals and groups
% clear all;
% clc;            
% addpath('D:\Project\StopSignal\analysis');
% 
% % Initialization
% % color
% %cl = [0 0 1; 1 0 0; 0.8516 0.6445 0.1250; 0.9297 0.5078 0.9297; 0.5 0.5 0]; 
% cl = [0 0 1; 1 0 0; 0.9297 0.5078 0.9297; 0.5 0.5 0]; 
% %c = [176,224,230; 1 0 0; 0.9297 0.5078 0.9297; 0.5 0.5 0]; 
% 
% mks = 3; lw = 1.4;
% cols = [0 0 0; 0 100 255; 255 0 0; 100 0 255]/256;
% 
% %DATA = readtable('tr-training.txt');
% load Ready_for_Use.mat;  % D
% 
% %head(DATA,3)  % check if data look right
% xplot= D.xplot;
% % deal with individual subject's data
% figure_individual_5 = figure('name','Individual_5');
% figure_individual_6 = figure('name','Individual_6');    
% for s = 1:length(D.sub_name)
%     data = [];
%     data = D.nogo{s};
%     if s <= 20
%         figure(figure_individual_5);
% 
%         subplot(4,5,s)
%         set(gcf,'color','w');
%         hold on
%         sub_tt = D.sub_name(s);
%         set(gca,'TickDir','out');
%         set(gca,'fontsize',10)
%         axis([0 0.52 0 0.6])
%         set(gca,'xTick',[0:0.1:0.5],'xTickLabel',[0:100:500],'FontSize',10, 'FontWeight','normal','FontName','Arial');
%         set(gca,'yTick',[0:0.1:0.6],'yTickLabel',[0:100:6800],'FontSize',10, 'FontWeight','normal','FontName','Arial');
%         title(sub_tt,'FontSize',12, 'FontWeight','normal');
%         xlabel('DT/ST(ms)','FontSize',12, 'FontWeight','normal');
%         ylabel('RT(ms)','FontSize',12, 'FontWeight','normal');
%         pbaspect([1 1 1])
% %         ind = find(data.t_choice ~= -1);
% %         plot(data.t_prep(ind),data.t_choice(ind),'ro','markersize',mks,'LineWidth',lw,'MarkerFaceColor','r');
%         ind = find(data.button ~= -99);
%         plot(data.t_prep(ind),data.button(ind) - 0.5 + data.t_prep(ind),'ro','markersize',mks,'LineWidth',lw,'MarkerFaceColor','r');
%         
%         % smooth t_choice
% %         t_prep = data.t_prep(ind);
% %         t_choice = smooth(t_prep,data.t_choice(ind),0.4,'rloess');
% %         [B,I] = sort(t_prep);
% %         plot(t_prep(I),t_choice(I),'k-','markersize',mks,'LineWidth',lw,'MarkerFaceColor','k');
%         
%         %find out first index when Resp Time enter the target zone
%         %t_prep_sort = t_prep(I); t_choice_sort = t_choice(I);
% 
%         %         for i = 2:(numel(t_choice_sort)-1)
% %             if t_choice_sort(i-1) > 0.53 && t_choice_sort(i) <= 0.53 && t_choice_sort(i+1) <= 0.53
% %                 cut_off = t_prep_sort(i);
% %                 break
% %             end
% %         end
% %         cut_off = mu_nogo(s);
% %         plot([cut_off, cut_off],[0,1],'k--','LineWidth',lw - 0.4)
%         
%         ind = find(data.button == -99);
%         data.button(ind) = 0.02;
%         rng('default'); s = rng; rng(s);
%         tmp1 = 0 + normrnd(0,0.01,numel(ind),1);
%         scatter(data.t_prep(ind),data.button(ind)+tmp1,36,'r','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
%         
%         
%         % plot off-diagonal lines
%         x = 0:0.1:0.6;
%         intercept = 0;
%         for jj = 1:20
%             y = repmat(intercept,1,numel(x));
%             plot(x,y,'b.-')
%             intercept = intercept + 0.05;
%         end
%     elseif s > 20
%         figure(figure_individual_6);
% 
%         subplot(4,5,s-20)
%         set(gcf,'color','w');
%         hold on
%         sub_tt = D.sub_name(s);
%         set(gca,'TickDir','out');
%         set(gca,'fontsize',10)
%         axis([0 0.52 0 0.6])
%         set(gca,'xTick',[0:0.1:0.5],'xTickLabel',[0:100:500],'FontSize',10, 'FontWeight','normal','FontName','Arial');
%         set(gca,'yTick',[0:0.1:0.6],'yTickLabel',[0:100:6800],'FontSize',10, 'FontWeight','normal','FontName','Arial');
%         title(sub_tt,'FontSize',12, 'FontWeight','normal');
%         xlabel('DT/ST(ms)','FontSize',12, 'FontWeight','normal');
%         ylabel('RT(ms)','FontSize',12, 'FontWeight','normal');
%         pbaspect([1 1 1])
%         
%         ind = find(data.button ~= -99);
%         plot(data.t_prep(ind),data.button(ind) - 0.5 + data.t_prep(ind),'ro','markersize',mks,'LineWidth',lw,'MarkerFaceColor','r');
%         
%         % smooth t_choice
% %         t_prep = data.t_prep(ind);
% %         t_choice = smooth(t_prep,data.t_choice(ind),0.4,'rloess');
% %         [B,I] = sort(t_prep);
% %         plot(t_prep(I),t_choice(I),'k-','markersize',mks,'LineWidth',lw,'MarkerFaceColor','k');
%         
%         %find out first index when Resp Time enter the target zone
%         %t_prep_sort = t_prep(I); t_choice_sort = t_choice(I);
% %         for i = 2:(numel(t_choice_sort)-1)
% %             if t_choice_sort(i-1) > 0.53 && t_choice_sort(i) <= 0.53 && t_choice_sort(i+1) <= 0.53
% %                 cut_off = t_prep_sort(i);
% %                 break
% %             end
% %         end
% %         cut_off = mu_nogo(s);
% %         plot([cut_off, cut_off],[0,1],'k--','LineWidth',lw - 0.4)
%         
%         ind = find(data.button == -99);
%         data.button(ind) = 0.02;
%         rng('default'); s = rng; rng(s);
%         tmp1 = 0 + normrnd(0,0.01,numel(ind),1);
%         scatter(data.t_prep(ind),data.button(ind)+tmp1,36,'r','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
%         
%         % plot off-diagonal lines
%         x = 0:0.1:0.6;
%         intercept = 0;
%         for jj = 1:20
%             y = repmat(intercept,1,numel(x));
%             plot(x,y,'b.-')
%             intercept = intercept + 0.05;
%         end
%     end
% end
% 
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
%      'PaperUnits', 'Inches', 'PaperSize', [5, 5])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal\analysis\Figure');
% print(figure_individual_5,'PrepTime_RT_individual_4', '-depsc','-r600');
% 
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
%      'PaperUnits', 'Inches', 'PaperSize', [5, 5])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal\analysis\Figure');
% print(figure_individual_6,'PrepTime_RT_individual_6', '-depsc','-r600');
% cd ..