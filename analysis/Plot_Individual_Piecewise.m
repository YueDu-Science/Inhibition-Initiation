%%% this code plot Phit for individuals and groups
clear all;
close all;
clc;            
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));

mks = 3; lw = 0.3; fs = 6;

load Init_Inhb_Analysis.mat;  % Analysis
D = Analysis.D;
A = Analysis.A;

load ModelFitting_PW_MLE;
mu = M_PW.Xfit(:,1);
Xfit = M_PW.Xfit;
st = M_PW.st;
INDX1 = M_PW.INDX1;
INDX2 = M_PW.INDX2;
xplot= D.xplot;

% deal with individual subject's data
figure_individual_1 = figure('name','Individual_1');
figure_individual_2 = figure('name','Individual_2');    
for s = 1:length(D.sub_name)
    data = D.nogo{s};
    indx1 = INDX1{s}; % index does satisfy condition 1: see 4_Fitting_Piecewise_Model
    indx2 = INDX2{s}; % index does satisfy condition 2: see 4_Fitting_Piecewise_Model
    if s <= 20
        figure(figure_individual_1);

        subplot(4,5,s)
        set(gcf,'color','w');
        hold on
        sub_tt = D.sub_name(s);
        set(gca,'TickDir','out');
        set(gca,'fontsize',10)
        axis([0 0.52 0.05 0.8])
        set(gca,'xTick',[0:0.1:0.5],'xTickLabel',[0:100:500],'FontSize',fs, 'FontWeight','normal','FontName','Arial');
        set(gca,'yTick',[0.1:0.1:0.8],'yTickLabel',[100:100:800],'FontSize',fs, 'FontWeight','normal','FontName','Arial');
        %title(sub_tt,'FontSize',fs, 'FontWeight','normal');
        if s == 1
            xlabel('Intended RT(ms)','FontSize',fs, 'FontWeight','normal');
            ylabel('Time of Response(ms)','FontSize',fs, 'FontWeight','normal');
        end
%         ind = find(data.t_choice ~= -1);
%         plot(data.t_prep(ind),data.t_choice(ind),'ro','markersize',mks,'LineWidth',lw,'MarkerFaceColor','r');
        ind = find(data.button ~= -99);
        f1 = plot(data.t_prep(ind),data.button(ind),'ro','markersize',mks,'LineWidth',lw,'MarkerFaceColor','w');
        f1.Color(4) = 0.5;
        g1 = plot(data.t_prep(indx1),data.button(indx1),'kx','markersize',mks+2,'LineWidth',lw,'MarkerFaceColor','w');
        g2 = plot(data.t_prep(indx2),data.button(indx2),'k+','markersize',mks+2,'LineWidth',lw,'MarkerFaceColor','w');
        
        ind = find(data.button == -99); % no response trial; just plot for reference
        data.button(ind) = 0.12;
        rng('default'); ss = rng; rng(ss);
        tmp1 = 0 + normrnd(0,0.01,numel(ind),1);
        f2 = scatter(data.t_prep(ind),data.button(ind)+tmp1,10,'k','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);

        % tolerance window (470 to 530 ms)
        y = 0:0.01:0.52;
        plot(y, repmat(0.5,numel(y),1),'k-','LineWidth',1)
        plot(y, repmat(0.62,numel(y),1),'k--','LineWidth',1)
        curve1 = repmat(0.47,1,length(y));
        curve2 = repmat(0.53,1,length(y));
        y2 = [y, fliplr(y)];
        inBetween = [curve1, fliplr(curve2)];
        h = fill(y2, inBetween, [0.5,0.5,0.5]);
        h.EdgeColor = [0.5 0.5 0.5];
        uistack(h,'bottom')

        % model prediction
        % y hat from the piecewise model
        % fit_nogo_rt_2 vs fit_nogo_rt_1_total: one is the normal component of the
        % exGaussian distribution, the other consistes of both normal and
        % exponential components
        st_sort = sort(st);
        fit_nogo_rt_1 = (linspace(0.05,Xfit(s,1),10) - Xfit(s,1))*Xfit(s,5) + Xfit(s,2);
        fit_nogo_rt_2 = (linspace(Xfit(s,1),0.5,10) - Xfit(s,1))*Xfit(s,7) + Xfit(s,2);
        fit_nogo_rt_3 = (linspace(0.05,Xfit(s,1),10) - Xfit(s,1))*Xfit(s,5) + Xfit(s,2) + Xfit(s,4);

       % plot model prediction
       f3 = plot(linspace(0.05,Xfit(s,1),10), fit_nogo_rt_1, '-', 'linewidth',1.2,'color',[30,144,255]/255);
       plot(linspace(Xfit(s,1),0.5,10), fit_nogo_rt_2, '-', 'linewidth',1.2,'color',[30,144,255]/255);
       f4 = plot(linspace(0.05,Xfit(s,1),10), fit_nogo_rt_3, '--', 'linewidth',1.2,'color',[30,144,255]/255);
       f4.Color(4) = 0.3;
       plot(Xfit(s,1), Xfit(s,2),'o','color',[30,144,255]/255,'linewidth',1.2,'MarkerFaceColor', [30,144,255]/255,'markersize',mks+5);

        if s == 20
            legend([f1, f3, f4, f2],{'Data' ,'Model (Gaussian component)','Model (Gaussian + Exponential components)','No Response Trials'},'Location',[0.05, 0.05, 0.05 , 0.05],'NumColumns',1,...
                 'fontsize',8,'textcolor','k');
                legend('boxoff');
        end
    elseif s > 20
        figure(figure_individual_2);

        subplot(4,5,s-20)
        set(gcf,'color','w');
        hold on
        sub_tt = D.sub_name(s);
        set(gca,'TickDir','out');
        set(gca,'fontsize',10)
        axis([0 0.52 0.05 0.8])
        set(gca,'xTick',[0:0.1:0.5],'xTickLabel',[0:100:500],'FontSize',fs, 'FontWeight','normal','FontName','Arial');
        set(gca,'yTick',[0.1:0.1:0.8],'yTickLabel',[100:100:800],'FontSize',fs, 'FontWeight','normal','FontName','Arial');
        %title(sub_tt,'FontSize',fs, 'FontWeight','normal');
        if s == 21
            xlabel('DT(ms)','FontSize',fs, 'FontWeight','normal');
            ylabel('Resp Time(ms)','FontSize',fs, 'FontWeight','normal');
        end
        ind = find(data.button ~= -99);
        f1 = plot(data.t_prep(ind),data.button(ind),'ro','markersize',mks,'LineWidth',lw,'MarkerFaceColor','w');
        f1.Color(4) = 0.5;
        g1 = plot(data.t_prep(indx1),data.button(indx1),'kx','markersize',mks+2,'LineWidth',lw,'MarkerFaceColor','w');
        g2 = plot(data.t_prep(indx2),data.button(indx2),'k+','markersize',mks+2,'LineWidth',lw,'MarkerFaceColor','w');
        
%         ind = find(data.block_key == min(data.block_key));
%         plot(data.t_prep(ind),data.button(ind),'bo','markersize',mks,'LineWidth',lw,'MarkerFaceColor','b');

        ind = find(data.button == -99);
        data.button(ind) = 0.12;
        rng('default'); ss = rng; rng(ss);
        tmp1 = 0 + normrnd(0,0.01,numel(ind),1);
        f2 = scatter(data.t_prep(ind),data.button(ind)+tmp1,10,'k','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);

        
        y = 0:0.01:0.52;
        plot(y, repmat(0.5,numel(y),1),'k-','LineWidth',1)
        plot(y, repmat(0.62,numel(y),1),'k--','LineWidth',1)
        curve1 = repmat(0.47,1,length(y));
        curve2 = repmat(0.53,1,length(y));
        y2 = [y, fliplr(y)];
        inBetween = [curve1, fliplr(curve2)];
        h = fill(y2, inBetween, [0.5,0.5,0.5]);
        h.EdgeColor = [0.5 0.5 0.5];
        uistack(h,'bottom')

        st_sort = sort(st);
        fit_nogo_rt_1 = (linspace(0.05,Xfit(s,1),10) - Xfit(s,1))*Xfit(s,5) + Xfit(s,2) ;
        fit_nogo_rt_2 = (linspace(Xfit(s,1),0.5,10) - Xfit(s,1))*Xfit(s,7) + Xfit(s,2);
        fit_nogo_rt_3 = (linspace(0.05,Xfit(s,1),10) - Xfit(s,1))*Xfit(s,5) + Xfit(s,2) + Xfit(s,4);

        f2 = plot(linspace(0.05,Xfit(s,1),10), fit_nogo_rt_1, '-', 'linewidth',1.2,'color',[30,144,255]/255)
        plot(linspace(Xfit(s,1),0.5,10), fit_nogo_rt_2, '-', 'linewidth',1.2,'color',[30,144,255]/255);
        f4 = plot(linspace(0.05,Xfit(s,1),10), fit_nogo_rt_3, '--', 'linewidth',1.2,'color',[30,144,255]/255);
        f4.Color(4) = 0.3;
        plot(Xfit(s,1), Xfit(s,2),'o','color',[30,144,255]/255,'linewidth',1.2,'MarkerFaceColor', [30,144,255]/255,'markersize',mks+5);

        if s == 21
            legend([f1, f3, f4, f2],{'Data' ,'Model (Gaussian component)','Model (Gaussian + Exponential components)','No Response Trials'},'Location',[0.05, 0.05, 0.05 , 0.05],'NumColumns',1,...
                 'fontsize',8,'textcolor','k');
                legend('boxoff');
        end
    end
    clear data
end
%%
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 5],...
     'PaperUnits', 'Inches', 'PaperSize', [8, 5])
set(gcf,'renderer','Painters');
cd('D:\Project\StopSignal_Analysis\analysis\Figure');
print(figure_individual_1,'PrepTime_individual_1_buttonclick', '-depsc','-r600');

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 5],...
     'PaperUnits', 'Inches', 'PaperSize', [8, 5])
set(gcf,'renderer','Painters');
cd('D:\Project\StopSignal_Analysis\analysis\Figure');
print(figure_individual_2,'PrepTime_individual_2_buttonclick', '-depsc','-r600');
