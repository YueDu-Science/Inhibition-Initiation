% plot individual speed-accuracy tradeoff
%% Figure S1
clear all;
clc;            
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));
mks = 6; lw = 0.8;

load Init_Inhb_Analysis.mat;  % Analysis
D = Analysis.D;
A = Analysis.A;
%head(A,3)  % check if data look right
xplot= D.xplot;

SAT_gono = [];
SAT_nogo = [];
Model_gono = [];
Model_nogo = [];

SUB = [1:15 16 17:size(D.sub_name,2)]; % also plot sub 16 who was excluded for analysis

for s = 1:length(SUB)
    SAT_gono(:,s) = A.p_gono_nolate{SUB(s)}; 
    SAT_nogo(:,s) = A.p_nogo_nolate{SUB(s)};
    Model_gono(:,s) = A.ycdf_gono_nolate{SUB(s)};
    Model_nogo(:,s) = A.ycdf_nogo_nolate{SUB(s)};
end

mks = 3; lw = 0.6;
figure_SAT_individual_1 = figure('name','STOP_Individual_1');
figure_SAT_individual_2 = figure('name','STOP_Individual_2');    

for s = 1:length(D.sub_name)
        if s <= 20
            figure(figure_SAT_individual_1);
            subplot(4,5,s)
            set(gcf,'color','w');
            hold on
            if s == 1
                xlabel('Allowed RT(ms)','FontSize',10, 'FontWeight','normal');
                ylabel('Probability','FontSize',10, 'FontWeight','normal');
            end
            axis([0 0.5 0 1])
            set(gca,'xTick',[0:0.25:0.5],'xTickLabel',[0:250:500],'FontSize',6, 'FontWeight','normal','FontName','Arial');
            set(gca,'TickDir','out');
            set(gca,'fontsize',6)
            
            plot(xplot,SAT_gono(:,s),'b-','LineWidth',lw,'MarkerFaceColor','w');
            plot(xplot,SAT_nogo(:,s),'r-','LineWidth',lw,'MarkerFaceColor','w');
            
            f3 = plot(xplot,Model_gono(:,s),'b-','markersize',mks,'LineWidth',lw + 0.6,'MarkerFaceColor','b');
            f4 = plot(xplot,Model_nogo(:,s),'r-','markersize',mks,'LineWidth',lw + 0.6,'MarkerFaceColor','b');
            f3.Color(4) = 0.6;
            f4.Color(4) = 0.6;

        elseif s > 20
            figure(figure_SAT_individual_2);
            subplot(4,5,s-20)
            set(gcf,'color','w');
            hold on

            axis([0 0.5 0 1])
            set(gca,'xTick',[0:0.25:0.5],'xTickLabel',[0:250:500],'FontSize',6, 'FontWeight','normal','FontName','Arial');
            set(gca,'TickDir','out');
            set(gca,'fontsize',6)
            
            f1 = plot(xplot,SAT_gono(:,s),'b-','LineWidth',lw,'MarkerFaceColor','w');
            f2 = plot(xplot,SAT_nogo(:,s),'r-','LineWidth',lw,'MarkerFaceColor','w');
            
            f3 = plot(xplot,Model_gono(:,s),'b-','markersize',mks,'LineWidth',lw + 0.6,'MarkerFaceColor','b');
            f4 = plot(xplot,Model_nogo(:,s),'r-','markersize',mks,'LineWidth',lw + 0.6,'MarkerFaceColor','b');
            f3.Color(4) = 0.6;
            f4.Color(4) = 0.6;
           
            if s == size(SAT_gono,2)
                legend([f1,f2,f3,f4],{'R-to-NR (data)' ,'NR-to-R (data)', 'R-to-NR (model)','NR-to-R (model)'},'Location',[0.5,0.12,0.05,0.05],'NumColumns',1,...
                'fontsize',14,'textcolor','k');
                legend('boxoff');
            end
            
        end
end

% figure(figure_SAT_individual_1);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 5],...
%      'PaperUnits', 'Inches', 'PaperSize', [8, 5])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal\analysis\Figure');
% print(figure_SAT_individual_1,'SAT_model_individual_1_nolate_adjust', '-depsc','-r600');
% 
% figure(figure_SAT_individual_2);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 5],...
%      'PaperUnits', 'Inches', 'PaperSize', [8, 5])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal\analysis\Figure');
% print(figure_SAT_individual_2,'SAT_model_individual_2_nolate_adjust', '-depsc','-r600');