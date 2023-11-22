%%% this code plot Phit and SAT for groups
clear all;
clc;            
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));

mks = 6; lw = 0.8;

load Init_Inhb_Analysis.mat;  % Analysis
D = Analysis.D;
A = Analysis.A;
%head(DATA,3)  % check if data look right
xplot= D.xplot;

%% Figure 2D
% plot model vs data; group level 
SAT_gono = [];
SAT_nogo = [];

Model_gono = [];
Model_nogo = [];

SUB = [1:15 17:size(D.sub_name,2)]; % subject 16 was excluded; see manuscript

for s = 1:length(SUB)
    SAT_gono(:,s) = A.p_gono_nolate{SUB(s)}; % Speed-accuracy tradeoff through sliding window; for visualization
    SAT_nogo(:,s) = A.p_nogo_nolate{SUB(s)};
    Model_gono(:,s) = A.ycdf_gono_nolate{SUB(s)}; % y hat predicted from model
    Model_nogo(:,s) = A.ycdf_nogo_nolate{SUB(s)};
end

figure_SAT = figure('name','SAT_model_nolate');  
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
set(gca,'fontsize',10)
xlabel('Allowed RT(ms)','FontSize',18, 'FontWeight','normal');
ylabel('Probability','FontSize',18, 'FontWeight','normal');
set(gca,'xTick',[0:0.1:0.5],'xTickLabel',[0:100:500],'FontSize',18, 'FontWeight','normal','FontName','Arial');
axis([0 0.5, 0, 1])
hAx=gca;                    % create an axes
hAx.LineWidth=1.2; 

% plot standard error bar for data
shadedErrorBar(xplot,nanmean(SAT_gono,2),seNaN(SAT_gono),{'-','color','b'});
shadedErrorBar(xplot,nanmean(SAT_nogo,2),seNaN(SAT_nogo),{'-','color','r'});

% data
f2 = plot(xplot,nanmean(SAT_nogo,2),'r--','markersize',mks+6,'LineWidth',1,'MarkerFaceColor','w','LineStyle','--');
f1 = plot(xplot,nanmean(SAT_gono,2),'b--','markersize',mks+6,'LineWidth',1,'MarkerFaceColor','w','LineStyle','--');

% model prediction
f3 = plot(xplot,nanmean(Model_gono,2),'b-','markersize',mks,'LineWidth',1.5,'MarkerFaceColor','b');
f4 = plot(xplot,nanmean(Model_nogo,2),'r-','markersize',mks,'LineWidth',1.5,'MarkerFaceColor','r');

f3.Color(4) = 0.9;
f4.Color(4) = 0.9;

legend([f1, f2,f3,f4],{'R-to-NR (data)' ,'NR-to-R (data)', 'R-to-NR (model)','NR-to-R (model)'},'Location','Northwest','NumColumns',1,...
             'fontsize',14,'textcolor','k');
            legend('boxoff');
            
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
%      'PaperUnits', 'Inches', 'PaperSize', [5, 5])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal\analysis\Figure');
% print(figure_SAT,'SAT_model_nolate_adjust', '-depsc','-r600');
% cd('D:\Project\StopSignal\analysis\Manuscript\NewAnalysis_NoLate');


%% Figure 2E
% plot mu_go vs mu_nogo
for s = 1:length(SUB)
    mu_gono(s) = A.model_gono_nolate{SUB(s)}(1); % estimated value for model parameter
    mu_nogo(s) = A.model_nogo_nolate{SUB(s)}(1); % (1) is \mu
end

figure_mu = figure('name','mu');  
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
axis([0.2 0.37 0.2 0.37])
hAx=gca;                    % create an axes
hAx.LineWidth=1.2; 
xlabel('\mu_{R} (ms)','FontSize',22, 'FontWeight','bold');
ylabel('\mu_{NR} (ms)','FontSize',22, 'FontWeight','bold');
pbaspect([1 1 1])
set(gca,'xTick',[0:0.05:0.5],'xTickLabel',[0:50:500],'FontSize',22, 'FontWeight','normal','FontName','Arial');
set(gca,'yTick',[0:0.05:0.5],'yTickLabel',[0:50:500],'FontSize',22, 'FontWeight','normal','FontName','Arial');

% diagonal line
plot([0.12 0.37],[0.12,0.37],'k-','markersize',mks,'LineWidth',1.5,'MarkerFaceColor','k')

% data
f1 = plot(mu_nogo, mu_gono,'ko','markersize',mks+8,'LineWidth',1.5,'MarkerFaceColor','w');

text(0.3,0.25,'\mu_{NR} - \mu_{R} = 5.5ms','FontSize',22,'FontWeight','normal');
text(0.3,0.22,'95% CI: [-2.8, 13.9]ms','FontSize',22,'FontWeight','normal');

% to get stats and ci
% nanmean(mu_gono' - mu_nogo')
% [h,p,ci,stats] = ttest(mu_gono,mu_nogo)

% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
%      'PaperUnits', 'Inches', 'PaperSize', [5, 5])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal\analysis\Figure');
% print(figure_mu,'mu_nolate_adjust', '-depsc','-r600');