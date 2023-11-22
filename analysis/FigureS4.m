%% Figure S3: how the way calculate accuracy of nogo trials affects the SAT
clear all;
clc;            
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));
mks = 6; lw = 0.8;

load Init_Inhb_Analysis.mat;  % Analysis
D = Analysis.D;
A = Analysis.A;
%head(DATA,3)  % check if data look right
xplot= D.xplot;

c = copper(100); % the 37 row is red
SUB = [1:15 17:size(D.sub_name,2)]; % sub 16 was excluded

figure_SAT = figure('name','SAT_model_nolate');  
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
set(gca,'fontsize',10)
%title(sub_tt,'FontSize',12, 'FontWeight','normal');
xlabel('Allowed RT (ms)','FontSize',18, 'FontWeight','normal');
ylabel('Probability','FontSize',18, 'FontWeight','normal');
set(gca,'xTick',[0:0.1:0.5],'xTickLabel',[0:100:500],'FontSize',18, 'FontWeight','normal','FontName','Arial');
axis([0 0.5, 0, 1])
hAx=gca;                    % create an axes
hAx.LineWidth=1.2; 
SUB = [1:15 17:size(D.sub_name,2)];
%1st criterion: original
for s = 1:length(SUB)
    SAT_gono(:,s) = A.p_gono_nolate{SUB(s)}; 
    SAT_nogo(:,s) = A.p_nogo_nolate{SUB(s)};
    Model_gono(:,s) = A.ycdf_gono_nolate{SUB(s)};
    Model_nogo(:,s) = A.ycdf_nogo_nolate{SUB(s)};
end

f(2) = plot(xplot,nanmean(SAT_nogo,2),'--','color','r','markersize',mks,'LineWidth',1,'MarkerFaceColor','w');
f(2) = plot(xplot,nanmean(Model_nogo,2),'-','color','r','markersize',mks,'LineWidth',2.5,'MarkerFaceColor','r');
f(2).Color(4) = 0.9;
f(1) = plot(xplot,nanmean(SAT_gono,2),'--','color','b','markersize',mks,'LineWidth',1,'MarkerFaceColor','w');
f(1) = plot(xplot,nanmean(Model_gono,2),'-','color','b','markersize',mks,'LineWidth',2.5,'MarkerFaceColor','r');
f(1).Color(4) = 0.9;

% 2nd criterion: no timing requirement
SAT_nogo = [];Model_nogo = [];

for s = 1:length(SUB)
    SAT_nogo(:,s) = A.p_nogo_button{SUB(s)};
    Model_nogo(:,s) = A.ycdf_nogo_button{SUB(s)};
end

f(6) = plot(xplot,nanmean(SAT_nogo,2),'--','color',c(1,:),'markersize',mks+3,'LineWidth',1,'MarkerFaceColor','w');
f(6) = plot(xplot,nanmean(Model_nogo,2),'-','color',c(1,:),'markersize',mks,'LineWidth',2.5,'MarkerFaceColor','r');
f(6).Color(4) = 0.9;

% 2nd criterion: no timing requirement
SAT_nogo = [];Model_nogo = [];

for s = 1:length(SUB)
    SAT_nogo(:,s) = A.p_nogo_key{SUB(s)};
    Model_nogo(:,s) = A.ycdf_nogo_key{SUB(s)};
end

f(7) = plot(xplot,nanmean(SAT_nogo,2),'--','color',c(10,:),'markersize',mks+3,'LineWidth',1,'MarkerFaceColor','w');
f(7) = plot(xplot,nanmean(Model_nogo,2),'-','color',c(10,:),'markersize',mks,'LineWidth',2.5,'MarkerFaceColor','r');
f(7).Color(4) = 0.9;

% 3rd criterion: loose 10 ms; 5ms for each half of the circular stimulus
SAT_nogo = []; Model_nogo = [];
for s = 1:length(SUB)
    SAT_nogo(:,s) = A.p_nogo_loose10{SUB(s)};
    Model_nogo(:,s) = A.ycdf_nogo_loose10{SUB(s)};
end

f(3) = plot(xplot,nanmean(SAT_nogo,2),'--','color',c(90,:),'markersize',mks+3,'LineWidth',1,'MarkerFaceColor','w');
f(3) = plot(xplot,nanmean(Model_nogo,2),'-','color',c(90,:),'markersize',mks,'LineWidth',2.5,'MarkerFaceColor','r');
f(3).Color(4) = 0.9;

% 4th criterion: loose 20 ms; 10ms for each half of the circular stimulus
SAT_nogo = []; Model_nogo = [];
for s = 1:length(SUB)
    SAT_nogo(:,s) = A.p_nogo_loose20{SUB(s)};
    Model_nogo(:,s) = A.ycdf_nogo_loose20{SUB(s)};
end

f(4) = plot(xplot,nanmean(SAT_nogo,2),'--','color',c(60,:),'markersize',mks+3,'LineWidth',1,'MarkerFaceColor','w');
f(4) = plot(xplot,nanmean(Model_nogo,2),'-','color',c(60,:),'markersize',mks,'LineWidth',2.5,'MarkerFaceColor','r');
f(4).Color(4) = 0.9;

% 5th criterion: loose 40 ms; 20ms for each half of the circular stimulus
SAT_nogo = []; Model_nogo = [];
for s = 1:length(SUB)
    SAT_nogo(:,s) = A.p_nogo_loose40{SUB(s)};
    Model_nogo(:,s) = A.ycdf_nogo_loose40{SUB(s)};
end

f(5) = plot(xplot,nanmean(SAT_nogo,2),'--','color',c(30,:),'markersize',mks+3,'LineWidth',1,'MarkerFaceColor','w');
f(5) = plot(xplot,nanmean(Model_nogo,2),'-','color',c(30,:),'markersize',mks,'LineWidth',2.5,'MarkerFaceColor','r');
f(5).Color(4) = 0.9;

legend(f,["R-to-NR","30ms (Original)","35ms (Matched)","40ms","50ms","Simple Reaction"],'Location','Northwest','NumColumns',1,...
             'fontsize',14,'textcolor','k');
            legend('boxoff');
            
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
%      'PaperUnits', 'Inches', 'PaperSize', [5, 5])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal_Analysis\analysis\Figure');
% print(figure_SAT,'SAT_model_nolate_timewindow_adjust', '-depsc','-r600');


% %% extra: to see the upper asymtotic for each criterion
% % get the upper aympstotic for each criterion
% for s = 1:length(SUB)
%     beta_gono_nolate(s) = A.model_gono_nolate{SUB(s)}(4); 
%     beta_nogo_nolate(s) = A.model_nogo_nolate{SUB(s)}(4);
%     beta_gono_loose10(s) = A.model_gono_loose10{SUB(s)}(4); 
%     beta_nogo_loose10(s) = A.model_nogo_loose10{SUB(s)}(4);
%     beta_gono_loose20(s) = A.model_gono_loose20{SUB(s)}(4); 
%     beta_nogo_loose20(s) = A.model_nogo_loose20{SUB(s)}(4);
%     beta_gono_loose40(s) = A.model_gono_loose40{SUB(s)}(4); 
%     beta_nogo_loose40(s) = A.model_nogo_loose40{SUB(s)}(4);
% end
% 
% [h,p,ci,stats] = ttest(beta_nogo_nolate,beta_gono_nolate)
% [h,p,ci,stats] = ttest(beta_nogo_loose10,beta_gono_nolate)
% [h,p] = lillietest(beta_gono_nolate)
% [h,p] = lillietest(beta_nogo_loose10')
% [p,h,stats] = ranksum(beta_nogo_loose10,beta_gono_nolate)
% [h,p,ci,stats] = ttest(beta_nogo_loose20,beta_gono_nolate)
% 
% for s = 1:length(SUB)
%     mu_gono_nolate(s) = A.model_gono_nolate{SUB(s)}(1); 
%     mu_nogo_nolate(s) = A.model_nogo_nolate{SUB(s)}(1);
%     
%     mu_gono_loose20(s) = A.model_gono_loose20{SUB(s)}(1); 
%     mu_nogo_loose20(s) = A.model_nogo_loose20{SUB(s)}(1);
%     mu_gono_loose10(s) = A.model_gono_loose10{SUB(s)}(1); 
%     mu_nogo_loose10(s) = A.model_nogo_loose10{SUB(s)}(1);
%     mu_gono_button(s) = A.model_gono_button{SUB(s)}(1); 
%     mu_nogo_button(s) = A.model_nogo_button{SUB(s)}(1);
% end
% 
% [h,p,ci,stats] = ttest(mu_gono_nolate,mu_nogo_loose10)
% [h,p] = lillietest(mu_gono_nolate)
% [h,p] = lillietest(mu_nogo_loose20')