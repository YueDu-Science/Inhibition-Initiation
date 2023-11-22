%%
%%% this code plot Phit for individuals and groups
clear all;
clc;            
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));
mks = 6; lw = 0.8;

load Init_Inhb_Analysis.mat;  % Analysis
D = Analysis.D;
A = Analysis.A;
%head(DATA,3)  % check if data look right
xplot= D.xplot;

s = rng(01092022);
s = RandStream('mlfg6331_64'); 

sub = 6;
% using data A is ok too; A already incorporate the adjustment of RT for
% no-go trials, so no need the "datasample" step
% but since our original figure use this code and random seed, so I keep it
% here

% these steps are the same in 3_save_data_for_analysis.m when generating
% data A

% for GONO trials
gogo = D.gogo{sub}(D.gogo{sub}.correct_choice == 1,:); % get all gogo trials
tmp = datasample(s,gogo.t_choice,length(D.gono{sub}.t_prep_nolate),'Replace',false) - 0.5; % sample RT from the gogo trials
D.gono{sub}.sim_t_prep = D.gono{sub}.t_prep_nolate;
ind = find(D.gono{sub}.correct_nolate == 1);
D.gono{sub}.sim_t_prep(ind) = D.gono{sub}.sim_t_prep(ind) + tmp(ind);
                
% for NOGO trials
tmp = datasample(s,gogo.t_choice,length(D.gono{sub}.t_prep_nolate),'Replace',false) - 0.5;
D.nogo{sub}.sim_t_prep = D.nogo{sub}.t_prep_nolate;
ind = find(D.nogo{sub}.correct_nolate == 0);
D.nogo{sub}.sim_t_prep(ind) = D.nogo{sub}.sim_t_prep(ind) + tmp(ind);

gono = D.gono{sub};
nogo = D.nogo{sub};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2B
% plot data
% NOGO trials
figure_PT_nogo = figure('name','PT_nogo');  
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
set(gca,'fontsize',10)
axis([0 0.5 -0.4 1.4])
title('NR-to-R Condition','FontSize',20, 'FontWeight','normal');
xlabel('Allowed RT (ms)','FontSize',20, 'FontWeight','normal');
%ylabel('Response','FontSize',16, 'FontWeight','normal');
set(gca,'xTick',[0:0.1:0.5],'xTickLabel',[0:100:500],'FontSize',20, 'FontWeight','normal','FontName','Arial');
set(gca,'yTick',[])
hAx=gca;                    % create an axes
hAx.LineWidth=1.2; 

% plot ST line
% these are the original allowed RT (# of 28 because of 60Hz monitor refresh rate)
for i = 1:length(D.xplot)
    f3 = plot([D.xplot(i) D.xplot(i)],[-0.4 1.4],'-','color',[0.8 0.8 0.8],'LineWidth',1);
    f3.Color(4) = 0.5;
end

rng('default'); s = rng; rng(s); % using the same random seed to recaptulate figures in maniscript
tmp1 = 0 + normrnd(0,0.05,numel(nogo.correct_nolate),1);
trial = find(nogo.correct_nolate == 1);
f1 = scatter(nogo.sim_t_prep(trial),nogo.correct_nolate(trial)+tmp1(trial),36,'r','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);

trial = find(nogo.correct_nolate == 0);
f2 = scatter(nogo.sim_t_prep(trial),nogo.correct_nolate(trial)+tmp1(trial),36,'r','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);

% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 3],...
%      'PaperUnits', 'Inches', 'PaperSize', [5, 3])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal\analysis\Figure');
% print(figure_PT_nogo,'PT_nogo_nolate_adjust', '-depsc','-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GONO trials
figure_PT_gono = figure('name','PT_gono');  
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
set(gca,'fontsize',10)
axis([0 0.5 -0.4 1.4])
title('R-to-NR Condition','FontSize',20, 'FontWeight','normal');
xlabel('Allowed RT(ms)','FontSize',20, 'FontWeight','normal');
%ylabel('Response','FontSize',16, 'FontWeight','normal');
set(gca,'xTick',[0:0.1:0.5],'xTickLabel',[0:100:500],'FontSize',20, 'FontWeight','normal','FontName','Arial');
set(gca,'yTick',[])
hAx=gca;                    % create an axes
hAx.LineWidth=1.2; 

% plot ST line
for i = 1:length(D.xplot)
    f3 = plot([D.xplot(i) D.xplot(i)],[-0.4 1.4],'-','color',[0.8 0.8 0.8],'LineWidth',1);
    f3.Color(4) = 0.5;
end

rng('default'); s = rng; rng(s);
tmp1 = 0 + normrnd(0,0.05,numel(gono.correct_nolate),1);
trial = find(gono.correct_nolate == 1);
f1 = scatter(gono.sim_t_prep(trial),gono.correct_nolate(trial) - 1 +tmp1(trial),36,'b','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);

trial = find(gono.correct_nolate == 0);
f2 = scatter(gono.sim_t_prep(trial),gono.correct_nolate(trial) + 1 +tmp1(trial),36,'b','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
          
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 3],...
%      'PaperUnits', 'Inches', 'PaperSize', [5, 3])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal\analysis\Figure');
% print(figure_PT_gono,'PT_gono_nolate_adjust', '-depsc','-r600');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2C
%%% this part plot Phit for individuals and groups
SAT_gono = A.p_gono_nolate{sub};   % speed-accuracy tradeoff
SAT_nogo = A.p_nogo_nolate{sub};

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

f2 = plot(xplot,SAT_nogo,'r--','markersize',mks+6,'LineWidth',1,'MarkerFaceColor','w','LineStyle','--');
f1 = plot(xplot,SAT_gono,'b--','markersize',mks+6,'LineWidth',1,'MarkerFaceColor','w','LineStyle','--');

legend([f1,f2],{'R-to-NR' ,'NR-to-R'},'Location','Northwest','NumColumns',1,...
             'fontsize',14,'textcolor','k');
            legend('boxoff');
            
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
%      'PaperUnits', 'Inches', 'PaperSize', [5, 5])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal\analysis\Figure');
% print(figure_SAT,'SAT_Individual_nolate_adjust', '-depsc','-r600');
% cd('D:\Project\StopSignal\analysis\Manuscript\NewAnalysis_NoLate');


