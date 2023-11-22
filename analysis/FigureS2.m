
% this code plot Phit for individuals and groups
%% Figue S2: non-switched trial
clear all;
clc;            
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));

mks = 8; lw = 1;

load Init_Inhb_Analysis.mat;  % Analysis
D = Analysis.D;
A = Analysis.A;
% for this analysis, A and D are the same because this analysis works on
% non-switched trials


Accuracy_gogo = [];
Accuracy_nono = [];
% individual data
for s = 1:size(D.nono,2)
    tmp_gogo = D.gogo{s}; Accuracy_gogo(s) = sum(tmp_gogo.correct_nolate) / numel(tmp_gogo.correct_nolate);
    tmp_nono = D.nono{s}; Accuracy_nono(s) = sum(tmp_nono.correct_nolate) / numel(tmp_nono.correct_nolate);
end

figure_accuracy = figure('name','Unchange_accuracy');  
hAx=gca;                    % create an axes
hAx.LineWidth=1.2;         
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
set(gca,'fontsize',10)
axis([0.5 2.5 0 1.2])
title('No-Switch Trials','FontSize',28, 'FontWeight','normal','FontName','Arial');
ylabel('Accuracy','FontSize',26, 'FontWeight','normal','FontName','Arial');
set(gca,'yTick',0:0.2:1,'yTickLabel',0:0.2:1,'FontSize',24, 'FontWeight','normal','FontName','Arial');

set(gca,'xTick',[1,2],'xTickLabel',{'Response','No Response'},'FontSize',24, 'FontWeight','normal','FontName','Arial');


bw = 0.6; lw = 1.2;
F1 = bar(1,nanmean(Accuracy_gogo),bw);
set(F1,'facecolor','w') % use color name
set(F1,'edgecolor','r','Linewidth',lw+1) % use color name

F2 = bar(2,nanmean(Accuracy_nono),bw);
set(F2,'facecolor','w') % use color name
set(F2,'edgecolor','b','Linewidth',lw+1) % use color name

rng('default'); s = rng; rng(s);
tmp1 = 1 + normrnd(0,0.12,1,numel(Accuracy_gogo));
plot(tmp1,Accuracy_gogo,'ro','MarkerSize',14,'LineWidth',1);

rng(s);
tmp2 = 2 + normrnd(0,0.12,1,numel(Accuracy_nono));
plot(tmp2,Accuracy_nono,'bo','MarkerSize',14,'LineWidth',1)

% identify the outlier
ind = find(Accuracy_nono < 0.6);
plot(tmp2(ind),Accuracy_nono(ind),'kx','MarkerSize',24,'LineWidth',1);
plot(tmp1(ind),Accuracy_gogo(ind),'kx','MarkerSize',24,'LineWidth',1);

% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
%      'PaperUnits', 'Inches', 'PaperSize', [5, 5])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal\analysis\Figure');
% print(figure_accuracy,'Unchange_Accuracy', '-depsc','-r600'); 