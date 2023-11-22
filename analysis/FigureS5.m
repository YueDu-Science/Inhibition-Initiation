% Figure 4: RT distribution for gogo trials
%% Figure S4A
clear all;
clc;            
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));
mks = 6; lw = 0.8;

load Init_Inhb_Analysis.mat;  % Analysis
D = Analysis.D;
A = Analysis.A;
%head(DATA,3)  % check if data look right
xplot= D.xplot;
mks = 8; lw = 1;

sub_name = D.sub_name;
for s = 1:length(sub_name)
        S.gogo{s} = D.gogo{s}(D.gogo{s}.correct_choice == 1,:);
end

figure_RTdis = figure('name','RTdis');  
hAx=gca;                    % create an axes
hAx.LineWidth=1.2;            
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
set(gca,'fontsize',10)
ylabel('Frequency','FontSize',27, 'FontWeight','normal','FontName','Arial');
xlabel('Time of Response (ms)','FontSize',27, 'FontWeight','normal','FontName','Arial');
set(gca,'yTick',0:0.02:1,'yTickLabel',0:0.02:1,'FontSize',27, 'FontWeight','normal','FontName','Arial');
set(gca,'ytick',[]);
set(gca,'xTick',400:100:600,'xTickLabel',400:100:600,'FontSize',27, 'FontWeight','normal','FontName','Arial');
axis([350 650 0 0.14])

% examplar participant s= 10 to match with other figures
for s = 10:10 
   histogram(S.gogo{s}.t_choice*1000,300:8:640,'facecolor','r','facealpha',0.5,'edgecolor','none','normalization','probability');
   plot([nanmean(S.gogo{s}.t_choice*1000) nanmean(S.gogo{s}.t_choice*1000)], [0, 0.14],'-')
end

% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
%      'PaperUnits', 'Inches', 'PaperSize', [5, 5])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal_Analysis\analysis\Figure');
% print(figure_RTdis,'RTDis_gogo', '-depsc','-r600'); 

%% Figure S4B
% RT confidence interval for each individual
z = 1.96; % 95% CI
RT_mu = [];
for s = 1:length(sub_name)
    n = length(S.gogo{s}.t_choice);
    mu_t = nanmean(S.gogo{s}.t_choice)*1000;
    se = nanstd(S.gogo{s}.t_choice)*1000/sqrt(n);
    ci_u = mu_t + z*se;
    ci_b = mu_t - z*se;
    RT_mu(s,:) = [mu_t ci_b ci_u];
end

%% group mean
figure_RT_mu = figure('name','RTdis_mu');  
hAx=gca;                    % create an axes
hAx.LineWidth=1.2;          
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
set(gca,'fontsize',10)
ylabel('RT (ms)','FontSize',18, 'FontWeight','normal','FontName','Arial');
xlabel('Participants','FontSize',18, 'FontWeight','normal','FontName','Arial');
set(gca,'yTick',460:20:540,'yTickLabel',460:20:540,'FontSize',18, 'FontWeight','normal','FontName','Arial');
set(gca, 'YDir','reverse')
set(gca,'xtick',[]);
axis([0 36 455 545])

plot([0 36],[500 500],'w-','Linewidth',lw+0.5)
x = 1:length(sub_name);
for s = 1:length(sub_name)
    f1 = plot([x(s) x(s)], RT_mu(s,[2 3]),'r-','Linewidth',lw,'Markersize',mks-2);
    f2 = plot(x(s), RT_mu(s,1),'ro','Linewidth',lw,'Markersize',mks-2,'Markerfacecolor','r');
end

% timing tolerance (470 to 530 ms)
y = 0.2:36;
curve1 = repmat(470,1,length(y));
curve2 = repmat(530,1,length(y));
y2 = [y, fliplr(y)];
inBetween = [curve1, fliplr(curve2)];
h = fill(y2, inBetween, [0.85,0.85,0.85]);
h.EdgeColor = [0.85 0.85 0.85];
uistack(h,'bottom')

legend([f1,f2],{'95% CI' ,'Mean'},'Location','North','NumColumns',2,...
             'fontsize',20,'textcolor','k');
            legend('boxoff');

% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 3],...
%      'PaperUnits', 'Inches', 'PaperSize', [5, 3])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal_Analysis\analysis\Figure');
% print(figure_RT_mu,'RTDis_mu', '-depsc','-r600'); 

