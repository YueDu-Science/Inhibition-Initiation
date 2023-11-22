clear all;
clc;            
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));

mks = 6; lw = 0.8;

load Init_Inhb_Analysis.mat;  % Analysis
D = Analysis.D;
A = Analysis.A;
%head(DATA,3)  % check if data look right
xplot= D.xplot;

load Init_Inhb_Clean.mat;  % load raw data
DATA = STOPSIG_CLEAN.EXP;
%% Figure S3A
% whether the order of nr-to-r or r-to-nr affect our results
% half participants had nr-to-r first
% the other half had r-to-nr first
% randomly assigned
SUB = [1:15 17:size(D.sub_name,2)]; % subject 16 was excluded; see manuscript
sub_name = unique(DATA.id);
for s = 1:length(SUB)
    ind_sub = DATA.id == sub_name(s);
    data = DATA(ind_sub == 1,:);
    ind_blk = find(data.block_key == 2,1);
    if data.initial(ind_blk) == 0
        sub_order(s) = 0; %nr-to-r first
    elseif data.initial(ind_blk) == 1
        sub_order(s) = 1; % r-to-nr first
    end
end

for s = 1:length(SUB)
    mu_gono(s) = A.model_gono_nolate{SUB(s)}(1); % estimated value for model parameter
    mu_nogo(s) = A.model_nogo_nolate{SUB(s)}(1); % (1) is \mu
end

mu_gono0 = mu_gono(sub_order == 0);
mu_gono1 = mu_gono(sub_order == 1);

mu_nogo0 = mu_nogo(sub_order == 0);
mu_nogo1 = mu_nogo(sub_order == 1);

figure_rt = figure('name','rt');  
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
axis([0 6 0.15 0.45])
hAx=gca;                    % create an axes
hAx.LineWidth=1.2; 
ylabel('\mu (ms)','FontSize',22, 'FontWeight','bold');


set(gca,'xTick',[1.5 4.5],'xTickLabel',{'\mu_{NR}','\mu_{R}'},'FontSize',22, 'FontWeight','normal','FontName','Arial');
set(gca,'yTick',[0:0.05:0.5],'yTickLabel',[0:50:500],'FontSize',22, 'FontWeight','normal','FontName','Arial');
ax = gca;
ax.XAxis.TickLength = [0,0];

f1 = bar(1,nanmean(mu_gono0),'FaceColor','b','EdgeColor','b','Linewidth',lw);
f2 = bar(2,nanmean(mu_gono1),'FaceColor','w','EdgeColor','b','Linewidth',lw);
alpha(f1,.3)
alpha(f2,.3)

xx = 1+randn(100,1)*0.12;
for s = 1:size(mu_gono0,2)
    plot(xx(s),mu_gono0(s),'o','color','b','markerfacecolor','b','Markersize',mks,'linewidth',lw-0.4)
end

xx = 2+randn(100,1)*0.12;
for s = 1:size(mu_gono1,2)
    plot(xx(s),mu_gono1(s),'o','color','b','markerfacecolor','w','Markersize',mks,'linewidth',lw-0.4)
end

h = [1];
hE = errorbar(h',nanmean(mu_gono0),nanstd(mu_gono0)/sqrt(numel(mu_gono0)),...
    'k','MarkerSize',3);
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;

h = [2];
hE = errorbar(h',nanmean(mu_gono1),nanstd(mu_gono1)/sqrt(numel(mu_gono1)),...
    'k','MarkerSize',3);
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;

f1 = bar(4,nanmean(mu_nogo0),'FaceColor','r','EdgeColor','r','Linewidth',lw);
f2 = bar(5,nanmean(mu_nogo1),'FaceColor','w','EdgeColor','r','Linewidth',lw);
alpha(f1,.3)
alpha(f2,.3)

xx = 4+randn(100,1)*0.12;
for s = 1:size(mu_nogo0,2)
    plot(xx(s),mu_nogo0(s),'o','color','r','markerfacecolor','r','Markersize',mks,'linewidth',lw-0.4)
end

xx = 5+randn(100,1)*0.12;
for s = 1:size(mu_nogo1,2)
    plot(xx(s),mu_nogo1(s),'o','color','r','markerfacecolor','w','Markersize',mks,'linewidth',lw-0.4)
end

h = [4];
hE = errorbar(h',nanmean(mu_nogo0),nanstd(mu_nogo0)/sqrt(numel(mu_nogo0)),...
    'k','MarkerSize',3);
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;

h = [5];
hE = errorbar(h',nanmean(mu_nogo1),nanstd(mu_nogo1)/sqrt(numel(mu_nogo1)),...
    'k','MarkerSize',3);
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;

legend([f1,f2],{'R-to-NR first','NR-to-R first'},'Location','NorthWest','NumColumns',2,...
              'fontsize',10,'textcolor','k');
             legend('boxoff');
% to get stats and ci
% nanmean(mu_gono(sub_order==1)' - mu_nogo(sub_order==1)')
% [h,p,ci,stats] = ttest(mu_gono(sub_order==1),mu_nogo(sub_order==1))

% nanmean(mu_gono(sub_order==0)' - mu_nogo(sub_order==0)')
% [h,p,ci,stats] = ttest(mu_gono(sub_order==0),mu_nogo(sub_order==0))

% nanmean(mu_gono(sub_order==1)) - nanmean(mu_gono(sub_order==0))
% [h,p,ci,stats] = ttest2(mu_gono(sub_order==1),mu_gono(sub_order==0))

% nanmean(mu_nogo(sub_order==1)) - nanmean(mu_nogo(sub_order==0))
% [h,p,ci,stats] = ttest2(mu_nogo(sub_order==1),mu_nogo(sub_order==0))

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
     'PaperUnits', 'Inches', 'PaperSize', [5, 5])
set(gcf,'renderer','Painters');
cd('D:\Project\StopSignal_Analysis\analysis\Figure');
print(figure_rt,'order_effect', '-depsc','-r600');

%% Figure 3B
% whether differen expecttion build up for nr-to-r and r-to-nr
% plot speed for each three blocks of each condition

SUB = [1:15 17:size(D.sub_name,2)]; % subject 16 was excluded; see manuscript

for s = 1:length(SUB)
    mu_gono_early(s) = A.model_gono_nolate_early{SUB(s)}(1); % estimated value for model parameter
    mu_nogo_early(s) = A.model_nogo_nolate_early{SUB(s)}(1); % (1) is \mu
    mu_gono_late(s) = A.model_gono_nolate_late{SUB(s)}(1); % estimated value for model parameter
    mu_nogo_late(s) = A.model_nogo_nolate_late{SUB(s)}(1); % (1) is \mu
end

figure_exp = figure('name','exp');  
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
axis([0 6 0.15 0.45])
hAx=gca;                    % create an axes
hAx.LineWidth=1.2; 
ylabel('\mu (ms)','FontSize',22, 'FontWeight','bold');


set(gca,'xTick',[1.5 4.5],'xTickLabel',{'Early','Late'},'FontSize',22, 'FontWeight','normal','FontName','Arial');
set(gca,'yTick',[0:0.05:0.5],'yTickLabel',[0:50:500],'FontSize',22, 'FontWeight','normal','FontName','Arial');
ax = gca;
ax.XAxis.TickLength = [0,0];

f1 = bar(1,nanmean(mu_gono_early),'FaceColor','w','EdgeColor','b','Linewidth',lw);
f2 = bar(2,nanmean(mu_nogo_early),'FaceColor','w','EdgeColor','r','Linewidth',lw);
alpha(f1,.3)
alpha(f2,.3)

xx = 1+randn(100,1)*0.12;
for s = 1:size(mu_gono_early,2)
    plot(xx(s),mu_gono_early(s),'o','color','b','markerfacecolor','w','Markersize',mks,'linewidth',lw-0.4)
end

xx = 2+randn(100,1)*0.12;
for s = 1:size(mu_nogo_early,2)
    plot(xx(s),mu_nogo_early(s),'o','color','r','markerfacecolor','w','Markersize',mks,'linewidth',lw-0.4)
end

h = [1];
hE = errorbar(h',nanmean(mu_gono_early),nanstd(mu_gono_early)/sqrt(numel(mu_gono_early)),...
    'k','MarkerSize',3);
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;

h = [2];
hE = errorbar(h',nanmean(mu_nogo_early),nanstd(mu_nogo_early)/sqrt(numel(mu_nogo_early)),...
    'k','MarkerSize',3);
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;

f1 = bar(4,nanmean(mu_gono_late),'FaceColor','w','EdgeColor','b','Linewidth',lw);
f2 = bar(5,nanmean(mu_nogo_late),'FaceColor','w','EdgeColor','r','Linewidth',lw);
alpha(f1,.3)
alpha(f2,.3)

xx = 4+randn(100,1)*0.12;
for s = 1:size(mu_gono_late,2)
    plot(xx(s),mu_gono_late(s),'o','color','b','markerfacecolor','w','Markersize',mks,'linewidth',lw-0.4)
end

xx = 5+randn(100,1)*0.12;
for s = 1:size(mu_nogo_late,2)
    plot(xx(s),mu_nogo_late(s),'o','color','r','markerfacecolor','w','Markersize',mks,'linewidth',lw-0.4)
end

h = [4];
hE = errorbar(h',nanmean(mu_gono_late),nanstd(mu_gono_late)/sqrt(numel(mu_gono_late)),...
    'k','MarkerSize',3);
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;

h = [5];
hE = errorbar(h',nanmean(mu_nogo_late),nanstd(mu_nogo_late)/sqrt(numel(mu_nogo_late)),...
    'k','MarkerSize',3);
set(hE(1),'LineWidth',0.5,'color','k')
hE.CapSize = 0;

legend([f1,f2],{'\mu_{NR}','\mu_{R}'},'Location','NorthWest','NumColumns',2,...
              'fontsize',10,'textcolor','k');
             legend('boxoff');
% to get stats and ci
% nanmean(mu_gono_early' - mu_nogo_early')
% [h,p,ci,stats] = ttest(mu_gono_early' - mu_nogo_early')

% nanmean(mu_gono_late - mu_nogo_late)
% [h,p,ci,stats] = ttest(mu_gono_late - mu_nogo_late)

% nanmean(mu_gono_early - mu_gono_late)
% [h,p,ci,stats] = ttest(mu_gono_early,mu_gono_late)

% nanmean(mu_nogo_early - mu_nogo_late)
% [h,p,ci,stats] = ttest(mu_nogo_early - mu_nogo_late)

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
     'PaperUnits', 'Inches', 'PaperSize', [5, 5])
set(gcf,'renderer','Painters');
cd('D:\Project\StopSignal_Analysis\analysis\Figure');
print(figure_exp,'practice_effect', '-depsc','-r600');