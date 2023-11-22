% This code reproduce figure 3
clear all;
clc;            
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));

% load data
load Init_Inhb_Analysis.mat;  % D
D = Analysis.D;
A = Analysis.A;

% load piecewice model fitting results
load ModelFitting_PW_MLE;
Xfit = M_PW.Xfit;
Xfit(16,:) = []; % sub 16 was excluded
INDX1 = M_PW.INDX1;
INDX2 = M_PW.INDX2;

st = M_PW.st;
xplot= D.xplot;
D.xplot = D.xplot;

% uses sub 10, consistent with Figure 2
s = 6;
data = D.nogo{s};
xfit = Xfit(s,:);
indx1 = INDX1{s}; % index does satisfy condition 1: see 4_Fitting_Piecewise_Model
indx2 = INDX2{s}; % index does satisfy condition 2: see 4_Fitting_Piecewise_Model
    
mks = 6; lw = 0.8;

%% Figure 3A
figure_RT_nogo = figure('name','RT_nogo');  
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
set(gca,'fontsize',10)
axis([0 0.52 0.3 1])
set(gca,'xTick',[0:0.1:0.5],'xTickLabel',[0:100:500],'FontSize',16, 'FontWeight','normal','FontName','Arial');
set(gca,'yTick',[0.3:0.1:1],'yTickLabel',[300:100:1000],'FontSize',16, 'FontWeight','normal','FontName','Arial');
xlabel('Intended RT (ms)','FontSize',16, 'FontWeight','normal');
ylabel('Time of Response (ms)','FontSize',16, 'FontWeight','normal');

hAx=gca;                    % create an axes
hAx.LineWidth=1.2; 

ind = find(data.button ~= -99); % trials that had a response no matter the time of that response
f = plot(data.t_prep(ind),data.button(ind),'ko','markersize',mks+1,'LineWidth',lw,'MarkerFaceColor','w');
f.Color(4) = 0.5;
g1 = plot(data.t_prep(indx1),data.button(indx1),'kx','markersize',mks+2,'LineWidth',lw,'MarkerFaceColor','w');
g2 = plot(data.t_prep(indx2),data.button(indx2),'k+','markersize',mks+2,'LineWidth',lw,'MarkerFaceColor','w');
        
%%% these lines plot trials without a response for reference
% ind = find(data.button == -99);
% data.button(ind) = 0.32;
% rng('default'); ss = rng; rng(ss);
% tmp1 = 0 + normrnd(0,0.01,numel(ind),1);
% scatter(data.t_prep(ind),data.button(ind)+tmp1,36,'r','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'MarkerFaceColor','r');

% the tolerance window for correct response (0.47 to 0.53 s)
y = 0:0.01:0.52;
plot(y, repmat(0.5,numel(y),1),'k-','LineWidth',1) % target line
plot(y, repmat(0.62,numel(y),1),'k--','LineWidth',1) % bottom of the monitor beyond which the stimulu is not visible
curve1 = repmat(0.47,1,length(y));
curve2 = repmat(0.53,1,length(y));
y2 = [y, fliplr(y)];
inBetween = [curve1, fliplr(curve2)];
h = fill(y2, inBetween, [0.5,0.5,0.5]);
h.EdgeColor = [0.5 0.5 0.5];
uistack(h,'bottom')

% y hat from the piecewise model
% fit_nogo_rt_2 vs fit_nogo_rt_1_total: one is the normal component of the
% exGaussian distribution, the other consistes of both normal and
% exponential components
st_sort = sort(st);
fit_nogo_rt_1 = (linspace(0.05,xfit(1),10) - xfit(1))*xfit(5) + xfit(2);
fit_nogo_rt_2 = (linspace(xfit(1),0.5,10) - xfit(1))*xfit(7) + xfit(2);
fit_nogo_rt_1_total = (linspace(0.05,xfit(1),10) - xfit(1))*xfit(5) + xfit(2) + xfit(4);

plot(linspace(0.05,xfit(1),10), fit_nogo_rt_1, '-','color',[30,144,255]/255, 'linewidth',2)
plot(linspace(xfit(1),0.5,10), fit_nogo_rt_2, '-','color',[30,144,255]/255, 'linewidth',2);
% the jonit point
plot(xfit(1), xfit(2),'o','color',[30,144,255]/255,'linewidth',2,'MarkerFaceColor', [30,144,255]/255,'markersize',mks+5);

text(0.3,0.8,'\tau','fontsize',18)

legend([h],{'60 ms window'},'Location','NorthWest','NumColumns',2,...
              'fontsize',10,'textcolor','k');
             legend('boxoff');
             
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 3],...
%      'PaperUnits', 'Inches', 'PaperSize', [5, 3])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal_Analysis\analysis\Figure');
% print(figure_RT_nogo,'RT_nogo_nolate_example', '-depsc','-r600');

%% Figure 3B
%DATA = readtable('tr-training.txt');
SUB = [1:15 17:size(D.sub_name,2)]; % sub 16 was excluded
mu_gono = [];
mu_nogo = [];
for s = 1:length(SUB) % load \mu estimated by SAT
    mu_gono(s) = A.model_gono_nolate{SUB(s)}(1); 
    mu_nogo(s) = A.model_nogo_nolate{SUB(s)}(1);
end

mks = 6; lw = 0.8;
figure_mu = figure('name','mu');  
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
%set(gca,'fontsize',10)
axis([0.2 0.4 0.2 0.4])
hAx=gca;                    % create an axes
hAx.LineWidth=1.2; 

%title(sub_tt,'FontSize',12, 'FontWeight','normal');
xlabel('\mu_{T} (ms)','FontSize',22, 'FontWeight','bold');
ylabel('\mu_{NR} (ms)','FontSize',22, 'FontWeight','bold');
pbaspect([1 1 1])
set(gca,'xTick',[0:0.05:0.5],'xTickLabel',[0:50:500],'FontSize',22, 'FontWeight','normal','FontName','Arial');
set(gca,'yTick',[0:0.05:0.5],'yTickLabel',[0:50:500],'FontSize',22, 'FontWeight','normal','FontName','Arial');

% plot diagnonal line
xl = get(gca, 'xlim');
plot(xl, xl,'k-','markersize',mks,'LineWidth',1.5,'MarkerFaceColor','k')

% individual data points
f1 = plot(Xfit(:,1), mu_gono,'ko','markersize',mks+8,'LineWidth',1.5,'MarkerFaceColor','w');

% [h,p,ci,stats] = ttest(mu_gono', Xfit(:,1))
% nanmean(mu_gono' - Xfit(:,1))
% [rho, p] = corr(mu_nogo', Xfit(:,1))

text(0.4,0.23,'\mu_{NR} - \mu_{T} = -6.57 ms','FontSize',22,'FontWeight','normal');
text(0.4,0.2,'95% CI: [-18.67, 5.54]ms','FontSize',22,'FontWeight','normal');

% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
%      'PaperUnits', 'Inches', 'PaperSize', [5, 5])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal_Analysis\analysis\Figure');
% print(figure_mu,'tau_mu_nolate_adjust', '-depsc','-r600');