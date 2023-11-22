% this code demonstrates that our estimation of speeds is not affected by
% trigger failure
% and the original method used in classic SRT task is affected by trigger
% failure
clear all;
clc;            
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));

%% Fig S6A
% simulation of two SAT with different upper asymptote

mks = 6; lw = 0.8;

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

xplot = 0.001:0.001:0.5;
c = copper(100); % the 37 row is red

para1 = [0.25, 0.05, 0.01, 1];
para2 = [0.25, 0.05, 0.01, 0.8];

ycdf1 = para1(4)*normcdf(xplot,para1(1),para1(2))+para1(3)*(1-normcdf(xplot,para1(1),para1(2)));
ycdf2 = para2(4)*normcdf(xplot,para2(1),para1(2))+para2(3)*(1-normcdf(xplot,para2(1),para2(2)));

f1 = plot(xplot,ycdf1,'-','color',c(70,:),'markersize',mks+3,'LineWidth',2,'MarkerFaceColor','w');
f2 = plot(xplot,ycdf2,'-','color',c(37,:),'markersize',mks,'LineWidth',2,'MarkerFaceColor','r');

legend([f1, f2],{'No Trigger Failure','Trigger Failure'},'Location','Northwest','NumColumns',1,...
             'fontsize',14,'textcolor','k');
            legend('boxoff');

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5],...
     'PaperUnits', 'Inches', 'PaperSize', [5, 5])
set(gcf,'renderer','Painters');
cd('D:\Project\StopSignal_Analysis\analysis\Figure');
print(figure_SAT,'Simulation_SAT', '-depsc','-r600');
%% Fig S6B
% plot estimated speeds by our model vs. original methods

% i.e., simulate dataset based on two sets of parameters; both with fixed
% but different upper asymptote
load Init_Inhb_Clean.mat;  % load raw data
DATA = STOPSIG_CLEAN.EXP;
sub_name = unique(DATA.id);
xplot = 0.001:0.001:0.5;
x_size = 0.05; 

tic;
MU = [0.247, 0.365];
SIGMA = [0.05, 0.05]; 
ALPHA = [0.01, 0.01];
BETA = [0.65, 0.99];

num_para = 4; % how many parameters the model has
num_start = 30; % for the simple sat fit, not need of multipe initial value
num_sim = 200; % how many sets of simulation data being fitted

% initialize outcome arrays
Xsim1 = nan(num_sim, num_para); % simulated para
Xsim2 = nan(num_sim, num_para); % all simulated value including alpha and beta
Xfit1 = nan(num_sim, num_para); % fitted para
Xfit2 = nan(num_sim, num_para); % fitted para
for t = 1:num_sim
     % simulate data
    data = [];
    rng shuffle 
    ind_r = datasample([1:15 17:36],1,'Replace',true); % randomly select a participant for each stimulation
    ind_sub = DATA.id == sub_name(ind_r);
    IND_R(t) = ind_r;
    
    data = DATA(ind_sub == 1,:);
    gono = data(data.initial == 1 & data.final == 0,:);
    gogo = data(data.initial == 1 & data.final == 1,:);
    gogo_resp = gogo(gogo.correct_choice == 1,:);
    st = gono.t_prep_nolate; % get thie participant's allowed rt
    
    st_sim = st;
    tmp = datasample(gogo_resp.t_choice,length(st),'Replace',false) - 0.5;
    ind = find(gono.correct_nolate == 1);
    st_sim(ind) = st_sim(ind) + tmp(ind);
        
    mu = MU(1) + (MU(2)-MU(1))*rand(1);
    sigma = SIGMA(1) + (SIGMA(2)-SIGMA(1))*rand(1);
    alpha = ALPHA(1) + (ALPHA(2)- ALPHA(1))*rand(1);
    beta1 = 1;
    beta2 = 0.8;
    
    Xsim1(t,:) = [mu, sigma, alpha, beta1];
    [rsim1] = simulation_sat(st_sim, Xsim1(t,:));

    Xsim2(t,:) = [mu, sigma, alpha, beta2];
    [rsim2] = simulation_sat(st_sim, Xsim2(t,:));
    
    [f1 N] = sliding_window(st_sim, rsim1,xplot,x_size);
    [f2 N] = sliding_window(st_sim, rsim2,xplot,x_size);
    
    % calculate using addictive method
    ind = find(f1>=0.5,1); speed1(t) = 0.5*ind/numel(xplot);
    ind = find(f2>=0.5,1); speed2(t) = 0.5*ind/numel(xplot);
    % fit model to simulated data
    for j = 1:num_start
        rng shuffle
        
        mu0 = MU(1) + (MU(2)-MU(1))*rand(1);
        sigma0 = SIGMA(1) + (SIGMA(2)-SIGMA(1))*rand(1);
        alpha0 = ALPHA(1) + (ALPHA(2)- ALPHA(1))*rand(1);
        beta0 = BETA(1) + (BETA(2)- BETA(1))*rand(1);

        X0 = [mu0, sigma0,alpha0, beta0];
        [para1(j,:), ycdf1(j,:), fval1(j)] = fit_SAT(st_sim, rsim1', X0, xplot);
        [para2(j,:), ycdf2(j,:), fval2(j)] = fit_SAT(st_sim, rsim2', X0, xplot);
    end
    Xfit1(t,:) = para1(find(fval1 == min(fval1),1,'last'),:);
    Xfit2(t,:) = para2(find(fval2 == min(fval2),1,'last'),:);
end

toc;

% plot results

figure_para_recovery = figure('name','P_Recovery');
mks = 4; lw = 1.2;

subplot(2,2,1)
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
set(gca,'fontsize',10)
title('Model Estimation','FontSize',12, 'FontWeight','normal');
xlabel('w/o Trigger Failure','FontSize',12, 'FontWeight','normal');
ylabel('w/ Trigger Failure','FontSize',12, 'FontWeight','normal');

axis([min(speed1) max(speed2) min(speed1) max(speed2)])

g1 = plot(Xfit1(:,1), Xfit2(:,1),'o','Markersize', mks)
xl = get(gca, 'xlim');
plot(xl, xl, 'k--','linewidth',lw)
pbaspect([1 1 1])    
    
rho = corr(Xfit1(:,1), Xfit2(:,1));
ax = gca;
txt = ['\rho =' num2str(rho)];
text((ax.XLim(2) - ax.XLim(1))*0.1 + ax.XLim(1), (ax.YLim(2) - ax.YLim(1))*0.9 + ax.YLim(1), txt, 'FontSize', 12,'color','b');

mse = sum((Xfit1(:,1)*1000 - Xfit2(:,1)*1000).^2)/num_sim;
mse = round(mse,2);
txt = ['MSE = ' num2str(mse)];
text((ax.XLim(2) - ax.XLim(1))*0.1 + ax.XLim(1), (ax.YLim(2) - ax.YLim(1))*0.8 + ax.YLim(1), txt, 'FontSize', 12,'color','b');

diff = round(nanmean(Xfit2(:,1)*1000 - Xfit1(:,1)*1000),2);
diff_se = round(nanstd(Xfit2(:,1)*1000 - Xfit1(:,1)*1000),2);
txt = ['Diff = ' num2str(diff) '(' num2str(diff_se) ') ms'];
text((ax.XLim(2) - ax.XLim(1))*0.1 + ax.XLim(1), (ax.YLim(2) - ax.YLim(1))*0.7 + ax.YLim(1), txt, 'FontSize', 12,'color','b');

subplot(2,2,2)
set(gcf,'color','w');
hold on
set(gca,'TickDir','out');
set(gca,'fontsize',10)
title('The Mean Method','FontSize',12, 'FontWeight','normal');
xlabel('w/o Trigger Failure','FontSize',12, 'FontWeight','normal');
ylabel('w/ Trigger Failure','FontSize',12, 'FontWeight','normal');
axis([min(speed1) max(speed2) min(speed1) max(speed2)])

rho = corr(speed1', speed2');
txt = ['\rho =' num2str(rho)];
text((ax.XLim(2) - ax.XLim(1))*0.1 + ax.XLim(1), (ax.YLim(2) - ax.YLim(1))*0.9 + ax.YLim(1), txt, 'FontSize', 12,'color','r');

g2 =  plot(speed1, speed2,'ro','Markersize', mks)
xl = get(gca, 'xlim');
plot(xl, xl, 'k--','linewidth',lw)
pbaspect([1 1 1])

mse = sum((speed1*1000 - speed2*1000).^2)/num_sim;
mse = round(mse, 2);
txt = ['MSE = ' num2str(mse)];
text((ax.XLim(2) - ax.XLim(1))*0.1 + ax.XLim(1), (ax.YLim(2) - ax.YLim(1))*0.8 + ax.YLim(1), txt, 'FontSize', 12,'color','r');

diff = round(nanmean(speed2*1000 - speed1*1000),2);
diff_se = round(nanstd(speed2*1000 - speed1*1000),2);
txt = ['Diff = ' num2str(diff) '(' num2str(diff_se) ') ms'];
text((ax.XLim(2) - ax.XLim(1))*0.1 + ax.XLim(1), (ax.YLim(2) - ax.YLim(1))*0.7 + ax.YLim(1), txt, 'FontSize', 12,'color','r');

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 8],...
     'PaperUnits', 'Inches', 'PaperSize', [10, 8])
set(gcf,'renderer','Painters');
cd('D:\Project\StopSignal_Analysis\analysis\Figure');
print(figure_para_recovery,'Para_Recovery_TriggerFailure', '-depsc','-r600');