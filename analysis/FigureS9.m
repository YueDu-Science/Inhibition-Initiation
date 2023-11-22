clear all;
close all;
clc;        
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));
addpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis\Auxiliary');
load Init_Inhb_Clean.mat;  % load raw data
DATA = STOPSIG_CLEAN.EXP;
sub_name = unique(DATA.id);
xplot = 0.001:0.001:0.5;
x_size = 0.05; 

load Init_Inhb_Analysis.mat;  % Analysis
D = Analysis.D;
A = Analysis.A;
%%
tic;
MU = [0.247, 0.365];
SIGMA = [0.02, 0.09]; %only two participants had atypical sigma > 0.8 in r-to-nr condition
                      % so not inclue the max sigma 0.15 here. see line 57
ALPHA = [0.001, 0.099];
BETA = [0.84, 0.999];

num_para = 4; % how many parameters the model has
num_start = 30; % for the simple sat fit, not need of multipe initial value
num_sim = 200; % how many sets of simulation data being fitted

% initialize outcome arrays
Xsim = nan(num_sim, num_para); % simulated para
Xfit = nan(num_sim, num_para); % fitted para

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
    beta = BETA(1) + (BETA(2)- BETA(1))*rand(1);
    
    Xsim(t,:) = [mu, sigma, alpha, beta];
    
    Xsim(t,:) = (A.model_gono_nolate{ind_r})'; % use this line to simulate
    %instead of using random sigma, here we use subjects' sigma
    % the parameter recovery results is similar.
    
    [rsim] = simulation_sat(st_sim, Xsim(t,:));

    % fit model to simulated data
    for j = 1:num_start
        rng shuffle
        
        mu0 = MU(1) + (MU(2)-MU(1))*rand(1);
        sigma0 = SIGMA(1) + (SIGMA(2)-SIGMA(1))*rand(1);
        alpha0 = ALPHA(1) + (ALPHA(2)- ALPHA(1))*rand(1);
        beta0 = BETA(1) + (BETA(2)- BETA(1))*rand(1);

        X0 = [mu0, sigma0,alpha0, beta0];
        [para(j,:), ycdf(j,:), fval(j)] = fit_SAT(st_sim, rsim', X0, xplot);  
    end
    Xfit(t,:) = para(find(fval == min(fval),1,'last'),:);
end
toc;

% plot true value vs. estimated value
Symbols = {'\mu' '\sigma' '\alpha' '\beta'};
%symbols = Symbols(para_ind);
figure_para_recovery = figure('name','P_Recovery');

plot_row = ceil(sqrt(7));
plot_col = floor(sqrt(7));

if plot_row*plot_col < num_para
    plot_col = plot_col + 1;
end
mks = 4; lw = 1.2;
for i = 1:num_para
    subplot(plot_row,plot_col,i)
    set(gcf,'color','w');
    hold on
    set(gca,'TickDir','out');
    set(gca,'fontsize',10)
%     axis([0 0.52 0 0.5])
%     set(gca,'xTick',[0:0.1:0.5],'xTickLabel',[0:100:500],'FontSize',10, 'FontWeight','normal','FontName','Arial');
%     set(gca,'yTick',[0:0.1:0.6],'yTickLabel',[0:100:6800],'FontSize',10, 'FontWeight','normal','FontName','Arial');
    title(Symbols{i},'FontSize',12, 'FontWeight','normal');
    xlabel('true','FontSize',12, 'FontWeight','normal');
    ylabel('fitted','FontSize',12, 'FontWeight','normal');
    
    axis([min(Xfit(:,i)) max(Xfit(:,i)) min(Xfit(:,i)) max(Xfit(:,i))])
    
    plot(Xsim(:,i), Xfit(:,i),'o','Markersize', mks)
    xl = get(gca, 'xlim');
    plot(xl, xl, 'k--','linewidth',lw)
    pbaspect([1 1 1])

    rho = corr(Xsim(:,i),Xfit(:,i));
    ax = gca;
    txt = ['\rho =' num2str(rho)];
    text((ax.XLim(2) - ax.XLim(1))*0.1 + ax.XLim(1), (ax.YLim(2) - ax.YLim(1))*0.9 + ax.YLim(1), txt, 'FontSize', 12,'color','b');
end

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 8],...
     'PaperUnits', 'Inches', 'PaperSize', [10, 8])
set(gcf,'renderer','Painters');
cd('D:\Project\StopSignal_Analysis\analysis\Figure');
print(figure_para_recovery,'Para_Recovery_SAT', '-depsc','-r600');