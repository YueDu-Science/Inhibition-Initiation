clear all;
close all;
clc;        
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));
addpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis\Auxiliary');
load Init_Inhb_Analysis.mat;  % Analysis
D = Analysis.D;
A = Analysis.A;

% set options for fmincon
options = optimoptions('fmincon','MaxIterations',1e5,'MaxFunctionEvaluations',1e4, ...
              'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-12,'Display','iter');
          
% need experimentl manipulated ST; this was fixed to all participants
% get CT that matched with experimental data
s = 1;
st = D.nogo{s}.t_prep; 

tic;
TAU = [0.2 0.4];
MU_RT = [0.47, 0.53];
SIGMA_RT = [0.005 0.055]; % too small values are not realistic for RT
TAU_RT = [0.001, 0.056]; % the mean of exponetial dist; the rate of expoential distribution is 1/tau_rt;
BETA_RT = [-1.43, -0.17];
SIGMA_TA = [0.018 0.05];
BETA_TA  = [-0.36 0.11];

% another parameters needed are mu_s, sigma_s, alpha_s, and beta_s
% import them from SAT for button data with t_prep.
for s = 1:size(D.sub_name,2)
    SAT_nogo(:,s) = D.model.nogo_button_sim{s}; 
end


lb=[0.1 0.4 0.005 0.001 -2 0 -2]; ub = [0.5 0.6 inf inf 2 inf 2];  
% 0.001 lower bound for tau_rt; because plotting with esGaussianpdf function messes up when it is less than 0.001
% theoretically, tau_rt > 0;
%lb(3) = 0.001; % + penalty = 0 in lik_nogo_rt_piecewise_mle for comparing parameter recovery

Aeq = []; beq = [];

%%
num_para = 7; % how many parameters the model has
num_start = 30; % 30 different starting positions for optimization
num_sim = 200; % how many sets of simulation data being fitted

% initialize outcome arrays
Xsim = nan(num_sim, num_para); % simulated para
Psim = nan(num_sim, num_para + 4); % all simulated value including alpha and beta
Xfit = nan(num_sim, num_para); % fitted para
Rsim = nan(num_sim, numel(st)); % simulated Response (0 or 1)
RTsim = nan(num_sim, numel(st)); % simulated RT
NegLL = nan(num_sim,1); % negative loglikelihood

for t = 1:num_sim
    % simulate data
    rng shuffle 
    
    % this step is to simulate whether participants respond or not respond
    % as well as if they respond, what is the time of response
    sub_v = randi(size(D.sub_name,2));
    SAT(1) = SAT_nogo(1, sub_v);
    SAT(2) = SAT_nogo(2, sub_v);
    SAT(3) = SAT_nogo(3, sub_v);
    SAT(4) = SAT_nogo(4, sub_v);

    tau = TAU(1) + (TAU(2)- TAU(1))*rand(1);    
    mu_rt = MU_RT(1) + (MU_RT(2)-MU_RT(1))*rand(1);
    sigma_rt = SIGMA_RT(1) + (SIGMA_RT(2)-SIGMA_RT(1))*rand(1);
    tau_rt = TAU_RT(1) + (TAU_RT(2)-TAU_RT(1))*rand(1);
    beta_rt = BETA_RT(1) + (BETA_RT(2)- BETA_RT(1))*rand(1);
    sigma_ta = SIGMA_TA(1) + (SIGMA_TA(2)-SIGMA_TA(1))*rand(1);
    beta_ta = BETA_TA(1) + (BETA_TA(2)- BETA_TA(1))*rand(1);
    
    [rsim, rtsim] = simulation_nogo_piecewise_mle(st, SAT, tau, mu_rt, sigma_rt, tau_rt, beta_rt,sigma_ta, beta_ta);
    Psim(t,:) = [tau, mu_rt, sigma_rt, tau_rt, beta_rt,sigma_ta, beta_ta, SAT(1), SAT(2),SAT(3),SAT(4)];
    Xsim(t,:) = [tau, mu_rt, sigma_rt, tau_rt, beta_rt,sigma_ta, beta_ta];
    Rsim(t,:) = rsim;
    RTsim(t,:) = rtsim;
    
    rt = rtsim;
    obFunc = @(x) lik_nogo_rt_piecewise_mle(rt, st, x(1), x(2), x(3), x(4), x(5), x(6),x(7));
    
    % fit model to simulated data
    for j = 1:num_start
        
        % random initial values
        rng shuffle
        tau0 = TAU(1) + (TAU(2)- TAU(1))*rand(1);
        mu_rt0 = MU_RT(1) + (MU_RT(2)-MU_RT(1))*rand(1);
        sigma_rt0 = SIGMA_RT(1) + (SIGMA_RT(2)-SIGMA_RT(1))*rand(1);
        tau_rt0 = TAU_RT(1) + (TAU_RT(2)-TAU_RT(1))*rand(1);
        beta_rt0 = BETA_RT(1) + (BETA_RT(2)- BETA_RT(1))*rand(1);
        sigma_ta0 = SIGMA_TA(1) + (SIGMA_TA(2)-SIGMA_TA(1))*rand(1);
        beta_ta0 = BETA_TA(1) + (BETA_TA(2)- BETA_TA(1))*rand(1);
        
        
        X0 = [tau0, mu_rt0, sigma_rt0, tau_rt0, beta_rt0, sigma_ta0, beta_ta0];
        init_negLL = lik_nogo_rt_piecewise_mle(rt, st, tau0, mu_rt0, sigma_rt0, tau_rt0, beta_rt0, sigma_ta0, beta_ta0);

        if ~isinf(init_negLL)
            [X(j,:), fval(j)] = fmincon(obFunc,X0,[],[],[],[],lb,ub,[],options);
        else
            X(j,:) = ones(1, num_para); fval(j) = 9999;
        end
    end
    Xf = X(find(fval == min(fval),1,'last'),:);
    Fval = fval(find(fval == min(fval),1,'last'));
    
    Xfit(t,:) = Xf;
    NegLL(t) = Fval;
end
toc;

% calculate mean squred error
% how fitted value deviated from the true value
sum((Xfit*1000-Xsim*1000).^2)/200

P_rec.Xfit = Xfit;
P_rec.NegLL = NegLL;
P_rec.Psim = Psim;
P_rec.Xsim = Xsim;
P_rec.Rsim = Rsim;
P_rec.RTsim = RTsim;
P_rec.st = st;
save('ParaRec_piecewise_mle.mat','P_rec');

% plot true value vs. estimated value
Symbols = {'\mu_{T}' '\mu_{0}' '\sigma_{rt}' '\delta' '\beta_{1}' '\sigma_{r}' '\beta_{2}'};
%symbols = Symbols(para_ind);
figure_para_recovery = figure('name','P_Recovery');

plot_row = ceil(sqrt(num_para));
plot_col = floor(sqrt(num_para));

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

% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 8],...
%      'PaperUnits', 'Inches', 'PaperSize', [10, 8])
% set(gcf,'renderer','Painters');
% cd('D:\Project\StopSignal_Analysis\analysis\Figure');
% print(figure_para_recovery,'Para_Recovery_WP', '-depsc','-r600');

%% FIGURE S9
% parameter recovery when there were outlier trials
num_para = 7; % how many parameters the model has
num_start = 30; % 30 different starting positions for optimization
num_sim = 200; % how many sets of simulation data being fitted

% initialize outcome arrays
Xsim = nan(num_sim, num_para); % simulated para
Psim = nan(num_sim, num_para + 4); % all simulated value including alpha and beta
Xfit = nan(num_sim, num_para); % fitted para
Rsim = nan(num_sim, numel(st)); % simulated Response (0 or 1)
RTsim = nan(num_sim, numel(st)); % simulated RT
NegLL = nan(num_sim,1); % negative loglikelihood

for t = 1:num_sim
    % simulate data
    rng shuffle 
    
    % this step is to simulate whether participants respond or not respond
    % as well as if they respond, what is the time of response
    sub_v = randi(size(D.sub_name,2));
    SAT(1) = SAT_nogo(1, sub_v);
    SAT(2) = SAT_nogo(2, sub_v);
    SAT(3) = SAT_nogo(3, sub_v);
    SAT(4) = SAT_nogo(4, sub_v);

    tau = TAU(1) + (TAU(2)- TAU(1))*rand(1);    
    mu_rt = MU_RT(1) + (MU_RT(2)-MU_RT(1))*rand(1);
    sigma_rt = SIGMA_RT(1) + (SIGMA_RT(2)-SIGMA_RT(1))*rand(1);
    tau_rt = TAU_RT(1) + (TAU_RT(2)-TAU_RT(1))*rand(1);
    beta_rt = BETA_RT(1) + (BETA_RT(2)- BETA_RT(1))*rand(1);
    sigma_ta = SIGMA_TA(1) + (SIGMA_TA(2)-SIGMA_TA(1))*rand(1);
    beta_ta = BETA_TA(1) + (BETA_TA(2)- BETA_TA(1))*rand(1);
    
    [rsim, rtsim] = simulation_nogo_piecewise_mle(st, SAT, tau, mu_rt, sigma_rt, tau_rt, beta_rt,sigma_ta, beta_ta);
    % randomly select a few small st and replace rtsim to ~0.5 to 0.55;
    ind_st = find(st < tau); ind_st1 = datasample(ind_st,5,'Replace',false);
    rtsim(ind_st1) = 0.45+(0.55-0.45).*rand(numel(ind_st1),1);
    
    Psim(t,:) = [tau, mu_rt, sigma_rt, tau_rt, beta_rt,sigma_ta, beta_ta, SAT(1), SAT(2),SAT(3),SAT(4)];
    Xsim(t,:) = [tau, mu_rt, sigma_rt, tau_rt, beta_rt,sigma_ta, beta_ta];
    Rsim(t,:) = rsim;
    RTsim(t,:) = rtsim;
    
    rt = rtsim;
    obFunc = @(x) lik_nogo_rt_piecewise_mle(rt, st, x(1), x(2), x(3), x(4), x(5), x(6),x(7));
    
    % fit model to simulated data
    for j = 1:num_start
        
        % random initial values
        rng shuffle
        tau0 = TAU(1) + (TAU(2)- TAU(1))*rand(1);
        mu_rt0 = MU_RT(1) + (MU_RT(2)-MU_RT(1))*rand(1);
        sigma_rt0 = SIGMA_RT(1) + (SIGMA_RT(2)-SIGMA_RT(1))*rand(1);
        tau_rt0 = TAU_RT(1) + (TAU_RT(2)-TAU_RT(1))*rand(1);
        beta_rt0 = BETA_RT(1) + (BETA_RT(2)- BETA_RT(1))*rand(1);
        sigma_ta0 = SIGMA_TA(1) + (SIGMA_TA(2)-SIGMA_TA(1))*rand(1);
        beta_ta0 = BETA_TA(1) + (BETA_TA(2)- BETA_TA(1))*rand(1);
        
        
        X0 = [tau0, mu_rt0, sigma_rt0, tau_rt0, beta_rt0, sigma_ta0, beta_ta0];
        init_negLL = lik_nogo_rt_piecewise_mle(rt, st, tau0, mu_rt0, sigma_rt0, tau_rt0, beta_rt0, sigma_ta0, beta_ta0);

        if ~isinf(init_negLL)
            [X(j,:), fval(j)] = fmincon(obFunc,X0,[],[],[],[],lb,ub,[],options);
        else
            X(j,:) = ones(1, num_para); fval(j) = 9999;
        end
    end
    Xf = X(find(fval == min(fval),1,'last'),:);
    Fval = fval(find(fval == min(fval),1,'last'));
    
    Xfit(t,:) = Xf;
    NegLL(t) = Fval;
end
toc;

% calculate mean squred error
% how fitted value deviated from the true value
sum((Xfit*1000-Xsim*1000).^2)/200

P_rec.Xfit = Xfit;
P_rec.NegLL = NegLL;
P_rec.Psim = Psim;
P_rec.Xsim = Xsim;
P_rec.Rsim = Rsim;
P_rec.RTsim = RTsim;
P_rec.st = st;
save('ParaRec_piecewise_mle.mat','P_rec');

% plot true value vs. estimated value
Symbols = {'\mu_{T}' '\mu_{0}' '\sigma_{rt}' '\delta' '\beta_{1}' '\sigma_{r}' '\beta_{2}'};
%symbols = Symbols(para_ind);
figure_para_recovery = figure('name','P_Recovery');

plot_row = ceil(sqrt(num_para));
plot_col = floor(sqrt(num_para));

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
print(figure_para_recovery,'Para_Recovery_WP_outlier', '-depsc','-r600');