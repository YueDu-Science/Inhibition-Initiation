clear all;
clc;        
addpath(genpath('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis'));

load Init_Inhb_Analysis.mat;  
D = Analysis.D; % D: This is for extracting CT from experimental data.
AA = Analysis.A;
% extract data for each individual
for s = 1:size(D.sub_name,2)
    SAT_nogo(:,s) = D.p_nogo_button{s};
    Model_nogo(:,s) = D.ycdf_nogo_button{s};
    Para_nogo(:,s) = D.model.nogo_button_sim{s};
end

% set options for fmincon
options = optimoptions('fmincon','MaxIterations',1e6,'MaxFunctionEvaluations',1e4, ...
              'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-12,'Display','iter');
          
% set parameter value range
TAU = [0.2 0.42];  % \mu(T) in the paper: the transition time when participant switch response mode
MU_RT = [0.47, 0.53]; % \mu when t = tau
SIGMA_RT = [0.005 0.055]; % \sigma RT when t < tau
TAU_RT = [0.001, 0.056]; % the mean of exponetial dist; the rate of expoential distribution is 1/tau_rt;
BETA_RT = [-2, -0.17]; % \beta(1) in the paper: linear rate how RT changes with allowed RT when simple responding (when t < tau)
SIGMA_TA = [0.01 0.05]; % \sigma(r) in the paper: when t > tau
BETA_TA  = [-0.5 0.2]; % \beta(2) in the paper: linear rate how RT changes with allowed RT (when t > tau)

% constraints on parameters
lb=[0.1 0.4 0.005 0.001 -inf 0 -inf]; ub = [0.5 0.6 inf inf inf inf inf];  

Aeq = []; beq = [];
A = []; b = [];

num_para = 7;
num_start = 150; % use # of starting value to avoid local solutions
X = nan(num_start, num_para);
fval = nan(num_start,1);
XXX = {};

col = size(D.nogo{1},1);
Xfit= nan(length(D.sub_name),num_para);

for s = 1:length(D.sub_name)

    data = [];
    data = D.nogo{s};
    col = size(data,1);
    
    data.button(data.button == -99) = nan;
    st = data.t_prep;
    r = data.button_choice;
    rt = data.button; % 0.5 is where the target line locates
    
    rt1 = rt - 0.5 + st;
    
    % exclude bad trials
    % condition 1: respone should reflects responding after seeing the stimulus changes instead of the
    % attemption to respond from the beginning of the trial (i.e., I will respond this trial no matter what)
    
    % remove trials that participants press a button even when it is
    % impossible.
    
    % this was determined in this way: if the stimulus change time is short
    % (shorter than the cutoff at which the SAT starts to rise in figure 2)
    % and still a response was made within 200 ms, which is much shorter
    % than the time needed
    indx1 = find((st <= AA.t_min(s)/1000 & rt1 <= 0.2));
    
    INDX1{s} = indx1;
    IND1(s) = numel(indx1);
    RT(s) = numel(~isnan(rt));
    RT1(s) = numel(find(st <= AA.t_min(s)/1000));
   
    % condition 2: not responding too late 
    % for some trials, participants lack of attention
    % Responding too late even when allowed RT was long
    % e.g. Responding after the stimulus was not on the screen even if the stimulus appears very early
    % this is determined in this way: if the stimlus change is long
    % (longer than the cutoff at which the SAT approached the upper aymptote in figure 2)
    % and participant did not make a response until the stimulus goes off
    % the screen (i.e., the time of response > 620 ms)
    indx2 = [];
    indx2 = find((st >= AA.t_max(s)/1000 & rt > 0.62));
    INDX2{s} = indx2;
    IND2(s) = numel(indx2);
    RT2(s) = numel(find(st >= AA.t_max(s)/1000));

    indx = [indx1; indx2];
     
    rt(indx) = [];
    st(indx) = [];

    num_sample = sum(~isnan(rt));
    
    obFunc = @(x) lik_nogo_rt_piecewise_mle(rt, st, x(1), x(2), x(3), x(4), x(5), x(6), x(7));
    for j = 1:num_start
        
        rng shuffle
        tau0 = TAU(1) + (TAU(2)- TAU(1))*rand(1);
        mu_rt0 = MU_RT(1) + (MU_RT(2)-MU_RT(1))*rand(1);
        sigma_rt0 = SIGMA_RT(1) + (SIGMA_RT(2)-SIGMA_RT(1))*rand(1);
        tau_rt0 = TAU_RT(1) + (TAU_RT(2)-TAU_RT(1))*rand(1);
        beta_rt0 = BETA_RT(1) + (BETA_RT(2)- BETA_RT(1))*rand(1);
        sigma_ta0 = SIGMA_TA(1) + (SIGMA_TA(2)-SIGMA_TA(1))*rand(1);
        beta_ta0 = BETA_TA(1) + (BETA_TA(2)- BETA_TA(1))*rand(1);
        
        X0 = [tau0, mu_rt0, sigma_rt0, tau_rt0, beta_rt0, sigma_ta0, beta_ta0]
        init_negLL = lik_nogo_rt_piecewise_mle(rt, st, tau0, mu_rt0, sigma_rt0, tau_rt0, beta_rt0, sigma_ta0, beta_ta0)

        if ~isinf(init_negLL) && ~isnan(init_negLL)
            [X(j,:), fval(j)] = fmincon(obFunc,X0,A,b,Aeq,beq,lb,ub,[],options);
        else
            X(j,:) = ones(1, num_para); fval(j) = 9999;
        end
    end
    Xf = X(find(fval == min(fval)),:);
    Fval = fval(find(fval == min(fval)));
    count = 0;
    if size(Xf,1) > 1
        XXX{s} = Xf;
        ind = find(Xf(:,1) == min(Xf(:,1)), 1, 'first');
        Xf = Xf(ind,:);
        Fval = Fval(ind);
        count = count + 1;
    end
    Xfit(s,:) = Xf;
    NEGLL(s) = Fval;
    NUM_SAMPLE(s) = num_sample;
end

BIC = length(X0) * log(NUM_SAMPLE) + 2*NEGLL;
AIC = 2*length(X0)+ 2*NEGLL;

M_PW.Xfit = Xfit;
M_PW.st = st;
M_PW.AIC = AIC;
M_PW.BIC = BIC;
M_PW.NEGLL = NEGLL;
M_PW.INDX1 = INDX1;
M_PW.INDX2 = INDX2;
M_PW.num_plan = IND1;
M_PW.num_lapse = IND2;
M_PW.num_total = RT;
M_PW.num_total_plan = RT1;
M_PW.num_total_lapse = RT2;

cd('D:\Project\Publication_Data_Code\Initiation-versus-Inhibition\analysis');
datafname = ['ModelFitting_PW_MLE.mat'];
save(datafname, 'M_PW');

