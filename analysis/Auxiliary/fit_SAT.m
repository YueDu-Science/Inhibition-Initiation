% This function fit the Speed-Accuracy tradeoff function

function [para, ycdf, fval] = fit_SAT(x, y, x0, xplot)
% x: preparation RT
% y: hit or not hit  0/1
% xplot: the time scale to generate predicted probability
% para: estimated para values
% ycdf: generated predicted probability using para.

options = optimoptions('fmincon','MaxIterations',1e5,...
          'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-6);

% set up model parameters with initial values
% 4 parameters in this model
mu = x0(1); % center 
sigma = x0(2); % var T~N(mu,sigma)
chance = x0(3); % lower asymt
asymt = x0(4); % upper asymt
lb=[0.05 0.019 0.001 0]; ub = [1 1 0.1 1]; % boundary of parameters
% we set sigma lower bound as 0.019 as if below that, one participant
% alwayrs has a unrealistic small sigma

% could also use a regularization term to constraint sigma
% it does not affect our result
alpha = 0; % regularization parameter
slope0 = .05; % slope prior

% input data
LL = @(params) -nansum(y.*log(params(3)+(params(4)-params(3))*normcdf(x,params(1),params(2)))...
        + (1-y).*log(1-(params(3)+(params(4)-params(3))*normcdf(x,params(1),params(2)))))...
        + alpha*(params(2)-slope0)^2; % control the
   % slope
    %  + alpha*sum(abs(params)); % lasso
    %+ alpha*(params*params'); % ridge
   % 
[para, fval] = fmincon(LL,[mu sigma chance asymt],[],[],[],[],lb,ub,[],options);
ycdf = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));