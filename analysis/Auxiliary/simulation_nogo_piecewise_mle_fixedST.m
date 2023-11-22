function [R, RT] = simulation_nogo_piecewise_mle_fixedST(CT, R, tau, mu_rt, sigma_rt, tau_rt, beta_rt,sigma_ta, beta_ta)

% % step 1
% % T_s is the random time where when t >= T_s, a response can be issued with a probability of alpha;
% % when t < T_s, a response can be issued with a probability of beta;
% % T_s ~ N(mu_s, (sigma_s)^2);
% % This was estimated from the button data with t_prep where no time window criterion
% % in order to mimic participants' data, where most no responses happens
% when st is short

% mu_s = SAT(1); sigma_s = SAT(2); alpha_s = SAT(3); beta_s = SAT(4);
% R = zeros(numel(CT),1);
% sample_s = normrnd(mu_s, sigma_s, numel(R), 1);
% sample_lapse = rand(numel(R), 1);
% for i = 1:numel(CT)
%    if CT(i) >= sample_s(i) && sample_lapse(i) < beta_s
%        R(i) = 1;
%    end
%    if CT(i) < sample_s(i) && sample_lapse(i) < alpha_s
%        R(i) = 1;
%    end
% end

% step 2
% this part is used to simulate Response time for each trial.
% T_r is the random time where when t >= T_r, RT = CT, while when t < T_r, RT ~ N(mu_rt, (sigma_rt)^2);
% T_r ~ N(mu_r, (sigma_r)^2);

RT = nan(numel(CT),1);
ind_response = find(R == 1); % trials with a response

sample_tau = tau; %normrnd(mu_tau, sigma_tau, numel(RT), 1);
sample_rt_react = normrnd(mu_rt, sigma_rt, numel(RT), 1) + exprnd(tau_rt, numel(RT), 1); % generate RT
sample_rt_deadline = normrnd(mu_rt,sigma_ta,numel(RT), 1); 
%sample_rt_startle = normrnd(st,sigma_st,numel(RT), 1); 
for i = 1:numel(ind_response)
   if CT(ind_response(i)) >= sample_tau   
      RT(ind_response(i)) = sample_rt_deadline(ind_response(i)) + beta_ta*(CT(ind_response(i)) - sample_tau);
      
%       tmp = rand(1,1);
%       if tmp < 0.4 && CT(ind_response(i)) >= 0.4
%           RT(ind_response(i)) = RT(ind_response(i)) - CT(ind_response(i)) * 0.05;
%       end
   end
   
   
   if CT(ind_response(i)) < sample_tau
      RT(ind_response(i)) = sample_rt_react(ind_response(i)) + beta_rt*(CT(ind_response(i)) - sample_tau);
   end
end