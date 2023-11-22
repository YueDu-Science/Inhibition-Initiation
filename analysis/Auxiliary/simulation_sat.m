function [R] = simulation_sat(CT, PARA)

mu = PARA(1); sigma = PARA(2); alpha = PARA(3); beta = PARA(4);
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

Phi = [1 - normcdf(CT,mu,sigma)';
       normcdf(CT,mu,sigma)'];
   
P = [alpha beta];
presponse = P*Phi;

r_sample = rand(1,numel(CT));
R = 1 - sum(repmat(r_sample,size(presponse,1),1) > cumsum(presponse,1),1);