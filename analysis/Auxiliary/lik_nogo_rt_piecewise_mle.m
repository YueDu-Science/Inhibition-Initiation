function negLL = lik_nogo_rt_piecewise_mle(RT, CT, tau, mu_rt, sigma_rt, tau_rt, beta_rt, sigma_ta, beta_ta)

penalty = 2000*(sigma_rt - 0.03)^2 + 2000*(mu_rt - 0.5)^2 + 2000*(tau_rt - 0.03)^2; % similar as accuracy analysis


negLL = -nansum(log(exGaussianpdf(RT,mu_rt + beta_rt*(CT - tau),sigma_rt, tau_rt).*(CT < (tau)) ...
    + normpdf(RT,mu_rt + beta_ta*(CT - tau),sigma_ta).*(CT >= (tau)))) + penalty;

% digitsOld = digits(50);
%negLL = double(vpa(negLL));