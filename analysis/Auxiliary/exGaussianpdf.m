function f = exGaussianpdf(x,mu, sigma, tau)

% tau = 1/lamda;

part_1 = (1./(2*tau)).*exp((1./tau).*(mu + (sigma^2)./(2*tau) - x));

part_2 = erfc((mu + (sigma^2)./tau - x)./(sqrt(2)*sigma));

f = part_1.*part_2;