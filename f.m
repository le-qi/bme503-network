function r = f( lambda, lambda_pref)
% Firing rate from incident lambda based on preferred lambda
% width_factor: width of tuning curve
% max_rate: max firing rate
global width_factor
global max_rate
% sigma = 30; % Of Gaussian tuning curve
% width_factor = sigma; % = 1/(sigma * sqrt(2pi))
% max_rate =  20; % Hz, arbitrarily

% Gaussian function
a = max_rate; 
b = lambda_pref; 
c = width_factor; 
d = 0; % asymptotic limit
x = lambda; 

r = a*exp(-(x-b).^2/(2*c.^2))+d; 


end

% test