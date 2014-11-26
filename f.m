function r = f( lambda, lambda_pref, width_factor, max_rate)
% Firing rate from incident lambda based on preferred lambda
% width_factor: width of tuning curve
% max_rate: max firing rate

% Gaussian function
a = max_rate; 
b = lambda_pref; 
c = width_factor; 
d = 0; % asymptotic limit
x = lambda; 

r = a*exp(-(x-b).^2/(2*c.^2))+d; 


end

