% Work script

%% Color setup

% Input layer settings:
l_low = 390; % nm
l_high = 700; % nm; range based on human visual spectrum

num_receptors = 4; % Number of input layer neurons
x = linspace(l_low, l_high, num_receptors+2); 
lambda_center = x(2:end-1); % Tuned lambdas 

% Colors classifier categories
color_names = {'Violet', 'Blue', 'Green', 'Yellow', 'Orange', 'Red'}; 
color_ranges = [[380 450]; 
    [450 495]
    [495 570]
    [570 590]
    [590 620]
    [620 750]]; 

color_center=bsxfun(@minus,color_ranges(:, 2), color_ranges(:,1)); 
color_center=color_ranges(:,1) + color_center/2; 

distance = abs(bsxfun(@minus, color_center, lambda_center));
% [color_category, dist to this neuron on lambda axis]
% Each row corresponds to a color category
% Each column corresponds to distance from one input neuron
tot = l_high - l_low; 
weights = (tot - distance)/tot-0.09; % Normalize distances from 0 to 1
% ***In the future: account for negative weights? add some linear -DC term***

% Use these as input weights from input layer to summing layer

%% Neuron setup 
% Parameters:
k = 0.1; 
theta = 0.2; % F(u=theta) = 0.5. Shifts F(u) to the right if bigger.
g = 0.5; % Delay variable gain  
weight = -0.8;


% Synaptic connection weights from layer 1 to 2. 
% Row: father neuron number in layer 1. 
% Column: receiving neuron number in layer 2. 

W12 = [weights(:, 1)';
    weights(:, 2)'; 
    weights(:, 3)';
    weights(:, 4)';]; 
W22 = 0; % No connections between neurons in second layer...though that's a good idea....
 
I = [NaN NaN]'; % Input current
tau_r = 1;   % Firing rate time constant
tau_a = 100; % Adaptation time constant

% Models:
F = @(u) (1+exp(-(u-theta)/k)).^-1;          % SCALAR FUNCTION. Sigmoidal current-to-firing rate: 
dr_dt = @(r, a) 1/tau_r * (-r+F(W*r-g*a+I)); % 1X2 OUT, 1X2 IN. Firing rate based on synaptic connection and DC current
% da_dt = @(r, a) 1/tau_a * (-a + r);          % 1X2 OUT, 1X2 IN. Adaptation rate (slows down as rate goes up)

dr_dt = @(r1,r2, W) 1/tau_r * (-r2+F(r1*W)); % Firing rate model without adaptation

%% Simulation
simtime = 60; %ms
dt = 0.1; %ms
t = 0:dt:simtime; 

sigma = 15; % Of Gaussian tuning curve
width_factor = sigma; % = 1/(sigma * sqrt(2pi))
max_Hz =  20; % arbitrarily 

noise_factor = 30; % lambda of poisson noise distribution
lambda_stimulus = 155*cos(2*pi*1/50*t) + 545 + poissrnd(20, 1, length(t)); 

lambda_stimulus = color_center(6)*ones(1, length(t))+poissrnd(20, 1, length(t)); % Dc input
lambda_stimulus = color_center(1).*(t<10) +...
    color_center(2).*(t>=10 & t<20) +...
    color_center(3).*(t>=20 & t<30) +...
    color_center(4).*(t>=30 & t<40) +...
    color_center(5).*(t>=40 & t<50) +...
    color_center(6).*(t>=50); 
lambda_stimulus = l_low + (l_high - l_low)/max(t) .* t + poissrnd(20, 1, length(t)); 


r1 = NaN(length(t), length(lambda_center)); % Store input layer rates 
r2 = zeros(length(t), length(color_names)); % Store middle layer rates
r3 = zeros(length(t), 2*length(color_names))'; % Store output layer rates
%r3 = [1 0.8]'; 
a3 = zeros(length(t), 2*length(color_names))'; 
%a3 = [0 0]'; 
W33 = zeros(2*length(color_names), 2*length(color_names)); 
weight = -0.8; 
for i = 1:2:2*length(color_names)
    W33(i, i+1) = weight;
    W33(i+1, i) = weight; % CPG wiring, mutual inhibition between pairs of neurons 
end
     
tau_a = 1; % Adaptation time constant
g = 0.1; % Delay variable gain  

dr_dt = @(r, W, I, a) 1/tau_r * (-r+F(W*r-g*a+I)); % 1X2 OUT, 1X2 IN. Firing rate based on synaptic connection and DC current
da_dt = @(r, a) 1/tau_a * (-a + r);          % 1X2 OUT, 1X2 IN. Adaptation rate (slows down as rate goes up)

for i = 1:length(t)
    % Input layer:
    for j = 1:length(lambda_center)
       r1(i, j) = f(lambda_stimulus(i), lambda_center(j), width_factor, max_Hz); 
    end
    
    % Current to middle layer: 
    cur12 = r1(i, :) * W12; 
    cur12 = cur12/100; % Scale current
    if(i ~= length(t)) % Avoid that last extra data point
        r2(i+1, :) = r2(i, :) + dt * dr_dt(r2(i, :), W22, cur12, 0); % Euler's forward method
        
        % Last layer; 6 CNG circuits
        cur23 = zeros(1, 2*length(color_names))'; % Current from middle to output layer
        DC_cur = 0.25; % Will need to tweak 
        for p = 1:2:length(cur23)
            cur23(p) =  r2(i, ceil(p/2)); 
            cur23(p+1) = DC_cur;  
        end
        %W33 = [0 -0.8; -0.8 0]; 
        %cur23 = [1 1]'; 
        r3(:, i+1) = r3(:, i) + dt * dr_dt(r3(:, i), W33, cur23, a3(:, i)); 
        a3(:, i+1) = a3(:, i) + dt * da_dt(r3(:, i), a3(:, i)); 
    end
end

%% Output

visualize_results