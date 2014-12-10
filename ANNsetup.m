%% Define parameters 
clear 
close all
% Tuned neuron characteristics
global width_factor
global max_rate
global tune_l
width_factor = 30; % sigma in Guassian distribution 
max_rate = 20; %Hz

% ANN settings 
num_inputs = 4; % 4 tuned neurons
num_hidden_units = 10; % minimum with good performance
num_outputs = 6;

% Input layer characteristics:
num_receptors = num_inputs; % Number of input layer neurons
l_low = 380; l_high = 750; % nm; range based on human visual spectrum
x = linspace(l_low, l_high, num_receptors+2);
tune_l = x(2:end-1); % Tuned lambdas for receptors

% Output layer characteristics - just a reference, already hardcoded into
% ann_MATLAB_toolbox.m
color_names = {'Violet', 'Blue', 'Green', 'Yellow', 'Orange', 'Red'}; % Categories, for reference
color_ranges = [[380 450]; [450 495]; % Category definitions
    [495 570]; [570 590];
    [590 620]; [620 750]];


%% Update ANN model
ann_MATLAB_toolbox % also updates function color_net.m
test_colornet % Plots results
