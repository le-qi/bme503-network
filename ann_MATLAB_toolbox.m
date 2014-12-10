% use MATLAB's neuron toolbox
num_inputs = 4; % 4 tuned neurons
num_hidden_units = 5; % minimum with good performance
num_outputs = 6;

% Input layer characteristics:
color_names = {'Violet', 'Blue', 'Green', 'Yellow', 'Orange', 'Red'};
color_ranges = [[380 450]; [450 495];
    [495 570]; [570 590];
    [590 620]; [620 750]];
num_receptors = num_inputs; % Number of input layer neurons
l_low = 390; l_high = 700; % nm; range based on human visual spectrum
x = linspace(l_low, l_high, num_receptors+2);
tune_l = x(2:end-1); % Tuned lambdas
sigma = 30; % Of Gaussian tuning curve
width_factor = sigma; % = 1/(sigma * sqrt(2pi))
max_Hz =  20; % arbitrarily

% Make training dataset:
num_training_pts = 300;
l_model = linspace(l_low, l_high, num_training_pts);
in_model = [f(l_model, tune_l(1), width_factor, max_Hz); ...
    f(l_model, tune_l(2), width_factor, max_Hz);...
    f(l_model, tune_l(3), width_factor, max_Hz);...
    f(l_model, tune_l(4), width_factor, max_Hz);];  % inputs x lambdas

out_model = [(l_model>380 & l_model<450);...
    (l_model>450 & l_model<495);...
    (l_model>495 & l_model<570);...
    (l_model>570 & l_model<590);...
    (l_model>590 & l_model<620);...
    (l_model>620 & l_model<750)];


%% MATLAB functionality
close all
net = fitnet(num_hidden_units); 
view(net)

x = in_model; 
t = out_model; 
[net, tr] = train(net, x, t); 
nntraintool
plotperform(tr)

testX = x(:,tr.testInd);
testT = t(:,tr.testInd);
testY = net(testX);
perf = mse(net,testT,testY)

y = net(x);
plotregression(t,y)

e = t - y;
%ploterrhist(e)

genFunction(net, 'color_net'); 

test_colornet
