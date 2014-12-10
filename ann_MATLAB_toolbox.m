% Make training dataset:
num_training_pts = 300;
l_model = linspace(l_low, l_high, num_training_pts);
in_model = [f(l_model, tune_l(1)); ...
    f(l_model, tune_l(2));...
    f(l_model, tune_l(3));...
    f(l_model, tune_l(4));];  % inputs x lambdas

out_model = [(l_model>380 & l_model<450);...
    (l_model>450 & l_model<495);...
    (l_model>495 & l_model<570);...
    (l_model>570 & l_model<590);...
    (l_model>590 & l_model<620);...
    (l_model>620 & l_model<750)];

%% Use MATLAB ANN toolbox
close all
net = fitnet(num_hidden_units); 

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

