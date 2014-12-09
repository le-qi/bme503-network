% This is a script that demonstrates the "wavelength_to_RGB" script, 
% where there is a conversion from a scalar
% to R, G, and B values. It visualizes the results in a handy figure.

% Michael Lee, Duke University, 12/9/2014

% Settings: 
lambda_low = 380; %nm; Lambda range. 
lambda_high = 780; %nm
n_points = 100; % Number of simulated lamdba points.


% Conversion:
l = linspace(lambda_low, lambda_high, n_points); 
[R, G, B] = wavelength_to_RGB(l); % Convert to RGB values

% Plotting: 
figure(1); clf
subplot(2, 1, 2)
hold on
plot(l, R, 'r', 'LineWidth', 2)
alpha(0.0001)
plot(l, G, 'g', 'LineWidth', 2)
plot(l, B, 'b', 'LineWidth', 2)
legend('R', 'G', 'B')
title('Lambda to RGB Encoding', 'FontSize', 15, 'FontWeight', 'bold')
xlabel('\lambda (nm)', 'FontSize', 12)

axis tight
ylim([-0.1 1.1])

xlow = l; % linspace(0, 1, length(l)); 
xhigh = [xlow(2:end) xlow(end)+1/length(l)]; 

figure(1); 
subplot(2, 1, 1)
hold on 
for i = 1:length(R)
    a = patch([xlow(i) xhigh(i) xhigh(i) xlow(i)], [0 0 1 1], [R(i), G(i), B(i)]); 
    set(a, 'LineStyle', 'none')
    set(gca, 'ytick', [])
end
axis tight 
title('Lambda to Color', 'FontSize', 15, 'FontWeight', 'bold')

xlabel('\lambda (nm)', 'FontSize', 12)
