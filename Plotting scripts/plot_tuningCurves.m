% Plot tuning curves


lambda_model = linspace(l_low, l_high, 1e4);
for j = 1:length(lambda_center)
    r1_plot{j} = f(lambda_model, lambda_center(j), width_factor, max_Hz);
end

close all;
figure(1); clf
hold on
for j = 1:length(lambda_center)
    [R, G, B] = wavelength_to_RGB(lambda_center(j))
    plot(lambda_model, r1_plot{j}, 'Color', [R, G, B], 'LineWidth', 1)
end

title('Input Neurons Tuning Curves', 'FontSize', 14, 'FontWeight', 'bold')
axis tight
xlabel('\lambda (nm)', 'FontSize', 12)
ylabel('Firing rate (Hz)', 'FontSize', 12)
