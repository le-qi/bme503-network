
figure(1); clf 
subplot(2, 1, 1)
plot(t, r1)
xlabel ('Time (ms)', 'FontSize', 15, 'FontWeight', 'normal')
ylabel ('Firing rate (Hz)', 'FontSize', 15, 'FontWeight', 'normal')
title('Firing rate vs. time for input layer', 'FontSize', 15, 'FontWeight', 'bold')
subplot(2, 1, 2)
plot(t, lambda_stimulus)
xlabel ('Time (ms)', 'FontSize', 15, 'FontWeight', 'normal')
ylabel ('Incident light wavelength (nm)', 'FontSize', 15, 'FontWeight', 'normal')
title('Light input vs. time', 'FontSize', 15, 'FontWeight', 'bold')

figure(2); clf
for i = 1:length(color_names)
    subplot(length(color_names), 1, i)
    plot(t, r2(:, i), 'k', 'LineWidth', 3)
    title(sprintf(color_names{ceil(i)}), 'FontSize', 15, 'FontWeight', 'bold')
    if(i == 1) title(['Middle layer firing rates:' sprintf(color_names{ceil(i)})], 'FontSize', 15, 'FontWeight', 'bold')
    end
    ylim([0 1])
end

figure(3); clf; hold on
for i = 1:2:2*length(color_names)
subplot(length(color_names)+1, 1, ceil(i/2))
plot(t, r3(i, :), t, r3(i+1, :))
title(sprintf(color_names{ceil(i/2)}), 'FontSize', 15, 'FontWeight', 'bold')
    xlabel ('Time (ms)', 'FontSize', 12, 'FontWeight', 'normal')
    ylim([0 1])
if(i == 1) ylabel ('Firing rate', 'FontSize', 12, 'FontWeight', 'normal') 
    legend('Color detected', 'No color detected')
    title('Output layer firing rates: Violet', 'FontSize', 15, 'FontWeight', 'bold')
end
end

subplot(length(color_names)+1, 1, length(color_names)+1) 
plot(t, lambda_stimulus)
xlabel ('Time (ms)', 'FontSize', 15, 'FontWeight', 'normal')
ylabel ('Incident light wavelength (nm)', 'FontSize', 15, 'FontWeight', 'normal')
title('Light input vs. time', 'FontSize', 15, 'FontWeight', 'bold')