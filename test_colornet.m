% test color_net

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

figure(1); clf
subplot(2, 1, 1)
hold on
[R, G, B] = wavelength_to_RGB(tune_l(1))
plot(l_model, in_model(1, :), 'Color', [R G B] )
[R, G, B] = wavelength_to_RGB(tune_l(2))
plot(l_model, in_model(2, :), 'Color', [R G B] )
[R, G, B] = wavelength_to_RGB(tune_l(3))
plot(l_model, in_model(3, :), 'Color', [R G B] )
[R, G, B] = wavelength_to_RGB(tune_l(4))
plot(l_model, in_model(4, :), 'Color', [R G B] )
axis tight
xlabel '\lambda (nm)'
ylabel 'Firing rate (Hz)'
title 'Input Layer'


subplot(2, 1, 2);
hold on
predictions = color_net(in_model);

color_ranges = [[380 450];
    [450 495]
    [495 570]
    [570 590]
    [590 620]
    [620 750]];

color_center=bsxfun(@minus,color_ranges(:, 2), color_ranges(:,1));
color_center=color_ranges(:,1) + color_center/2;

for i = 1:length(color_center)
    [R, G, B] = wavelength_to_RGB(color_center(i));
    plot(l_model, predictions(i, :), 'Color', [R, G, B])
end


vline(380, 'k-', color_names{1})
for i = 1:length(color_ranges)-1
    vline(color_ranges(i, 2), 'k-', color_names{i+1})
end

axis tight
title 'Output Layer'
xlabel '\lambda (nm)'
ylabel 'Output'
hline(0.5, 'k:')
legend(['Number of hidden neurons:' sprintf(' %g', num_hidden_units)])

