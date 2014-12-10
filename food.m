%This function is to make food with a random color at a 
%random location. It converts a random wavelength to the 
%corresponding RGB value and then creates a food particle
%with the color associated with that RGB value

food_lambda = (750-380)*rand() + 380; % food color
[R, G, B] = wavelength_to_RGB(food_lambda); % for the purposes of drawing
f_center = tgt;
f_width = 5;
food = rectangle('Position',[f_center(1)-f_width(1)/2,f_center(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5],'EdgeColor',[R, G, B], 'LineWidth', 4);
set(food,'Position',[f_center(1)-f_width(1)/2,f_center(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5]);

