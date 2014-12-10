%This function is to make food with a random color at a 
%random location. It converts a random wabelength to the 
%corresponding RGB values and then creates a food particle
%with the color associated with that RGB value

function [] = food()

lambda = (750-380)*rand() + 380;

[R,G,B] = wavelength_to_RGB(lambda);

tgt = [150*rand(1,2)-[75. 75.]];
f_center = tgt;
f_width = 5;
f = rectangle('Position',[f_center(1)-f_width(1)/2,f_center(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5],'EdgeColor',[R,G,B]);
set(f,'Position',[f_center(1)-f_width(1)/2,f_center(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5]);

end
