% Linear separation figure
color_names = {'Violet', 'Blue', 'Green', 'Yellow', 'Orange', 'Red'}; 

color_ranges = [[380 450]; 
    [450 495]
    [495 570]
    [570 590]
    [590 620]
    [620 750]]; 

color_center=bsxfun(@minus,color_ranges(:, 2), color_ranges(:,1)); 
color_center=color_ranges(:,1) + color_center/2; 


vline(380, 'k-', color_names{1})
for i = 1:length(color_ranges)-1
    vline(color_ranges(i, 2), 'k-', color_names{i+1})
end