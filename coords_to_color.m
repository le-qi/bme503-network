function [ color ] = coords_to_color(x, y, xpatches, ypatches, coords)
% Converts a set of coordinates to the color found in the cell
% Right now it's red or green, can change to wavelengths as necessary

    xerror = abs(xpatches - x)
    x_center = xpatches(find(min(xerror) == xerror))
    
    yerror = abs(ypatches - y)
    y_center = ypatches(find(min(yerror) == yerror));
    
    for i = 1:length(coords)
        if (coords(i, 1) == x_center && coords(i, 2) == y_center)
            color = 'green';
            return;
        end
    end
    
    color = 'red';

