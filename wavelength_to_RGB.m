function [R, G, B] = wavelength_to_RGB(lambda)
% This converts a scalar, lambda (nm) into RGB values. 
% Lambda can be a nx1 vector. 

% Adapted by Michael Lee, Duke University, 12/9/2014
% From Dan Bruton's Work

R = zeros(size(lambda)); 
G = zeros(size(lambda)); 
B = zeros(size(lambda)); 

for i = 1:length(lambda)
    l = lambda(i);
    
    if 380<=l && l<=440
        R(i) = -(l-440)/(440-380);
        G(i) = 0;
        B(i) = 1;
    end
    if 440<=l && l<=490
        R(i) = 0;
        G(i) = (l-440)/(490-440);
        B(i) = 1;
    end
    if 490<=l && l<=510
        R(i) = 0;
        G(i) = 1;
        B(i) = -(l-510)/(510-490);
    end
    if 510<=l && l<=580
        R(i) = (l-510)/(580-510);
        G(i) = 1;
        B(i) = 0;
    end
    if 580<=l && l<=645
        R(i) = 1;
        G(i) = -(l-645)/(645-580);
        B(i) = 0;
    end
    if 645<=l && l<=780
        R(i) = 1;
        G(i) = 0;
        B(i) = 0;
    end
end
