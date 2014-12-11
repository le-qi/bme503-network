function movebug7
clear
clf
close all

global new_verts fverts
x=0;
y=0;

arena_size=75; %bigger arena makes things more interesting.
draw_robot() % show robot for the first time (no arguments)
axis([-arena_size arena_size -arena_size arena_size]);  % create dummy arena
set(gca,'color','k');
axis square
set(gca,'Visible','on');

t = 0;
DT=.1;
tstop = 200;
tvec = t:DT:tstop;



%Food color
food_lambda = (750-380)*rand() + 380; % Food color
% Conversion of color to RGB values
[R, G, B] = wavelength_to_RGB(food_lambda);      

%Make Food
tgt = [150*rand(1,2)-[75. 75.]]; %Food target location
f_center = tgt;
f_width = 5;
food = rectangle('Position',[f_center(1)-f_width(1)/2,f_center(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5],'EdgeColor',[R, G, B], 'LineWidth', 4);
set(food,'Position',[f_center(1)-f_width(1)/2,f_center(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5]);

heading_angle=0;


% Identify the X and Y coordinate of the R and L sensor from the
% robot description
xR = 10;   
yR = -5;
xL = 10;
yL = 5;


% Sense color with added noise
food_lambda_noise = poissrnd(food_lambda);
colorAll = color_sense(food_lambda); 

% Compute the magnitude of the sensor info
sensorL = sqrt((xL - f_center(1)).^2 + (yL - f_center(2)).^2);
sensorR = sqrt((xR - f_center(1)).^2 + (yR - f_center(2)).^2);
draw_robot(x,y, heading_angle);

last_food = 0; % counter for changing food location and color
for k=1:length(tvec) 
    
    %Brain function implementation
    [motorL, motorR] = brain(sensorL, sensorR, colorAll);

    %Compute left and right velocities
    vel_left = motorL / 40;   
    vel_right = motorR / 40;
    
    heading_angle = heading_angle + atan(((vel_right - vel_left)) / 10); % update heading angle and x,y location of robot
    
    x = x + (vel_left + vel_right)./2 * DT * cos(heading_angle); %change x direction
    y = y + (vel_left + vel_right)./2 * DT * sin(heading_angle); %change y direction

% Compute the magnitude of the sensor info 
    sensorL = sqrt((xL - f_center(1)).^2 + (yL - f_center(2)).^2);
    sensorR = sqrt((xR - f_center(1)).^2 + (yR - f_center(2)).^2);

    if x>arena_size % have arena wrap on itself
        x=-arena_size;
    end
    if y>arena_size
        y=-arena_size;
    end
    if x<-arena_size
        x=arena_size;
    end
    if y<-arena_size
        y=arena_size;
    end
 
    draw_robot(x,y, heading_angle);
    drawnow
    
    
    % Identify the new X and Y coordinate of the R and L sensor from the
    % robot description
    xR = new_verts(1);   
    yR = new_verts(14);
    xL = new_verts(11);
    yL = new_verts(24);
    
    DL = sqrt((x - f_center(1)).^2 + (y - f_center(2)).^2);
    
    %Get new food particle if eaten or
    %too much time has elapsed
    if DL < 10 || k-last_food>200 
        
        last_food = k; % Timer is updated 
        tgt = [150*rand(1,2)-[75. 75.]]; % New food position
        food_lambda = (750-380)*rand() + 380; % New food color 
        [R, G, B] = wavelength_to_RGB(food_lambda); 

        f_center = tgt;
        set(food,'Position',[f_center(1)-f_width(1)/2,f_center(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5], 'EdgeColor', [R, G, B]);
        
        food_lambda_noise = poissrnd(food_lambda);
        colorAll = color_sense(food_lambda_noise); % New neural input 
                
    end
end



function draw_robot(x,y,heading_angle)
persistent nverts npts verts linewidth color hrobot
global new_verts fverts 

%Persistent variables are local to the function but kept in memory for the
%next function call


if(nargin < 3)
    % initialize
    [verts linewidth color] = robot_coords;
    fverts=verts;
    nverts = size(verts,1);
    npts = size(verts, 2);
    figure(1);
    hold on;
    % Show the robot's current location
    for i=1:nverts
        hrobot(i) = plot(verts(i,1:(npts/2)), verts(i,(npts/2+1):npts),'color',color(i,:),'linewidth',linewidth(i));
    end
    
else
    % redraw the robot in the new configuration
    tmp_verts = verts;
    
    % adjust heading - rotate the entire robot
    cosa = cos(heading_angle);
    sina = sin(heading_angle);
    new_verts(:,1:(npts/2)) =  x + cosa * tmp_verts(:,1:(npts/2)) - sina * tmp_verts(:,(npts/2+1):npts);
    new_verts(:,(npts/2+1):npts) =  y + sina * tmp_verts(:,1:(npts/2)) + cosa * tmp_verts(:,(npts/2+1):npts);
    
    figure(1);
    hold on
    for i=1:nverts
        set(hrobot(i),'xdata',new_verts(i,1:(npts/2)),'ydata',new_verts(i,(npts/2+1):npts));
    end
    

   
end
%--------------------------------------------------------------------------
% end of draw_robot
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% robot_coords
%--------------------------------------------------------------------------
function [verts, width, color] = robot_coords

height=10;
width=5;

% the first set of coordinates are for x, the next set are for y.  All
% hard-coded in current iteration..

c=1;

  
xc=[10. 7.5   2.5  0. -7.5 -10 -7.5 0. 2.5 7.5 10. 7.5 7.5];
yc=[-5.  -2.5 -2.5 -5. -5. 0 5. 5. 2.5 2.5  5.   2.5 -2.5];

verts(c, :) = [xc'; yc'];



width(c) = 2; color(c,:) = [1 1 0]; c=c+1;
%--------------------------------------------------------------------------
% end of robot_coords
%--------------------------------------------------------------------------      
function output = color_sense(lambda)
    global tune_l
    n1 = f(lambda, tune_l(1)); 
    n2 = f(lambda, tune_l(2)); 
    n3 = f(lambda, tune_l(3)); 
    n4 = f(lambda, tune_l(4)); 
    output = color_net([n1; n2; n3; n4;]); 
    
function [motorL, motorR] = brain(sensorL, sensorR, colorAll)
  color_names = {'Violet', 'Blue', 'Green', 'Yellow', 'Orange', 'Red'}; % Categories, for reference
  thresh = 0.5; % threshold for detection
  
  % Compute a relationship of sensor input to motor output
  
  %Behaviour on detection danger aka red - Coward
  if(colorAll(6)>thresh)
   motorL = 250 - sensorL; 
   motorR = 250 - sensorR; 
   title('Ahh!! DANGER', 'FontSize', 20)
   return
  end
  
  % Behaviour on detecting food - Aggressor
  title('Yay!! Food.', 'FontSize', 16)
  motorL = 300 - sensorR;   
  motorR = 300 - sensorL;
