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
tstop = 300;
tvec = t:DT:tstop;


food_lambda = [510, 700];
[R1,G1,B1]= wavelength_to_RGB(food_lambda(1));
[R2, G2, B2]= wavelength_to_RGB(food_lambda(2));     


%Make Food
tgt = [150*rand(2)-[75. 75.; 75. 75.]]; %food target location
f_center = tgt;
f_width = 5;

%Green Food
food_1 = rectangle('Position',[f_center(1,1)-f_width(1)/2,f_center(1,2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5],'EdgeColor',[R1, G1, B1],'LineWidth', 4);

%Red Food
food_2 = rectangle('Position',[f_center(2,1)-f_width(1)/2,f_center(2,2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5],'EdgeColor',[R2, G2, B2], 'LineWidth', 4);

heading_angle=0;


% Identify the X and Y coordinate of the R and L sensor from the
% robot description
xR = 10;   
yR = -5;
xL = 10;
yL = 5;

%Sense the color of the food particles with noise
food_lambda_noise = poissrnd(food_lambda);
colorAll = color_sense(food_lambda_noise);  
 
% Compute the magnitude of the sensor info
sensorL(1) = sqrt((xL - f_center(1,1)).^2 + (yL - f_center(1,2)).^2);
sensorR(1) = sqrt((xR - f_center(1,1)).^2 + (yR - f_center(1,2)).^2);

sensorL(2)= sqrt((xL - f_center(2,1)).^2 + (yL - f_center(2,2)).^2);
sensorR(2)= sqrt((xR - f_center(2,1)).^2 + (yR - f_center(2,2)).^2);

draw_robot(x,y, heading_angle);

for k=1:length(tvec) 
    
    %Brain function
    [motorL, motorR] = brain(sensorL, sensorR, colorAll);

    %Compute the ;eft and right velocities
    vel_left = motorL / 50;  
    vel_right = motorR / 50;
    
    heading_angle = heading_angle + atan(((vel_right - vel_left)) / 10); % update heading angle and x,y location of robot
    
    x = x + (vel_left + vel_right)./2 * DT * cos(heading_angle); %change x direction
    y = y + (vel_left + vel_right)./2 * DT * sin(heading_angle); %change y direction

  % Identify the new X and Y coordinate of the R and L sensor from the
    % robot description
    xR = new_verts(1);   
    yR = new_verts(14);
    xL = new_verts(11);
    yL = new_verts(24);
     
% Compute the magnitude of the sensor info for both foods 
    sensorL(1) = sqrt((xL - f_center(1,1)).^2 + (yL - f_center(1,2)).^2);
    sensorR(1) = sqrt((xR - f_center(1,1)).^2 + (yR - f_center(1,2)).^2);

    sensorL(2)= sqrt((xL - f_center(2,1)).^2 + (yL - f_center(2,2)).^2);
    sensorR(2)= sqrt((xR - f_center(2,1)).^2 + (yR - f_center(2,2)).^2);

    %Arena wraps on itself
    if x>arena_size
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
    
  
    DL = sqrt((x - f_center(1,1)).^2 + (y - f_center(1,2)).^2);
    
    %Creating new green food on it getting eaten
    if DL < 10  
        
        tgt = [150*rand(1,2)-[75. 75.]]; % New food position
        f_center(1,:) = tgt;
        set(food_1,'Position',[f_center(1,1)-f_width(1)/2,f_center(1,2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5]);
        
        food_lambda_noise(1) = poissrnd(food_lambda(1)); 
        colorAll = color_sense(food_lambda_noise); 
    end
    
    % Creating red food at a new location after a certain time duration
    if mod(k,600) == 0 
        
        tgt = [150*rand(1,2)-[75. 75.]]; % New food position
        f_center(2,:) = tgt;
        set(food_2,'Position',[f_center(2,1)-f_width(1)/2,f_center(2,2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5]);
        
        food_lambda_noise(2) = poissrnd(food_lambda(2));
        colorAll = color_sense(food_lambda_noise); 
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
   
  %Average out sensed information from left and right sensors
  sensor_avg(1) = (sensorL(1) + sensorR(1))/2;
  sensor_avg(2) = (sensorL(2) + sensorR(2))/2;
  
  %Check for closest food particle and accordingly act
  if sensor_avg(1) > sensor_avg(2)
       if sensorL(2) < sensorR(2)
           motorL = 300 - sensorL(2); 
           motorR = 300 - sensorR(2);
       else
           motorL = 300 - sensorL(2); 
           motorR = 300 - sensorR(2);
       end
       title('Get away from Red!!', 'FontSize', 16);
  else
       motorL = 300 - sensorR(1);
       motorR = 300 - sensorL(1);
       title('Go for Green!!', 'FontSize', 16);
  end
  
  
  
  

  