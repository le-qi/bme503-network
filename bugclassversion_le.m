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

tgt = [150*rand(1,2)-[75. 75.]]; %food target location

%Make Food
f_center = tgt;
f_width = 5;
food = rectangle('Position',[f_center(1)-f_width(1)/2,f_center(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5],'EdgeColor','w');
set(food,'Position',[f_center(1)-f_width(1)/2,f_center(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5]);

heading_angle=0;



xR = 10;   % Identify the X and Y coordinate of the R and L sensor from the
          % robot description
yR = -5;
xL = 10;
yL = 5;

% Compute the magnitude of the sensor info
%    
l = 700; % nm; test
colorAll = color_sense(l); % single sensor; senses one value 
colorAll 

sensorL = sqrt((xL - f_center(1)).^2 + (yL - f_center(2)).^2);
sensorR = sqrt((xR - f_center(1)).^2 + (yR - f_center(2)).^2);
draw_robot(x,y, heading_angle);

for k=1:length(tvec) 
    
    [motorL, motorR] = brain(sensorL, sensorR, colorAll);

    vel_left = motorL / 10;   % relate the left and right velocities to motor
                                % outputs
    vel_right = motorR / 10;
    
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
 
    % heading_angle = pi/2;
    draw_robot(x,y, heading_angle);
    drawnow
    
    xR = new_verts(1);   % Identify the X and Y coordinate of the R and L sensor from the
    %         robot description
    yR = new_verts(14);
    xL = new_verts(11);
    yL = new_verts(24);
    
    DL = sqrt((x - f_center(1)).^2 + (y - f_center(2)).^2);
    
    if DL < 10
        tgt = [150*rand(1,2)-[75. 75.]];
        f_center = tgt;
        set(food,'Position',[f_center(1)-f_width(1)/2,f_center(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5]);
    end
end



function draw_robot(x,y,heading_angle)
persistent nverts npts verts linewidth color hrobot
global new_verts fverts 

%persistent variables are local to the function but kept in memory for the
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
  
  % COWARD (uncomment for action); % coward if it detects red
  if(colorAll(6)>thresh)
   motorL = 250 - sensorL; 
   motorR = 250 - sensorR; 
   return
  end
  
  % AGGRESSOR
  motorL = 250 - sensorR;   % compute a relationship of sensor input to motor output
  motorR = 250 - sensorL;
  

 