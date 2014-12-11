function movebug7
clear
clf
close all

global new_verts fverts
x=0;
y=0;

% Code to generate the obstacle course.

% Field is divided into cells of block_size length and width. Each field is
% assigned either red or green as determined by points, which hold the
% centers of the green (clear) cells. Red cells are bad, green cells are
% good.

global block_size;              % size of patch separation
global x_coords y_coords;       % coordinates of all cells
global points;                  % coordinates of green cells

clf
arena_size = 75;
axis([-arena_size arena_size -arena_size arena_size]);
set(gca,'color','r');           % background red

block_size = 25;

x_coords = -arena_size:block_size:arena_size; % x-coords for all cells
y_coords = -arena_size:block_size:arena_size; % y-coords for all cells

arena = figure(1);

% Create points for the centers of cell along path - change as necessary to
% make different paths
points = [-50 75;
    -50 50;
    -50 25;
    -50 0;
    -25 0;
    0 0;
    25 0;
    50 0;
    50 -25;
    50 -50;
    50 -75];

% Code actually generates patches and displays them on plot
for i = 1:length(points)
    center = repmat(points(i, :), 4, 1);
    delta = [-block_size/2 -block_size/2;...
        -block_size/2 block_size/2;...
        block_size/2 block_size/2;...
        block_size/2 -block_size/2];
    vert = delta + center; % x and y vertex coordinates
    fac = [1 2 3 4]; % vertices to connect to make square
    patch('Faces',fac,'Vertices',vert,'FaceColor','green')
end

draw_robot() % show robot for the first time (no arguments)

t = 0;
DT=.1;
tstop = 200;
tvec = t:DT:tstop;


heading_angle=0;


% Identify the X and Y coordinate of the R and L sensor from the
% robot description
xR = 10;   
yR = -5;
xL = 10;
yL = 5;

draw_robot(x,y, heading_angle);



% Sense color with added noise
lambda(1) = coords_to_color(xL, yL, x_coords, y_coords, points);
lambda(2) = coords_to_color(xR, yR, x_coords, y_coords, points);

lambda_noise = poissrnd(lambda);
colorAll = color_sense(lambda_noise); 

% Compute the magnitude of the sensor info
sensorL = colorAll(:,1);
sensorR = colorAll(:,1);


for k=1:length(tvec) 
    
    %Brain function implementation
    [motorL, motorR] = brain(sensorL, sensorR);

    %Compute left and right velocities
    vel_left = motorL / 60;   
    vel_right = motorR / 60;
    
    heading_angle = heading_angle + atan(((vel_right - vel_left)) / 10); % update heading angle and x,y location of robot
    
    x = x + (vel_left + vel_right)./2 * DT * cos(heading_angle); %change x direction
    y = y + (vel_left + vel_right)./2 * DT * sin(heading_angle); %change y direction

    % Identify the new X and Y coordinate of the R and L sensor from the
    % robot description
    xR = new_verts(1);   
    yR = new_verts(14);
    xL = new_verts(11);
    yL = new_verts(24);
    
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
    
    % Sense color with added noise
    lambda(1) = coords_to_color(xL, yL, x_coords, y_coords, points);
    lambda(2) = coords_to_color(xR, yR, x_coords, y_coords, points);

    lambda_noise = poissrnd(lambda);
    colorAll = color_sense(lambda_noise); 

    % Compute the magnitude of the sensor info
    sensorL = colorAll(:,1);
    sensorR = colorAll(:,2);
    
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

%Function to sense color from the surrounding    
 function output = color_sense(lambda)
    global tune_l
    n1 = f(lambda, tune_l(1)); 
    n2 = f(lambda, tune_l(2)); 
    n3 = f(lambda, tune_l(3)); 
    n4 = f(lambda, tune_l(4)); 
    output = color_net([n1; n2; n3; n4;]); 

% Converts a set of coordinates to the color found in the cell
% Right now it's red or green, can change to wavelengths as necessary
    function [color] = coords_to_color(x, y, xpatches, ypatches, coords)


    xerror = abs(xpatches - x);
    x_center = xpatches(find(min(xerror) == xerror));
    
    yerror = abs(ypatches - y);
    y_center = ypatches(find(min(yerror) == yerror));
    
    for i = 1:length(coords)
        if (coords(i,1) == x_center && coords(i,2)== y_center)
            color = 510;
            return;
        end
    end
    
    color = 700;
    
function [motorL, motorR] = brain(sensorL, sensorR)
  color_names = {'Violet', 'Blue', 'Green', 'Yellow', 'Orange', 'Red'}; % Categories, for reference
  thresh = 0.5; % threshold for detection
  
  % Compute a relationship of sensor input to motor output
  
  %Both sensors sense red
  if(sensorL(6)>thresh && sensorR(6)>thresh)
    motorL = 50 - sensorL(6);
    motorR = 50 - 20*sensorR(6);
  end
  
  %Left sensor senses red and right senses green
  if(sensorL(6)>thresh && sensorR(6)<thresh)
    motorL = 250;
    motorR = 250 - 150*sensorL(6);
  end
  
  %Left sensor senses green and right senses red
  if (sensorL(6)<thresh && sensorR(6)>thresh)
    motorL = 250 - 150*sensorR(6);
    motorR = 250;
  end
  
  %Both sensors sense green
  if (sensorL(6)<thresh && sensorR(6)<thresh)
      motorL = 250 - 150*sensorR(6);
      motorR = 250 - 150*sensorL(6);
  end
  

  
  
  
  
 