function movebug7
<<<<<<< HEAD
=======
% Michael Lee, based on 
% work of Dr. Craig Henriquez
%% Initialize
>>>>>>> 965b0a9c1bc04ad8c056660b9a30673a8db97aea
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

<<<<<<< HEAD
t = 0;
=======
%% Time vector
t = 0
>>>>>>> 965b0a9c1bc04ad8c056660b9a30673a8db97aea
DT=.1;
tstop = 200;
tvec = t:DT:tstop;

<<<<<<< HEAD
tgt = [150*rand(1,2)-[75. 75.]]; %food target location

%Make Food
=======


% Make Food
tgt = [150*rand(1,2)-[75. 75.]]; %food target location
>>>>>>> 965b0a9c1bc04ad8c056660b9a30673a8db97aea
f_center = tgt;
f_width = 5;
food = rectangle('Position',[f_center(1)-f_width(1)/2,f_center(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5],'EdgeColor','w');
set(food,'Position',[f_center(1)-f_width(1)/2,f_center(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5]);


<<<<<<< HEAD
heading_angle=0;



xR = 10;   % Identify the X and Y coordinate of the R and L sensor from the
          % robot description
yR = -5;
xL = 10;
yL = 5;

% Compute the magnitude of the sensor info
%                                                                                                                                                                                                                                   
sensorL = sqrt((xL - f_center(1)).^2 + (yL - f_center(2)).^2);
sensorR = sqrt((xR - f_center(1)).^2 + (yR - f_center(2)).^2);
draw_robot(x,y, heading_angle);

for k=1:length(tvec) 
    
    [motorL, motorR] = brain(sensorL, sensorR);

    vel_left = motorL / 10;   % relate the left and right velocities to motor
                                % outputs
    vel_right = motorR / 10;
    
    heading_angle = heading_angle + atan(((vel_right - vel_left)) / 10); % update heading angle and x,y location of robot
    
    x = x + (vel_left + vel_right)./2 * DT * cos(heading_angle); %change x direction
    y = y + (vel_left + vel_right)./2 * DT * sin(heading_angle); %change y direction

% Compute the magnitude of the sensor info 
    sensorL = sqrt((xL - f_center(1)).^2 + (yL - f_center(2)).^2);
    sensorR = sqrt((xR - f_center(1)).^2 + (yR - f_center(2)).^2);

=======
% Initial bug heading and antenna position
heading_angle=0;
xR = 15;
yR = -3;
xL = 15;
yL = 3;

% Compute the magnitude of the sensor info
sensorL= sense(xL, yL, tgt(1), tgt(2));
sensorR= sense(xR, yR, tgt(1), tgt(2));
draw_robot(x,y, heading_angle);


% Construct a questdlg specifying bug type
bugtype = questdlg('Choose your bug type.', ...
    'Braitenbug',...
	'Aggressor', ...
	'Coward', 'Aggressor');

% Running simulation:
for k=1:length(tvec)
    
    
    % motor output = f(sensor readings)
    [motorL, motorR] = brain(sensorL, sensorR, bugtype);
    vel_left= motorL
    vel_right= motorR
    
    % update heading angle and x,y location of robot
    heading_angle = heading(heading_angle, motorL, motorR);
    x = x +(motorL+motorR)*cos(heading_angle); %change x direction
    y = y + (motorL+motorR)*sin(heading_angle); %change y direction
    
    % plot path of bug, for fun.
    % plot(x, y)
    
    %Identify the X and Y coordinate of the R and L sensor from the robot description
    
    alpha = atan(-3/15); % offset of sensor from heading angle
    r = dist(0, 0, 3, 15);
    
    dx=r*cos(heading_angle+alpha);
    dy=r*sin(heading_angle+alpha);
    xR = x+dx;
    yR = y+dy;
    
    dx=r*cos(heading_angle-alpha);
    dy=r*sin(heading_angle-alpha);
    xL = x+dx;
    yL = y+dy;
    
    %Plots sensor locations
    %plot(xR, yR, 'ro')
    %plot(xL, yL, 'ro')
    
    % Compute the magnitude of the sensor info
    sensorL= sense(xL, yL, tgt(1), tgt(2));
    sensorR= sense(xR, yR, tgt(1), tgt(2));
    
    
>>>>>>> 965b0a9c1bc04ad8c056660b9a30673a8db97aea
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
<<<<<<< HEAD
 
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
=======
    
    draw_robot(x,y, heading_angle);
    drawnow
    
    DL = dist(x, y, tgt(1), tgt(2)); % distance to food
    
    if DL <10
        tgt = [150*rand(1,2)-[75. 75.]];
        f_center = tgt;
        set(food,'Position',[f_center(1)-f_width(1)/2,f_center(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5]);
        
        
>>>>>>> 965b0a9c1bc04ad8c056660b9a30673a8db97aea
    end
end



function draw_robot(x,y,heading_angle)
persistent nverts npts verts linewidth color hrobot
<<<<<<< HEAD
global new_verts fverts 
=======
global new_verts fverts
>>>>>>> 965b0a9c1bc04ad8c056660b9a30673a8db97aea

%persistent variables are local to the function but kept in memory for the
%next function call


if(nargin < 3)
    % initialize
    [verts linewidth color] = robot_coords;
    fverts=verts;
<<<<<<< HEAD
    nverts = size(verts,1);
=======
    nverts = size(verts,1)
>>>>>>> 965b0a9c1bc04ad8c056660b9a30673a8db97aea
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
    
<<<<<<< HEAD

   
=======
    
    
>>>>>>> 965b0a9c1bc04ad8c056660b9a30673a8db97aea
end
%--------------------------------------------------------------------------
% end of draw_robot
%--------------------------------------------------------------------------

<<<<<<< HEAD
=======


>>>>>>> 965b0a9c1bc04ad8c056660b9a30673a8db97aea
%--------------------------------------------------------------------------
% robot_coords
%--------------------------------------------------------------------------
function [verts, width, color] = robot_coords

<<<<<<< HEAD
=======

>>>>>>> 965b0a9c1bc04ad8c056660b9a30673a8db97aea
height=10;
width=5;

% the first set of coordinates are for x, the next set are for y.  All
% hard-coded in current iteration..

c=1;

<<<<<<< HEAD
  
xc=[10. 7.5   2.5  0. -7.5 -10 -7.5 0. 2.5 7.5 10. 7.5 7.5];
yc=[-5.  -2.5 -2.5 -5. -5. 0 5. 5. 2.5 2.5  5.   2.5 -2.5];

verts(c, :) = [xc'; yc'];
=======
%xc=[10. 7.5   0.  0. -10. -10. 0. 0. 7.5 10. 7.5 7.5];
%yc=[-5.  -2.5 -2.5 -5. -5.  5. 5. 2.5 2.5  5.   2.5 -2.5];

xy=[0 1;
    
2 1;
-1 3;
2 1;

5 1;
3 2.5;
5 1;

8 1;
6 2.5;
8 1;

10 1;

10 0.8;
15 3;
10 0.7;

10 -0.8;
15 -3;
10 -0.7;

10 -1;

8 -1;
6 -2.5;
8 -1;

5 -1;
3 -2.5;
5 -1;

2 -1;
-1 -3;
2 -1;

0 -1;
0 1;];

xc=xy(:, 1)';
yc=xy(:, 2)';

verts(c, :) = [xc'; yc'];1
>>>>>>> 965b0a9c1bc04ad8c056660b9a30673a8db97aea



width(c) = 2; color(c,:) = [1 1 0]; c=c+1;
%--------------------------------------------------------------------------
% end of robot_coords
<<<<<<< HEAD
%--------------------------------------------------------------------------      

function [motorL, motorR] = brain(sensorL, sensorR)
  
  % COWARD (uncomment for action)
  motorL = 250 - sensorL
  motorR = 250 - sensorR
  
  % AGGRESSOR
  % motorL = 250 - sensorR;   % compute a relationship of sensor input to motor output
  % motorR = 250 - sensorL;
 
=======
%--------------------------------------------------------------------------

function [motorL, motorR] = brain(sensorL, sensorR, bugtype)

switch bugtype
    case 'Aggressor' % Minimize difference between sensors. 
        gain=0.5;
        
        % Logic: some default speed for both, each adjusted by difference in
        % sensors.
        diff=sensorL-sensorR; 
        motorL=gain*(1-diff);
        motorR=gain*(1+diff);
        
    case 'Coward' % If close to food, run away fast. If far, slow down. 
       gain=0.1;
        motorL=gain*sensorL;
        motorR=gain*sensorR;
end


% distance function, for convenience. Returns distance between two sets of (x, y) points
function distance = dist(x, y, xt, yt)
distance = sqrt((x-xt)^2+(y-yt)^2);

% transducing function for the antennae.
function mag = sense(xs, ys, xt, yt)
% xs: x location of sensor
% ys: y location of sensor
% xt: x location of food
% yt: y location of food

distance = sqrt((xs-xt)^2+(ys-yt)^2);
mag = sqrt(distance);

% Adjusting the heading angle based on motors.
function theta= heading(theta0, motorL, motorR)
gain=1;

% If gain is too low, the bug can only make incremental changes to path.
% If its too high, it will spazz out and spin at the slightest stimuli.
% A good balance is found at gain=1.

theta = theta0+(motorL-motorR)*gain;



>>>>>>> 965b0a9c1bc04ad8c056660b9a30673a8db97aea
