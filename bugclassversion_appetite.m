function movebug7
clear
%ANNsetup
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

tgt1 = [150*rand(1,2)-[75. 75.]]; %food target location
tgt2 = [150*rand(1,2)-[75. 75.]]; %food target location
%Make Food
food_lambda1 = (680-380)*rand() + 380; % food color
[R1, G1, B1] = wavelength_to_RGB(food_lambda1); % for the purposes of drawing

food_lambda2 = (680-380)*rand() + 380; % food color
[R2, G2, B2] = wavelength_to_RGB(food_lambda2); % for the purposes of drawing


f_center1 = tgt1;
f_center2 = tgt2;

f_width = 5;
food1 = rectangle('Position',[f_center1(1)-f_width(1)/2,f_center1(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5],'EdgeColor',[R1, G1, B1], 'LineWidth', 4);
set(food1,'Position',[f_center1(1)-f_width(1)/2,f_center1(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5]);

food2 = rectangle('Position',[f_center2(1)-f_width(1)/2,f_center2(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5],'EdgeColor',[R2, G2, B2], 'LineWidth', 4);
set(food2,'Position',[f_center2(1)-f_width(1)/2,f_center2(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5]);


heading_angle=0;



xR = 10;   % Identify the X and Y coordinate of the R and L sensor from the
% robot description
yR = -5;
xL = 10;
yL = 5;

% Compute the magnitude of the sensor info
%
colorAll1 = color_sense(food_lambda1); % single sensor; senses one value
colorAll2 = color_sense(food_lambda2);

sensorL1 = sqrt((xL - f_center1(1)).^2 + (yL - f_center1(2)).^2);
sensorR1 = sqrt((xR - f_center1(1)).^2 + (yR - f_center1(2)).^2);

sensorL2 = sqrt((xL - f_center2(1)).^2 + (yL - f_center2(2)).^2);
sensorR2 = sqrt((xR - f_center2(1)).^2 + (yR - f_center2(2)).^2);

draw_robot(x,y, heading_angle);

last_food1 = 0;
last_food2 = 0;

% Set up stomach system
stomach = 1;
tau_stomach = 100;
a = 1/300; % movement food cost
dstomach_dt = @(stomach, eaten, v) -stomach/tau_stomach + eaten/DT - v*a;

% Set up thirst system; depreciates faster than food but less expended when moving
thirst = 1;
tau_thirst = 50;
b = 1/500;
dthirst_dt = @(thirst, drank, v) -thirst/tau_thirst + drank/DT - v*b;
eaten = 0;
drank = 0;

for k=1:length(tvec)
    
    [motorL, motorR] = brain(sensorL1, sensorR1, sensorL2, sensorR2, colorAll1, colorAll2, stomach, thirst);
    
    vel_left = motorL / 10;   % relate the left and right velocities to motor
    % outputs
    vel_right = motorR / 10;
    
    heading_angle = heading_angle + atan(((vel_right - vel_left)) / 10); % update heading angle and x,y location of robot
    
    x = x + (vel_left + vel_right)./2 * DT * cos(heading_angle); %change x direction
    y = y + (vel_left + vel_right)./2 * DT * sin(heading_angle); %change y direction
    
    % Compute the magnitude of the sensor info
    sensorL1 = sqrt((xL - f_center1(1)).^2 + (yL - f_center1(2)).^2);
    sensorR1 = sqrt((xR - f_center1(1)).^2 + (yR - f_center1(2)).^2);
    
    sensorL2 = sqrt((xL - f_center2(1)).^2 + (yL - f_center2(2)).^2);
    sensorR2 = sqrt((xR - f_center2(1)).^2 + (yR - f_center2(2)).^2);
    
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
    
    DL1 = sqrt((x - f_center1(1)).^2 + (y - f_center1(2)).^2);
    
    if DL1 < 10 || k-last_food1>200 % get new food after certain time, or after eating
        last_food1 = k; % Timer
        tgt1 = [150*rand(1,2)-[75. 75.]]; % New food position
        food_lambda1 = (680-380)*rand() + 380; % New food color
        f_center1 = tgt1;
        colorAll1 = color_sense(food_lambda1); % New neural input
        [R1, G1, B1] = wavelength_to_RGB(food_lambda1);
        
        set(food1,'Position',[f_center1(1)-f_width(1)/2,f_center1(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5], 'EdgeColor', [R1, G1, B1]);
    end
    
    DL2 = sqrt((x - f_center2(1)).^2 + (y - f_center2(2)).^2);
    if DL2 < 10 || k-last_food2>200 % get new food after certain time, or after eating
        last_food2 = k; % Timer
        tgt2 = [150*rand(1,2)-[75. 75.]]; % New food position
        food_lambda2 = (680-380)*rand() + 380; % New food color
        f_center2 = tgt2;
        colorAll2 = color_sense(food_lambda2); % New neural input
        [R2, G2, B2] = wavelength_to_RGB(food_lambda2);
        set(food2,'Position',[f_center2(1)-f_width(1)/2,f_center2(2)-f_width(1)/2,f_width(1),f_width(1)],'Curvature',[0.5,0.5], 'EdgeColor', [R2, G2, B2]);
    end
    
    % Compute food and statute:
    b1 = colorAll1(2)>0.5 || colorAll1(1)>0.5;
    g1 = colorAll1(3)>0.5;
    
    b2 = colorAll2(2)>0.5 || colorAll2(1)>0.5;
    g2 = colorAll2(3)>0.5;
    
    if(DL1 < 10)  % If food 1 was eaten
        if(g1)
            eaten = 1;
        end
        if(b1) % If water
            drank = 1;
        end
    end
    
    if(DL2 < 10)  % If food 2 was eaten
        if(g2)
            eaten = 1;
        end
        if(b2)
            drank = 1;
        end
    end
    
    stomach = stomach + DT*dstomach_dt(stomach, eaten, norm(vel_left, vel_right));
    if(stomach<0.01) stomach = 0.01; end
    thirst = thirst + DT*dthirst_dt(thirst, drank, norm(vel_left, vel_right));
    if(thirst<0.01) thirst =0.01; end
    subplot(1, 2, 2)
    if(mod(k, 4) == 0)
        cla reset
        hold on
        bar(0, stomach, 'g')
        bar(1, thirst, 'b')
        ylim([0 2])
        set(gca, 'XTick', 0:1)
        set(gca, 'XTickLabel', {'Stomach', 'Thirst'})
    end
    subplot(1, 2, 1)
    eaten = 0;
    drank = 0;
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
    subplot(1, 2, 1)
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
    subplot(1, 2, 1)
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

function [motorL, motorR] = brain(sensorL1, sensorR1, sensorL2, sensorR2, colorAll1, colorAll2, stomach, thirst)
color_names = {'Violet', 'Blue', 'Green', 'Yellow', 'Orange', 'Red'}; % Categories, for reference
thresh = 0.5; % threshold for detection

colorAll1 = colorAll1 > 0.5;
b1 = colorAll1(2) || colorAll1(1);
g1 = colorAll1(3);

colorAll2 = colorAll2 > 0.5;
b2 = colorAll2(2) || colorAll2(1);
g2 = colorAll2(3);

% AGGRESSOR % get resources
if(stomach>thirst)
    if(g1 && g2)
        motorL = 250-sensorR1; %min([250 - sensorR1; 250 - sensorR2;]);
        motorR = 250-sensorL1; %min([250 - sensorL1; 250 - sensorL2;]);
        return
    end
    if(g1)
        motorL = 250 - sensorR1;
        motorR = 250 - sensorL1;
        return
    end
    if(g2)
        motorL = 250 - sensorR2;
        motorR = 250 - sensorL2;
        return
    end
    
    if(b1 && b2)
        motorL = 250-sensorR1; %min([250 - sensorR1; 250 - sensorR2;]);
        motorR = 250-sensorL1; %min([250 - sensorL1; 250 - sensorL2;]);
        return
    end
    if(b1)
        motorL = 250 - sensorR1;
        motorR = 250 - sensorL1;
        return
    end
    if(b2)
        motorL = 250 - sensorR2;
        motorR = 250 - sensorL2;
        return
    end
    
end

if(thirst>=stomach)
    if(b1 && b2)
        motorL = 250-sensorR1;%min([250 - sensorR1; 250 - sensorR2;]);
        motorR = 250-sensorL1;%min([250 - sensorL1; 250 - sensorL2;]);
        return
    end
    if(b1)
        motorL = 250 - sensorR1;
        motorR = 250 - sensorL1;
        return
    end
    if(b2)
        motorL = 250 - sensorR2;
        motorR = 250 - sensorL2;
        return
    end
    
    if(g1 && g2)
        motorL = 250-sensorR1;%min([250 - sensorR1; 250 - sensorR2;]);
        motorR = 250-sensorL1;%min([250 - sensorL1; 250 - sensorL2;]);
        return
    end
    if(g1)
        motorL = 250 - sensorR1;
        motorR = 250 - sensorL1;
        return
    end
    if(g2)
        motorL = 250 - sensorR2;
        motorR = 250 - sensorL2;
        return
    end
end


% NEUTRAL % neither water nor food; don't move
motorL = (10+(-5 + 10*rand(1, 1)))/10; % Just jiggle in place
motorR = (10+(-5 + 10*rand(1, 1)))/10;


