% global arena_size;
% global block_size;
% global arena;

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arena_size = 75; %bigger arena makes things more interesting.
axis([-arena_size arena_size -arena_size arena_size]);  % create dummy arena
set(gca,'color','r');

block_size = 25;

x_coords = -arena_size:block_size:arena_size;
y_coords = -arena_size:block_size:arena_size;

arena = figure(1);

for i = 1:length(x_coords)
    for j = 1:length(y_coords)
        center = repmat([x_coords(i) y_coords(j)], 4, 1);
        delta = [-block_size/2 -block_size/2;...
            -block_size/2 block_size/2;...
            block_size/2 block_size/2;...
            block_size/2 -block_size/2];
        vert = delta + center; % x and y vertex coordinates
        fac = [1 2 3 4]; % vertices to connect to make square
        patch('Faces',fac,'Vertices',vert,'FaceColor','red')
    end
end

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


