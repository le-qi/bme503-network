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

% Code was supposed to be for generation of red patches, but actually if
% you make the background red, you can visualize it just as easily.
% Uncomment to make red patches

% for i = 1:length(x_coords)
%     for j = 1:length(y_coords)
%         center = repmat([x_coords(i) y_coords(j)], 4, 1);
%         delta = [-block_size/2 -block_size/2;...
%             -block_size/2 block_size/2;...
%             block_size/2 block_size/2;...
%             block_size/2 -block_size/2];
%         vert = delta + center; % x and y vertex coordinates
%         fac = [1 2 3 4]; % vertices to connect to make square
%         patch('Faces',fac,'Vertices',vert,'FaceColor','red')
%     end
% end

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

% Test cases for coords_to_color
% c = coords_to_color(32, 17, x_coords, y_coords, points)
% c = coords_to_color(-34, 18, x_coords, y_coords, points)
% c = coords_to_color(19, 24, x_coords, y_coords, points)
% c = coords_to_color(67, -68, x_coords, y_coords, points)
% c = coords_to_color(20, -56, x_coords, y_coords, points)
% c = coords_to_color(-56, 1, x_coords, y_coords, points)
% c = coords_to_color(47, -55, x_coords, y_coords, points)




