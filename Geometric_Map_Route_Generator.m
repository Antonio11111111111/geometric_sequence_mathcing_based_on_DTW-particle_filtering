function route = Geometric_Map_Route_Generator(map, route_length)
%GEOMETRIC_MAP_ROUTE_GENERATOR Generates a random route on a given map.
%   route = GEOMETRIC_MAP_ROUTE_GENERATOR(map, route_length) generates a
%   random route on the input 'map' using a random walk algorithm.
%
%   Parameters:
%       map: A 1D or 2D matrix representing the geometric map.
%
%       route_length: (Optional) The number of points in the route. 
%                     If not provided, it defaults to 15% of the total 
%                     number of points in the map.
%
%   Output:
%       route: A structure with two fields:
%          .path: An [N x D] matrix of coordinates for the route path,
%                 where N is route_length and D is the map dimension (1 or 2).
%          .intensity: An [N x 1] vector of map values at each point
%                      in the path.

% This is the function of generating a random route, which has 2 elements
% in this structure: location(loc), and its geometric intensity value based on the given map, 
% the details are given below:
% name: geometric_map_route_generator
% param: map: it can be given either in 1 or 2 dimension, the infomation
% will be detected in the function.
% output: route: it is suppoed to be a structure, with the location and the
% geometric intensity value.

% --- Input Validation & Default Parameters ---
if nargin < 1
    error('A map must be provided.');
end

if nargin < 2
    % Set a default route length if not provided
    route_length = floor(numel(map) * 0.15); 
end

% --- Dimension Detection ---
[rows, cols] = size(map);
if rows == 1 || cols == 1
    is_1D = true;
    map_dim = 1;
    map_bounds = numel(map);
else
    is_1D = false;
    map_dim = 2;
    map_bounds = [rows, cols];
end

% Initialize the path matrix
path = zeros(route_length, map_dim);

% --- Generate Route using a Random Walk ---
if is_1D
    % 1D Random Walk
    % Possible moves: left (-1) or right (+1)
    moves = [-1; 1];
    
    % Start at a random position
    path(1) = randi(map_bounds); 
    
    for i = 2:route_length
        current_pos = path(i-1);
        is_valid = false;
        while ~is_valid
            % Pick a random move
            move = moves(randi(numel(moves)));
            next_pos = current_pos + move;
            
            % Check if the move is within bounds
            if next_pos >= 1 && next_pos <= map_bounds
                is_valid = true;
            end
        end
        path(i) = next_pos;
    end
    
    % Get intensity values for the 1D path
    intensity = map(path)';
    
else
    % 2D Random Walk
    % Possible moves: up, down, left, right
    moves = [-1, 0;  % Up
              1, 0;  % Down
              0, -1; % Left
              0, 1]; % Right
           
    % Start at a random [row, col] position
    start_pos = [randi(rows), randi(cols)];
    path(1,:) = start_pos;
    
    for i = 2:route_length
        current_pos = path(i-1, :);
        is_valid = false;
        while ~is_valid
            % Pick a random move
            move = moves(randi(size(moves, 1)), :);
            next_pos = current_pos + move;
            
            % Check if the move is within bounds
            if (next_pos(1) >= 1 && next_pos(1) <= rows && ...
                next_pos(2) >= 1 && next_pos(2) <= cols)
                is_valid = true;
            end
        end
        path(i,:) = next_pos;
    end
    
    % Convert [row, col] subscripts to linear indices to efficiently get values
    linear_indices = sub2ind(size(map), path(:,1), path(:,2));
    intensity = map(linear_indices);
end

% --- Final Output Structure ---
route.path = path;
route.intensity = intensity;

end