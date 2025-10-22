function next_location = Get_Next_Step(map, last_location)
%
% function description:
%   Get_Next_Step simulates a single step of a random walk on a map. 
%   It takes the last known location and returns a new, valid location
%   that is one step away (up, down, left, right, or stay).
%
% parameter:
%   map:           [NxM double] The 1D or 2D map.
%   last_location: [1xD double] The starting coordinate(s). 
%                  For 1D, this is a single index.
%                  For 2D, this is a [row, col] vector.
%
% output:
%   next_location: [1xD double] The new, valid coordinate(s) after one step.
%
% example of the function(usage):
%
%   % --- 1D Example ---
%   map_1d = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
%   current_pos_1d = 5;
%   next_pos_1d = Get_Next_Step(map_1d, current_pos_1d);
%   % >> next_pos_1d will be 4, 5, or 6
%
%   % --- 2D Example ---
%   map_2d = magic(10);
%   current_pos_2d = [5, 5]; % Start at [row 5, col 5]
%   next_pos_2d = Get_Next_Step(map_2d, current_pos_2d);
%   % >> next_pos_2d will be a valid neighbor like [4,5] or [5,6] etc.
%

    % --- Dimension Detection ---
    [rows, cols] = size(map);
    if rows == 1 || cols == 1
        is_1D = true;
        map_bounds = numel(map);
    else
        is_1D = false;
    end

    % --- Generate Next Step using a Random Walk ---
    is_valid = false;
    while ~is_valid
        if is_1D
            % 1D moves: left (-1), stay (0), or right (+1)
            move = randi([-1, 1]);
            next_location = last_location + move;
            
            % Check if the move is within bounds
            if next_location >= 1 && next_location <= map_bounds
                is_valid = true;
            end
        else % 2D
            % 2D moves: up, down, left, right, stay
            moves = [-1, 0;  % Up
                      1, 0;  % Down
                      0, -1; % Left
                      0, 1;  % Right
                      0, 0]; % Stay
            
            % Pick a random move
            move = moves(randi(size(moves, 1)), :);
            next_location = last_location + move;
            
            % Check if the move is within bounds
            if (next_location(1) >= 1 && next_location(1) <= rows && ...
                next_location(2) >= 1 && next_location(2) <= cols)
                is_valid = true;
            end
        end
    end
end