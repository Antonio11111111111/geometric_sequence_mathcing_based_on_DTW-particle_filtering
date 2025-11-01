function [min_distance, target_location] = Geomectric_Map_Matching_DTW_1D(geo_map, live_reading, show_figure)
%
% function description:
%   Geomectric_Map_Matching_DTW_1D finds the best-matching segment 
%   for a short sequence (`live_reading`) within a longer sequence 
%   (`geo_map`) using Subsequence Dynamic Time Warping (DTW).
%
%   This algorithm allows the match to start at any point in the 
%   `geo_map` and finds the path with the minimum total cost that 
%   aligns the *entire* `live_reading`.
%
% parameter:
%   geo_map:      [1xN double] The long 1D reference sequence (database/map).
%   live_reading: [1xM double] The short 1D query sequence to be matched.
%   show_figure:  (Optional) [logical] If true, generates a plot of the
%                 match. Defaults to false.
%
% output:
%   min_distance:    [double] The minimum accumulated DTW cost of the best match.
%   target_location: [1xK double] The actual vector segment from `geo_map` 
%                    that corresponds to the best match.
%
% example of the function(usage):
%
%   % 1. Define the map and the reading
%   geo_map = [1.2, 1.1, 1.3, 2.0, 5.0, 8.0, 9.0, 8.5, 4.0, 2.1, 1.4, 1.3, 1.5];
%   live_reading = [5.2, 7.8, 9.1, 9.0, 8.2, 3.5];
%
%   % --- Usage 1: Get outputs only (no figure) ---
%   [dist, loc] = Geomectric_Map_Matching_DTW_1D(geo_map, live_reading);
%   fprintf('Best match distance: %f\n', dist);
%   disp(loc);
%
%   % --- Usage 2: Get outputs AND show the figure ---
%   [dist_fig, loc_fig] = Geomectric_Map_Matching_DTW_1D(geo_map, live_reading, true);
%

    % --- 0. Handle Optional Input ---
    % 'nargin' is the number of inputs provided to the function
    if nargin < 3
        show_figure = false; % Default: do not show the figure
    end

    m = length(live_reading); 
    n = length(geo_map);      
    fprintf('Matching short sequence (len %d) against long map (len %d)...\n', m, n);

    % --- 2. Build the Cost Matrices ---
    local_cost = (live_reading' - geo_map).^2;
    D = zeros(m, n);
    
    % Initialize first row
    D(1, :) = local_cost(1, :);
    
    % Initialize first column
    for i = 2:m
        D(i, 1) = local_cost(i, 1) + D(i-1, 1);
    end
    
    % Fill the rest of the matrix
    for i = 2:m
        for j = 2:n
            D(i, j) = local_cost(i, j) + min([D(i-1, j), D(i, j-1), D(i-1, j-1)]);
        end
    end
    
    % --- 3. Find Best Match ---
    [min_distance, end_index] = min(D(m, :));
    fprintf('Best match found with total distance: %f\n', min_distance);
    
    
    % --- 4. Backtrack to Find the Start Index ---
    path = [m, end_index]; 
    i = m;
    j = end_index;
    
    while i > 1
        if j == 1
            i = i - 1;
        else
            [~, step] = min([D(i-1, j), D(i, j-1), D(i-1, j-1)]);
            switch step
                case 1 % Up
                    i = i - 1;
                case 2 % Left
                    j = j - 1;
                case 3 % Diagonal
                    i = i - 1;
                    j = j - 1;
            end
        end
        path = [i, j; path]; 
    end
    
    start_index = path(1, 2); 
    fprintf('Match *starts* at index %d in the geo_map.\n', start_index);
    
    % --- 5. Prepare Output ---
    target_location = start_index:end_index;
    
    % --- 6. (Optional) Visualize the Result ---
    if show_figure
        figure;
        plot(1:n, geo_map, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Geomagnetic Map (Long)');
        hold on;
        
        plot(start_index:end_index, geo_map(start_index:end_index), 'r-s', 'LineWidth', 2.5, 'DisplayName', 'Best Match');
        
        plot((start_index:start_index+m-1), live_reading, 'g--*', 'LineWidth', 1.5, 'DisplayName', 'Live Reading (Short)');
        
        title('Subsequence DTW for Geomagnetic Matching (Corrected)');
        xlabel('Map Index');
        ylabel('Geomagnetic Value');
        legend;
        grid on;
    end
end