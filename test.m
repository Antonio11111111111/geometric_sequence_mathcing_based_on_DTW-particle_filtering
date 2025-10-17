%%
% All the test fo the function will be put HERE
%%
% 1. Define parameters for a 2D map
dimension = 2;
map_size = [50, 75]; % 50 rows, 75 columns

% 2. Generate the geometric map
fprintf('Generating a %dD map of size %dx%d...\n', dimension, map_size(1), map_size(2));
geometric_map = Geometric_Map_Generator(dimension, map_size);

% 3. Display some information about the generated map
fprintf('Map generated successfully. Size: %s\n', mat2str(size(geometric_map)));
fprintf('Minimum value: %.4f\n', min(geometric_map(:)));
fprintf('Maximum value: %.4f\n', max(geometric_map(:)));

% 4. Visualize the map
figure; % Create a new figure window
imagesc(geometric_map); % Display the matrix as an image
colorbar; % Add a color bar to indicate values
colormap(parula); % Use the 'parula' colormap (good for continuous data)
title(sprintf('Generated 2D Geometric Map (%dx%d)', map_size(1), map_size(2)));
xlabel('Column Index');
ylabel('Row Index');
axis tight; % Adjust axis limits to fit the data

%%
% 1. Generate a base map
map_size = [100, 150];
geo_map = Geometric_Map_Generator(2, map_size);

% 2. Generate a random route on the map
% We'll ask for a route with 200 points
route_points = 2000;
my_route = Geometric_Map_Route_Generator(geo_map, route_points);

% 3. Display the contents of the generated structure
disp('Generated Route Structure:');
disp(my_route);

% 4. Visualize the map and the route
figure;
% Display the map itself as a background image
imagesc(geo_map);
colormap(parula); % A nice looking colormap
colorbar;
hold on; % Keep the map image and plot on top of it

% Plot the route path. We plot the columns against the rows.
% 'path(:,2)' is the x-coordinate (columns)
% 'path(:,1)' is the y-coordinate (rows)
plot(my_route.path(:,2), my_route.path(:,1), ...
     'r-', 'LineWidth', 2); % A thick red line

% Add a marker for the start point
plot(my_route.path(1,2), my_route.path(1,1), ...
    'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
hold on
plot(my_route.path(route_points,2), my_route.path(route_points,1), ...
    'go', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
title(sprintf('Random Route with %d Steps on a %dx%d Map', ...
              route_points, map_size(1), map_size(2)));
xlabel('Column Index');
ylabel('Row Index');
legend('Route Path', 'Start Point', 'Location', 'northeast');
axis tight;
hold off;