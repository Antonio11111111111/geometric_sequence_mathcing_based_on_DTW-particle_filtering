%% there is no need to have a function, but here is a possible example of using the DTW function 
%  the matching result itself is not well, conversely, it is BAD
clc
clear
geo_map = Geometric_Map_Generator(1, [100, 0]);
route = Geometric_Map_Route_Generator(geo_map, 10);

[~, loc] = Geomectric_Map_Matching_DTW_1D(geo_map, route.intensity',true); % Always remember that the route.intensity should be transformed

disp(loc)
disp(route.path)


