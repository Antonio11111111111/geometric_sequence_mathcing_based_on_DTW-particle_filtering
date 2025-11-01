clc
clear
close all

%% Test for the function of 1D classical dtw on 2 sequences(without figure output)
geo_map = [1.2, 1.1, 1.3, 2.0, 5.0, 8.0, 9.0, 8.5, 4.0, 2.1, 1.4, 1.3, 1.5];
live_reading = [5.2, 7.8, 9.1, 9.0, 8.2, 3.5];
[dist, loc] = Geomectric_Map_Matching_DTW_1D(geo_map, live_reading);
fprintf("The minimum distanse is %d . \n", dist);
fprintf('The target sequence on the geo map is:\n');
disp(loc);
%% Test for the function of 1D classical dtw on 2 sequences(with figure output)
geo_map = [1.2, 1.1, 1.3, 2.0, 5.0, 8.0, 9.0, 8.5, 4.0, 2.1, 1.4, 1.3, 1.5];
live_reading = [5.2, 7.8, 9.1, 9.0, 8.2, 3.5];
[dist, loc] = Geomectric_Map_Matching_DTW_1D(geo_map, live_reading, true);