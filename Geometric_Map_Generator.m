function geometric_map = Geometric_Map_Generator(D, size)
%GEOMETRIC_MAP_GENERATOR Generates a random geometric map.
%   map = GEOMETRIC_MAP_GENERATOR(D, size) generates a random map with
%   dimensions specified by D. The values in the map are uniformly
%   distributed random numbers in the interval (0,1).
%
%   Parameters:
%       D (dimension): A scalar value, must be 1 or 2.
%           - If D = 1, a 1-dimensional map (row vector) is generated.
%           - If D = 2, a 2-dimensional map (matrix) is generated.
%
%       size (size vector): A vector specifying the map dimensions.
%           - If D = 1, the length of the map is determined by size(1).
%           - If D = 2, the map size is [size(1), size(2)].
%
%   Output:
%       geometric_map: The randomly generated map (vector or matrix).
%
%   Examples:
%       % Generate a 1D map of length 100
%       map1D = geometric_map_generator(1, [100]);
%
%       % Generate a 2D map of size 50x75
%       map2D = geometric_map_generator(2, [50, 75]);

% This is the generation of the geometric map. 
% The design of this is below:
% name: geometric_map_generator
% param: D(dimension): 
%        for this particular project, D has only 1 or 2, with
%        1, only 1 dimensional map will be generated, and 2 dimensional map will
%        be generated with the parameter is 2.
%        size: 
%        input a vector [x, y], if the parameter D is 1, x is the
%        length while the length and the width will be determinated if D is
%        two.
% output: geometric_map, the shape will be decided with the parameter D,
%         and we will get a randomly generated geometric map.


% Use a switch statement to handle the different dimensions
switch D
    case 1
        % For 1D, generate a row vector of random numbers
        % The length is the first element of the size vector
        if numel(size) < 1
            error('For D=1, size_vec must have at least one element specifying the length.');
        end
        map_length = size(1);
        geometric_map = rand(1, map_length);
        
    case 2
        % For 2D, generate a matrix of random numbers
        % The dimensions are the first two elements of the size vector
        if numel(size) < 2
            error('For D=2, size_vec must have at least two elements for length and width.');
        end
        map_rows = size(1);
        map_cols = size(2);
        geometric_map = rand(map_rows, map_cols);
        
    otherwise
        % If D is not 1 or 2, throw an error
        error('Invalid dimension specified. D must be 1 or 2.');
end

end