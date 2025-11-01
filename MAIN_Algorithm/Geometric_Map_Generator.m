function geometric_map = Geometric_Map_Generator(D, map_size)

switch D
    case 1
        if numel(map_size) < 1
            error('For D=1, size_vec must have at least one element specifying the length.');
        end
        map_length = map_size(1);
        geometric_map = rand(1, map_length);
        
    case 2
        if numel(map_size) < 2
            error('For D=2, size_vec must have at least two elements for length and width.');
        end
        map_rows = map_size(1);
        map_cols = map_size(2);
        geometric_map = rand(map_rows, map_cols);
        
    otherwise
        error('Invalid dimension specified. D must be 1 or 2.');
end
end
