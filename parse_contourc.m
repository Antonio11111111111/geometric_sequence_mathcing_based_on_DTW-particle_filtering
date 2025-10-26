% =========================================================================
% File 6:  parse_contourc.m 
% =========================================================================

function contour_map = parse_contourc(C)
    % This function now returns a containers.Map
    % Format: contour_map(level) = [ (Nx2) point coordinates ]
    
    contour_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
    if isempty(C)
        return;
    end
    
    idx = 1;
    while idx < size(C, 2)
        % 1. Read segment header
        level = C(1, idx); % Magnetic value for this contour
        N = C(2, idx);     % Number of points in this segment
        
        % 2. Determine coordinate indices for this segment
        start_col = idx + 1;
        end_col = idx + N;
        
        % 3. Extract coordinates
        x_coords = C(1, start_col:end_col);
        y_coords = C(2, start_col:end_col);
        
        segment_points = [x_coords', y_coords'];
        
        % 4. Add these points to the Map
        if isKey(contour_map, level)
            % If this magnetic value (key) already exists, append points
            contour_map(level) = [contour_map(level); segment_points];
        else
            % Otherwise, create a new entry
            contour_map(level) = segment_points;
        end
        
        % 5. Move to the next segment
        idx = end_col + 1;
    end
end