% =========================================================================
% File 5:  find_closest_contour_points.m 
% =========================================================================

function Y = find_closest_contour_points(P, M, X_map, Y_map, Z_map)

    N = size(P, 1);
    Y = zeros(N, 2);
    
    % --- Optimization Start ---
    % 1. Get all unique magnetic values needed
    M_unique = unique(M);
    
    % 2. Calculate all required contours at once
    %    (This is the *only* call to contourc in this function)
    C = contourc(X_map(1,:), Y_map(:,1), Z_map, M_unique);
    
    % 3. Parse all contours into a fast-lookup Map
    %    contour_map(magnetic_value) = [point_coordinates]
    contour_map = parse_contourc(C);
    % --- Optimization End ---

    for i = 1:N
        p_i = P(i, :);
        m_i = M(i);
        
        % 4. (Fast) Look up pre-calculated contour points from Map
        if ~isKey(contour_map, m_i)
            % If no contour was found for this value 
            % (e.g., map edge or noise), skip
            Y(i, :) = p_i;
            continue;
        end
        contour_points = contour_map(m_i);
    
        % 5. (Fast) Find the closest point on this contour
        distances_sq = sum((contour_points - p_i).^2, 2);
        [~, min_idx] = min(distances_sq);
        
        if isempty(min_idx)
           Y(i, :) = p_i; 
        else
           Y(i, :) = contour_points(min_idx, :);
        end
    end
end