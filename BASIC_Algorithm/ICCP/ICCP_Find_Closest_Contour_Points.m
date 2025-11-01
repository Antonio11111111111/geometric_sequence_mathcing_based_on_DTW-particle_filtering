function Y = ICCP_Find_Closest_Contour_Points(P, M, X_map, Y_map, Z_map)
% Name: ICCP_Find_Closest_Contour_Points
% Description: Finds the closest point (Y) on a geomagnetic contour map
%              for each point in a track (P), given the measured magnetic
%              value (M) at that point.
% Parameters:
%   P: (N x 2) Input track points.
%   M: (N x 1) Measured magnetic values for each point in P.
%   X_map, Y_map, Z_map: Geomagnetic reference map.
% Output:
%   Y: (N x 2) The set of closest points on the contours.
% Example:
%   Y_targets = ICCP_Find_Closest_Contour_Points(P_current, M_measured, X, Y, Z);

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
contour_map = ICCP_Parse_Contourc(C);
% --- Optimization End ---

for i = 1:N
    p_i = P(i, :);
    m_i = M(i);
    
    % 4. (Fast) Look up pre-calculated contour points from Map
    if ~isKey(contour_map, m_i)
        % If no contour was found for this value 
        % (e.g., map edge or noise), use the point itself as target
        Y(i, :) = p_i;
        continue;
    end
    contour_points = contour_map(m_i);

    % 5. (Fast) Find the closest point on this contour
    % Calculate squared Euclidean distances
    distances_sq = sum((contour_points - p_i).^2, 2);
    [~, min_idx] = min(distances_sq);
    
    if isempty(min_idx)
       % Should not happen if key exists, but as a fallback
       Y(i, :) = p_i; 
    else
       % Assign the closest point
       Y(i, :) = contour_points(min_idx, :);
    end
end

end
