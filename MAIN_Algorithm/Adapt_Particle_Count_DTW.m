function M_new = Adapt_Particle_Count_DTW(M_curr, dist, Thresh_Low, Thresh_High, M_min, M_max)
    
    if dist > Thresh_High 
        % Error is high (bad match), increase particles
        M_new = min(M_curr * 2, M_max);
    elseif dist < Thresh_Low 
        % Error is low (good match), decrease particles
        M_new = max(M_curr / 2, M_min);
    else
        % Error is acceptable, no change
        M_new = M_curr;
    end
    
    M_new = round(M_new); % Ensure integer count
end
