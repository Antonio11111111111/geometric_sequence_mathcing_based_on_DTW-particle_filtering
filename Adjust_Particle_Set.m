function P_out = Adjust_Particle_Set(P_in, M_new)
    
    M_curr = size(P_in, 1); % Current particle count
    
    if M_new == M_curr 
        P_out = P_in; % No change needed
        return;
        
    elseif M_new > M_curr 
        % --- Add Particles ---
        N_add = M_new - M_curr;
        
        % Randomly sample with replacement
        indices = randi(M_curr, N_add, 1); 
        P_add = P_in(indices, :);
        
        % Add jitter to new particles for diversity
        jitter_pos = randn(N_add, 2) * 0.1; 
        jitter_ang = randn(N_add, 1) * 0.01;
        P_add(:, 1:2) = P_add(:, 1:2) + jitter_pos;
        P_add(:, 3) = P_add(:, 3) + jitter_ang;
        
        P_out = [P_in; P_add]; % Combine old and new
        
    else 
        % --- Remove Particles (M_new < M_curr) ---
        % Randomly sample without replacement
        indices = randperm(M_curr, M_new); 
        P_out = P_in(indices, :); % Keep only selected particles
    end
end