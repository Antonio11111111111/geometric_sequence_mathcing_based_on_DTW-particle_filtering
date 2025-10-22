function [new_particles, best_guess] = Particle_Filter_DTW_Step_2D(old_particles, live_sequence, pdr_history, ...
                                    geo_map, process_noise, dtw_noise_std)
% 扩展到2D的粒子滤波 + DTW 步骤

    N = size(old_particles, 1);
    K = length(live_sequence); 
    
    predicted_particles = zeros(N, 3);
    weights = zeros(N, 1); 

    % --- 1. PREDICT (预测) ---
    last_pdr_step = pdr_history(end, :); 
    
    for i = 1:N
        old_x = old_particles(i, 1);
        old_y = old_particles(i, 2);
        old_theta = old_particles(i, 3);
        
        % 应用 "粒子噪声" (process_noise)
        noisy_step = last_pdr_step(1) + randn() * process_noise.step_std;
        if noisy_step < 0, noisy_step = 0; end 
        
        noisy_d_theta = last_pdr_step(2) + randn() * process_noise.theta_std;
        new_theta = old_theta + noisy_d_theta;
        
        new_x = old_x + sin(new_theta) * noisy_step;
        new_y = old_y + cos(new_theta) * noisy_step;
        
        predicted_particles(i, :) = [new_x, new_y, new_theta];
    end
    
    map_X_max = geo_map.X_grid(1, end);
    map_Y_max = geo_map.Y_grid(end, 1);
    predicted_particles(:, 1) = max(1, min(map_X_max, predicted_particles(:, 1)));
    predicted_particles(:, 2) = max(1, min(map_Y_max, predicted_particles(:, 2)));

    % --- 2. UPDATE (更新权重) ---
    for i = 1:N
        particle_state = predicted_particles(i, :); 
        map_sequence = zeros(1, K);
        current_pos = particle_state(1:2);
        current_theta = particle_state(3);
        
        for k = K:-1:1
            map_sequence(k) = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                                      current_pos(1), current_pos(2), 'linear', 0); 
                                  
            step_k = pdr_history(k, 1);
            d_theta_k = pdr_history(k, 2);
            prev_theta = current_theta - d_theta_k;
            
            prev_pos_x = current_pos(1) - sin(current_theta) * step_k;
            prev_pos_y = current_pos(2) - cos(current_theta) * step_k;
            
            current_pos = [prev_pos_x, prev_pos_y];
            current_theta = prev_theta;
        end
        
        dist = simple_global_dtw(live_sequence, map_sequence);
        % 使用更新后的 dtw_noise_std
        weights(i) = exp( -(dist^2) / (2 * dtw_noise_std^2) );
    end
    
    if sum(weights) == 0
        weights = ones(N, 1); 
    end
    weights = weights / sum(weights); 

    % --- 3. RESAMPLE (重采样) ---
    new_particles = zeros(N, 3);
    c = cumsum(weights);
    u = (0:N-1)/N + rand()/N; 
    
    i_ptr = 1;
    for j = 1:N
        while u(j) > c(i_ptr)
            i_ptr = i_ptr + 1;
        end
        new_particles(j, :) = predicted_particles(i_ptr, :);
    end

    % --- 4. ESTIMATE (估计) ---
    best_guess_x = sum(weights .* predicted_particles(:, 1));
    best_guess_y = sum(weights .* predicted_particles(:, 2));
    best_guess = [best_guess_x, best_guess_y]; 
end

