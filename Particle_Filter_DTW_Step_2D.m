function [particles_out, best_guess, d_min] = Particle_Filter_DTW_Step_2D(particles_in, live_sequence, ...
                                            pdr_history, geo_map, process_noise, DTW_NOISE_STD)
%
% Function: Particle_Filter_DTW_Step_2D
%
% Summary: Performs one step of a 2D Particle Filter (Propagation, Weighting, Resampling).
%          This function is designed for PDR (Pedestrian Dead Reckoning) aided
%          by geomagnetic field matching.
%
% Weighting: Particle weights are determined based on the Dynamic Time Warping (DTW)
%            distance between the sensor's live sequence and the sequence
%            reconstructed from the geomagnetic map along the particle's path.
%
% Note: This function requires the Signal Processing Toolbox (for the 'dtw' function).
%
% -------------------------------------------------------------------------
% Inputs:
%   particles_in:   [Mx3] matrix. The set of M particles from the previous step.
%                   Each row is [x, y, theta], where (x, y) is the position
%                   and theta is the heading (in radians).
%
%   live_sequence:  [Lx1] or [1xL] vector. The "live" geomagnetic sequence
%                   of length L collected by the sensor, corresponding to the
%                   steps in pdr_history.
%
%   pdr_history:    [Lx2] matrix. A history of the last L PDR steps.
%                   Each row is [step_length, d_theta], where d_theta is the
%                   change in heading (radians) *for that step*.
%
%   geo_map:        Struct. A structure containing the geomagnetic map data:
%                   .X_grid:  [RxC] X-coordinates grid (from meshgrid)
%                   .Y_grid:  [RxC] Y-coordinates grid (from meshgrid)
%                   .Mag_map: [RxC] Magnetic field strength at each (X, Y) point
%
%   process_noise:  Struct. Standard deviations for the motion model noise.
%                   .step_std:  Std dev for step length noise.
%                   .theta_std: Std dev for heading change (d_theta) noise.
%
%   DTW_NOISE_STD:  Scalar. The standard deviation of the measurement noise.
%                   Used in the Gaussian kernel to convert DTW distance
%                   to a particle weight. (e.g., exp(-dist^2 / (2*STD^2))).
%
% -------------------------------------------------------------------------
% Outputs:
%   particles_out: [Mx3] matrix. The new set of M particles after resampling.
%                  Format is [x, y, theta].
%
%   best_guess:    [1x3] vector. The best estimate of the current state.
%                  Format is [mean_x, mean_y, mean_theta].
%                  The heading (theta) is calculated using a circular mean.
%
%   d_min:         Scalar. The minimum DTW distance found among all
%                  particles in this step. Useful for monitoring filter health.
%
% -------------------------------------------------------------------------
% Usage Example (Conceptual):
% -------------------------------------------------------------------------
%{
    % This is a conceptual example. You need real data for it to run.
    % To use this, save the function above in its own .m file named
    % "Particle_Filter_DTW_Step_2D.m" and run this script from a
    % separate file or the command window.
    
    clc; clear;
    
    % --- 1. Setup ---
    M = 1000; % Number of particles
    L = 20;   % Length of history/sequence
    
    % a. Create a dummy geomagnetic map (e.g., 100m x 100m)
    [X, Y] = meshgrid(1:100, 1:100);
    Z = sin(X/10) + cos(Y/10); % A simple magnetic map
    geo_map.X_grid = X;
    geo_map.Y_grid = Y;
    geo_map.Mag_map = Z;
    
    % b. Initialize particles (e.g., centered around (50, 50))
    particles_in = [50 + randn(M,1), 50 + randn(M,1), randn(M,1) * 0.1];
    
    % c. Define noise parameters
    process_noise.step_std = 0.1;  % 10cm step noise
    process_noise.theta_std = deg2rad(3); % 3 degrees heading noise
    DTW_NOISE_STD = 5.0; % DTW distance noise
    
    % --- 2. Simulation Step ---
    
    % Assume at this step, we received:
    % i. A "live" magnetic sequence
    live_sequence = Z(50, 50:69) + randn(1, L)*0.1; % Dummy sequence
    
    % ii. A history of PDR steps that led to this sequence
    %    (L steps, each is [step_length, d_theta])
    pdr_history = [repmat(0.7, L, 1), zeros(L, 1)]; % Dummy: 20 steps, 0.7m, 0 heading change
    pdr_history(end, 2) = deg2rad(1); % Last step had a 1-degree turn
    
    fprintf('Running one step of 2D Particle Filter...\n');
    
    % --- 3. Call the Function ---
    % Make sure Particle_Filter_DTW_Step_2D.m is in your MATLAB path
    [particles_out, best_guess, d_min] = Particle_Filter_DTW_Step_2D(...
                                            particles_in, ...
                                            live_sequence, ...
                                            pdr_history, ...
                                            geo_map, ...
                                            process_noise, ...
                                            DTW_NOISE_STD);
    
    fprintf('Step complete.\n');
    fprintf('Best Guess (X, Y, Theta): %.2f, %.2f, %.2f rad\n', ...
            best_guess(1), best_guess(2), best_guess(3));
    fprintf('Minimum DTW distance: %.4f\n', d_min);
    
    % In a real application, you would loop:
    % particles_in = particles_out;
    % and get new live_sequence and pdr_history for the next time step.
    
    % Visualize (optional)
    figure;
    hold on;
    % Draw the map
    contourf(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, 20);
    colormap(jet);
    colorbar;
    
    % Draw new particle cloud
    scatter(particles_out(:, 1), particles_out(:, 2), 10, 'r', 'filled', 'MarkerFaceAlpha', 0.3);
    
    % Draw best guess
    plot(best_guess(1), best_guess(2), 'gx', 'MarkerSize', 15, 'LineWidth', 3);
    
    axis equal;
    title('Particle Filter State After One Step');
    xlabel('X Coordinate (m)');
    ylabel('Y Coordinate (m)');
    legend('Magnetic Map', 'Particle Cloud (Out)', 'Best Guess');
%}
% -------------------------------------------------------------------------
    
    M = size(particles_in, 1);
    L = size(pdr_history, 1);
    
    % --- 1. Propagation (Prediction) ---
    % Propagate each particle forward by one step based on the *last*
    % PDR measurement, adding process noise.
    
    last_pdr_step = pdr_history(end, :); % Get the most recent [step_len, d_theta]
    
    % Add noise to the PDR measurement for each particle
    step_noise = randn(M, 1) * process_noise.step_std;
    theta_noise = randn(M, 1) * process_noise.theta_std;
    
    steps = last_pdr_step(1) + step_noise;
    d_thetas = last_pdr_step(2) + theta_noise;
    
    % Get old states
    old_x = particles_in(:, 1);
    old_y = particles_in(:, 2);
    old_theta = particles_in(:, 3);
    
    % Calculate new states
    new_theta = old_theta + d_thetas;
    new_x = old_x + sin(new_theta) .* steps; % Note: sin for X
    new_y = old_y + cos(new_theta) .* steps; % Note: cos for Y (assuming 0 rad = North/Y-axis)
    
    % Boundary check: Keep particles within the map limits
    new_x = max(1, min(geo_map.X_grid(1,end), new_x));
    new_y = max(1, min(geo_map.Y_grid(end,1), new_y));
    
    propagated_particles = [new_x, new_y, new_theta];

    % --- 2. Weighting ---
    weights = zeros(M, 1);
    dtw_variance = DTW_NOISE_STD^2;
    
    %% [Modification] Add an array to store all DTW distances
    all_distances = zeros(M, 1); 
    
    for m = 1:M
        % --- 2a. Path Reconstruction ---
        % For each particle, reconstruct its *past path* (L steps long)
        % by back-propagating from its new position using the PDR history.
        
        particle_path = zeros(L, 2); % To store [x, y] coords for this particle's path
        
        % Start from the particle's *newly propagated* position
        current_pos = propagated_particles(m, 1:2);
        current_theta = propagated_particles(m, 3);
        
        particle_path(end, :) = current_pos; % The last point in the path is the current position
        
        % Work backwards from L-1 down to 1
        for k = L-1:-1:1
            % Get the PDR step that *led to* the (k+1)-th position
            pdr_step_len = pdr_history(k+1, 1);
            pdr_d_theta = pdr_history(k+1, 2);
            
            % Revert the state update
            prev_theta = current_theta - pdr_d_theta;
            prev_x = current_pos(1) - sin(current_theta) * pdr_step_len;
            prev_y = current_pos(2) - cos(current_theta) * pdr_step_len;
            
            particle_path(k, :) = [prev_x, prev_y];
            
            % Update for next iteration
            current_pos = [prev_x, prev_y];
            current_theta = prev_theta;
        end
        
        % --- 2b. Generate Map Sequence ---
        % Sample the geomagnetic map along the reconstructed particle path
        map_sequence = interp2(geo_map.X_grid, geo_map.Y_grid, geo_map.Mag_map, ...
                              particle_path(:, 1), particle_path(:, 2), 'linear', 0);
                          
        % --- 2c. Calculate DTW Distance and Convert to Weight ---
        % Compare the sensor's sequence with the map's sequence
        distance = dtw(live_sequence(:), map_sequence(:));
        
        %% [Modification] Store this particle's distance
        all_distances(m) = distance;
        
        % Use a Gaussian kernel to convert distance to weight
        weights(m) = exp(-distance^2 / (2 * dtw_variance));
    end
    
    % --- 2d. Normalize Weights ---
    sum_weights = sum(weights);
    if sum_weights > 1e-15 % Check for numerical stability
        weights = weights / sum_weights;
    else
        % All particles have zero weight (lost track)
        % Reset to uniform distribution to prevent failure
        weights = ones(M, 1) / M;
    end
    
    % --- 3. Resampling ---
    % Use systematic resampling to "kill" low-weight particles and
    % "duplicate" high-weight particles.
    
    particles_out = zeros(M, 3);
    cdf = cumsum(weights); % Cumulative Distribution Function
    r_0 = rand() / M; % Initial random offset
    idx_j = 1;
    
    for m = 1:M
        U = r_0 + (m-1)/M; % Get the m-th sample point
        while U > cdf(idx_j)
            idx_j = idx_j + 1; % Find the particle index this sample corresponds to
        end
        particles_out(m, :) = propagated_particles(idx_j, :);
    end
    
    % --- 4. Estimation ---
    % Calculate the best guess (state estimate) from the new particle set.
    
    % Simple mean for position
    best_guess(1:2) = mean(particles_out(:, 1:2), 1);
    
    % Circular mean for heading (theta) to handle 0/2*pi wraparound
    best_guess(3) = atan2(mean(sin(particles_out(:, 3))), mean(cos(particles_out(:, 3))));
    
    d_min = min(all_distances);
end