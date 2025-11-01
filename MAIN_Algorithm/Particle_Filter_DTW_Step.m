function [new_particles, best_guess] = Particle_Filter_DTW_Step(old_particles, live_sequence, ...
                                    geo_map, window_radius, process_noise, dtw_noise_std)
%
% Function: Performs one step of a Particle Filter (Predict, Update, Resample).
% Core:   Uses Dynamic Time Warping (DTW) distance as the measurement model
%         to calculate particle weights.
% Summary: This function is primarily intended for indoor localization based on
%          geomagnetic sequences.
%          "Particles" represent possible locations (indices) on the complete
%          geomagnetic map (geo_map).
%          We compare the "real-time sequence from the sensor" with the
%          "sequence corresponding to a particle's location on the map".
%          DTW is used to calculate their similarity, which assigns a
%          "weight" (likelihood) to that particle.
%
% -------------------------------------------------------------------------
% param:
%   old_particles:   [Nx1] vector. The N particle locations (indices on the map)
%                    from the previous time step.
%   live_sequence:   [1xK] vector. The "live" geomagnetic sequence collected
%                    by the sensor.
%   geo_map:         [1xM] vector. The complete geomagnetic map database
%                    (a very long sequence).
%   window_radius:   Scalar. The radius of the window used to extract
%                    subsequences from the map.
%                    The extracted map sequence length will be 2*window_radius + 1.
%                    This length (K') must be compatible with the
%                    live_sequence length (K).
%   process_noise:   Scalar. Standard deviation of the process noise (motion noise).
%                    Used in the prediction step.
%   dtw_noise_std:   Scalar. Standard deviation of the measurement noise.
%                    Used to convert DTW distance into a weight.
%                    (Weight formula: exp(-dist^2 / (2*std^2)) ).
%
% -------------------------------------------------------------------------
% Outputs:
%   new_particles: [Nx1] vector. The new generation of particle locations
%                  after resampling.
%   best_guess:    Scalar. The best estimate of the current location
%                  (mean of the new particles).
% -------------------------------------------------------------------------
%
% This is the core of your solution: Particle Filter combined with DTW
%
%   old_particles:   [Nx1] Particles from the previous step
%   live_sequence:   [1xK] The "real" geomagnetic sequence (from sensor)
%   geo_map:         [1xM] The complete geomagnetic map
%   window_radius:   Window radius (e.g., 2)
%   process_noise:   Motion noise
%   dtw_noise_std:   Noise variance for DTW distance (for Gaussian function)
%
    N = length(old_particles);
    map_len = length(geo_map);
    window_len = 2*window_radius + 1;
    
    % --- 1. PREDICT ---
    % "Jitter" all particles (using a simplified -1, 0, +1 motion model)
    % The process_noise determines the magnitude of the jitter.
    noise = round(randn(N, 1) * process_noise);
    predicted_particles = old_particles + noise;
    
    % Boundary check: Ensure particles do not fall off the "measurable"
    % range of the map.
    % Since we extract from (pos - r) to (pos + r), pos must be
    % >= r+1 and <= map_len - r.
    predicted_particles(predicted_particles <= window_radius) = window_radius + 1;
    predicted_particles(predicted_particles >= (map_len - window_radius)) = map_len - window_radius - 1;
    
    % --- 2. UPDATE (Update Weights) ---
    weights = zeros(N, 1); % Store the weight for each particle
    
    for i = 1:N
        % a. Get the predicted location of the i-th particle
        particle_loc = predicted_particles(i);
        
        % b. Extract the sequence this particle "should" see from the "map"
        %    Location is particle_loc, radius is window_radius
        map_sequence_indices = (particle_loc - window_radius) : (particle_loc + window_radius);
        map_sequence = geo_map(map_sequence_indices);
        
        % c. [Core] Calculate the DTW distance between the "real sequence"
        %    and the "map sequence".
        %    The smaller the dist, the more likely this particle is correct.
        dist = simple_global_dtw(live_sequence, map_sequence);
        
        % d. [Core] Use a Gaussian function to convert "distance" to "weight"
        %    Smaller distance -> higher weight (closer to 1)
        %    Larger distance -> lower weight (closer to 0)
        weights(i) = exp( -(dist^2) / (2 * dtw_noise_std^2) );
    end
    
    % Prevent all-zero weights (e.g., if sensor data is completely off-map)
    if sum(weights) == 0
        % If all particles are lost, reset weights to a uniform distribution
        % to avoid NaN.
        weights = ones(N, 1); 
    end
    % Normalize weights
    weights = weights / sum(weights);

    % --- 3. RESAMPLE ---
    % "Kill" low-weight particles, "duplicate" high-weight particles
    % (Using standard systematic resampling)
    new_particles = zeros(N, 1);
    c = cumsum(weights);         % Calculate the cumulative distribution
    u = (0:N-1)/N + rand()/N;    % Generate N evenly spaced random points
    
    i_ptr = 1;
    for j = 1:N
        while u(j) > c(i_ptr)
            i_ptr = i_ptr + 1;
        end
        % Select the i_ptr-th particle
        new_particles(j) = predicted_particles(i_ptr);
    end
    
    % --- 4. ESTIMATE ---
    % Our best guess is the mean position of all new particles
    best_guess = mean(new_particles);
end


%% --- Usage Example ---
%{
    % To run this example, copy the function above and all the code below
    % into the same .m file, then run.
    
    clc; clear; close all;
    
    % --- 1. Define Simulated Data ---
    
    % a. Create a long "geomagnetic map" (simulated with a noisy sine wave)
    map_len = 500;
    geo_map = sin( (1:map_len) / 20) * 10 + randn(1, map_len) * 0.5 + 50;
    
    % b. Define particle and window parameters
    N = 1000;            % Number of particles
    window_radius = 10;  % Window radius
    window_len = 2*window_radius + 1; % Total window length (21)
    
    % c. Define noise parameters
    process_noise = 2.0;   % Motion noise (moves 1 step, but with 2.0 std dev jitter)
    dtw_noise_std = 5.0;   % Measurement noise (how quickly weights drop off with DTW distance)
    
    % d. Initialize particles (assume we don't know where we are,
    %    distribute uniformly)
    current_particles = round(rand(N, 1) * (map_len - 2*window_radius) + window_radius + 1);
    
    % e. Define a "true" path
    %    Assume the "true" location starts at 50 and moves 1 step forward each time
    true_path = 50:1:150;
    
    % --- 2. Simulate the Filtering Process ---
    estimations = zeros(size(true_path)); % Store the estimate at each step
    
    fprintf('Starting simulation...\n');
    
    for t = 1:length(true_path)
        
        % a. Get the current "true location"
        current_true_loc = true_path(t);
        
        % b. "Sensor" collects the "live sequence"
        %    (Extracted from the map at the true location, plus a little noise)
        live_indices = (current_true_loc - window_radius) : (current_true_loc + window_radius);
        live_sequence = geo_map(live_indices) + randn(1, window_len) * 0.2; % Add sensor noise
        
        % c. [Key Step] Call the particle filter function
        [new_particles, best_guess] = Particle_Filter_DTW_Step(...
                                            current_particles, ...
                                            live_sequence, ...
                                            geo_map, ...
                                            window_radius, ...
                                            process_noise, ...
                                            dtw_noise_std);
                                        
        % d. Save results and update particles
        estimations(t) = best_guess;
        current_particles = new_particles;
        
        if mod(t, 20) == 0
            fprintf('Time step %d: True Loc = %.1f, Estimated Loc = %.1f\n', t, current_true_loc, best_guess);
        end
    end
    
    fprintf('Simulation finished.\n');
    
    % --- 3. Visualize Results ---
    figure;
    hold on;
    grid on;
    
    % Plot the geomagnetic map
    plot(1:map_len, geo_map, 'c-', 'LineWidth', 1, 'DisplayName', 'Geomagnetic Map (Geo Map)');
    
    % Plot the true path
    plot(true_path, true_path, 'k-', 'LineWidth', 2, 'DisplayName', 'True Path');
    
    % Plot the particle filter's estimated path
    plot(true_path, estimations, 'r--', 'LineWidth', 2, 'DisplayName', 'PF Estimate');
    
    % Plot the final particle distribution
    y_vals = max(geo_map) * ones(N, 1);
    scatter(current_particles, y_vals, 15, 'b', 'filled', 'MarkerFaceAlpha', 0.2, 'DisplayName', 'Final Particle Distribution');
    
    xlabel('Map Position (Index)');
    ylabel('Geomagnetic Strength / Position');
    title('Particle Filter (PF) + DTW Localization Simulation');
    legend('show');
    hold off;
%}


%% -----------------------------------------------------------------------
%                 DTW Helper Function
% -----------------------------------------------------------------------
function [dist] = simple_global_dtw(A, B)
% A very basic implementation of global DTW (no windowing/constraints)
% A and B must be row or column vectors

    n = length(A);
    m = length(B);
    
    % Cost Matrix
    D = zeros(n, m);
    
    % Calculate the Euclidean distance (here, 1D, so just squared difference)
    % between every pair of points in the two sequences.
    cost_matrix = zeros(n, m);
    for i = 1:n
        for j = 1:m
            cost_matrix(i, j) = (A(i) - B(j))^2;
        end
    end

    % --- Dynamic Programming ---
    
    % Initialize the first element
    D(1, 1) = cost_matrix(1, 1);
    
    % Initialize the first row
    for i = 2:n
        D(i, 1) = cost_matrix(i, 1) + D(i-1, 1);
    end
    
    % Initialize the first column
    for j = 2:m
        D(1, j) = cost_matrix(1, j) + D(1, j-1);
    end
    
    % Fill the rest of the cost matrix
    for i = 2:n
        for j = 2:m
            min_prev_cost = min([D(i-1, j),    % Top
                                 D(i, j-1),    % Left
                                 D(i-1, j-1)]); % Diagonal
                                 
            D(i, j) = cost_matrix(i, j) + min_prev_cost;
        end
    end
    
    % Final distance (usually sqrt is taken to return to original magnitude)
    dist = sqrt(D(n, m));
end