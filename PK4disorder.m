% PURPOSE: Take track data and determine the disorder for each particle in
% each frame. Disorder is determined by variance of nearest neighbor
% distances for each particle. 
% 
% calculate the average disorder for each frame
% import track data
% use kd-tree to find neighbors of track data at different timesteps
% in a some radius around a single particle, define neighbors, calculate the distance between the neighbors
% then calculate the variance in this distance
% store this infor for each and every particle
% average all the variances over the whole system (i.e. for every particle)
% normalize by max variationn defined in bottom section

% INPUT: [frame, x position, y position]

% Load data
data = C1705Pa035mAFov14x2mm;
trajectory = data{:,1};
frame = data{:,2};
x_pos = data{:,3};
y_pos = data{:,4};

% Conversion: pixel to micrometers
dx = 14.20;  % Pixel size in micrometers

x_pos = x_pos * dx;
y_pos = y_pos * dx;

max_frame = max(frame);

% Initialize
positions_by_frame = cell(max_frame + 1, 1);
trajectory_counts = zeros(max_frame + 1, 1);

for i = 0:max_frame
    idx = (frame == i);
    positions = [trajectory(idx), x_pos(idx), y_pos(idx)];
    positions_by_frame{i + 1} = positions;

    trajectory_counts(i + 1) = numel(unique(positions(:, 1)));
end

% ---- Compute Averages ----
avg_particles_per_frame = mean(trajectory_counts);

% Combine all positions into one array for spacing estimate
all_coords = [x_pos, y_pos];
N = size(all_coords, 1);

% Compute nearest neighbor distances in X and Y
Mdl = KDTreeSearcher(all_coords);
[idxs, dists] = knnsearch(Mdl, all_coords, 'K', 2);  % 2: self and nearest neighbor
nearest_neighbors = all_coords(idxs(:,2), :);  % second column = nearest neighbor

dxs = abs(all_coords(:,1) - nearest_neighbors(:,1));
dys = abs(all_coords(:,2) - nearest_neighbors(:,2));

avg_dx_um = mean(dxs);
avg_dy_um = mean(dys);

% ---- Print Summary ----
fprintf('--- Dusty Plasma Frame Statistics ---\n');
fprintf('Average number of particles per frame: %d\n', round(avg_particles_per_frame));


sample_frame = positions_by_frame{100};  % or any non-empty frame
x_pixels = sample_frame(:,2);
y_pixels = sample_frame(:,3);

frame_width_pixels = max(x_pixels) - min(x_pixels);
frame_height_pixels = max(y_pixels) - min(y_pixels);

fprintf('Frame size: %.2f pixels wide × %.2f pixels tall\n', ...
        frame_width_pixels/dx, frame_height_pixels/dx);


% OUTPUT: {frame}[x position, y position, disorder]
% 
% disorder value calculated as variance of displacement around a
% particle for each particle. This will give you a map of all of the
% disorder spatially around each particle for each frame.

% Plot disorder over all positions at some frame or animation or 2D
% histogram plot

% total/average disorder of all particle at each frame and a plot of total
% disorder as it changes with time/frames

% total disorder over all frames 


% Create KD-Tree for Fast Neighbor Search
% KD-Tree Overview:
% -----------------
% The KD-Tree (K-dimensional tree) is a spatial data structure that 
% efficiently organizes points in multi-dimensional space for fast 
% nearest-neighbor searches.

% What KD-Tree Does:
% ------------------
% 1. It **divides the data space into hierarchical partitions** to 
%    allow efficient searches instead of brute-force comparisons.
% 2. It **speeds up nearest-neighbor queries** significantly compared to 
%    computing pairwise distances for all points.
% 3. It **supports fast searching operations**, such as:
%    - **knnsearch**: Finds the k-nearest neighbors to a query point.
%    - **rangesearch**: Finds all points within a given radius of a query point.

% MATLAB's KD-Tree Tool:
% ----------------------
% MATLAB provides **KDTreeSearcher**, which creates a KD-tree from data.
% Once the tree is built, you can perform fast nearest-neighbor searches.

% KD-Tree Functions:
% ------------------
% - `KDTreeSearcher(data)`: Creates a KD-tree model from input data.
% - `knnsearch(Mdl, query_point)`: Finds the nearest neighbor to the query point.
% - `rangesearch(Mdl, query_point, radius)`: Finds all points within a given radius.

%% Uniform Max Disorder with Multiple Realizations 

% ---- Parameters ----
num_realizations = 20;
num_particles = round(avg_particles_per_frame);
samples_per_particle = 60;
alpha = 0.6;
window_radius = 30*dx;  % in µm

% Frame size from real data (already in µm)
sample_frame = positions_by_frame{100};
x_real = sample_frame(:,2);
y_real = sample_frame(:,3);
frame_size = [max(x_real) - min(x_real), max(y_real) - min(y_real)];

% Store disorder averages
all_disorder_means = zeros(num_realizations, 1);

for rep = 1:num_realizations
    % --- Generate Uniform Random Particle Positions ---
    x = rand(num_particles, 1) * frame_size(1);
    y = rand(num_particles, 1) * frame_size(2);
    coords = [x, y];
    Mdl = KDTreeSearcher(coords);

    % --- Density-Based Disorder ---
    density_disorder = zeros(num_particles, 1);
    for i = 1:num_particles
        pt = coords(i,:);
        N_counts = zeros(samples_per_particle, 1);
        for s = 1:samples_per_particle
            angle = 2 * pi * rand;
            r = window_radius * sqrt(rand);
            cx = pt(1) + r * cos(angle);
            cy = pt(2) + r * sin(angle);
            if cx < 0 || cx > frame_size(1) || cy < 0 || cy > frame_size(2)
                N_counts(s) = NaN; continue;
            end
            dists = sqrt((coords(:,1) - cx).^2 + (coords(:,2) - cy).^2);
            N_counts(s) = sum(dists <= window_radius) - 1;
        end
        valid = ~isnan(N_counts);
        if sum(valid) > 1
            density_disorder(i) = var(N_counts(valid));
        else
            density_disorder(i) = NaN;
        end
    end

    % --- Structure-Based Disorder (G6) ---
    structure_disorder = zeros(num_particles, 1);
    for i = 1:num_particles
        pt = coords(i,:);
        [idxs, ~] = knnsearch(Mdl, pt, 'K', 7);  % self + 6 neighbors
        neighbor_coords = coords(idxs(2:end), :);
        angles = atan2(neighbor_coords(:,2) - pt(2), neighbor_coords(:,1) - pt(1));
        g6 = mean(exp(1i * 6 * angles));
        structure_disorder(i) = 1 - abs(g6);  % Invert to get disorder
    end

    % --- Combine and Store ---
    combined_disorder = alpha * density_disorder + (1 - alpha) * structure_disorder;
    all_disorder_max(rep) = max(combined_disorder);

    % Store last realization for plot
    if rep == num_realizations
        last_coords = coords;
        last_combined_disorder = combined_disorder;
    end
end

% --- Output Results ---
avg_combined_disorder = mean(all_disorder_max);
std_combined_disorder = std(all_disorder_max);

fprintf('--- Max Disorder (Uniform) ---\n');
fprintf('Mean of all Max disorders over %d runs: %.4f\n', num_realizations, avg_combined_disorder);
fprintf('Standard deviation: %.4f\n', std_combined_disorder);

% --- Plot Last Realization ---
figure('Units', 'normalized', 'Position', [0.1, 0.2, 0.9, 0.6]);
scatter(last_coords(:,1), last_coords(:,2), 15, last_combined_disorder, 'filled');
colormap(parula); 
colorbar;
colorbar_handle = colorbar;
colorbar_handle.Label.String = 'Disorder';
colorbar_handle.Label.FontSize = 18;
title(sprintf('Uniform Max Disorder Sample (\\alpha = %.2f), Rep #%d', alpha, num_realizations), 'FontSize', 24);
xlabel('X (µm)', 'FontSize', 18); 
ylabel('Y (µm)', 'FontSize', 18);
set(gca, 'FontSize', 18);
axis equal; 
grid on;

  

%% Generate Ideal Hexagonal Grid for Minimal Disorder Case

% Use frame size from a real frame
sample_frame = positions_by_frame{100};
x_real = sample_frame(:,2);
y_real = sample_frame(:,3);
max_x = max(x_real);
max_y = max(y_real);
min_x = min(x_real);
min_y = min(y_real);
frame_size = [max_x - min_x, max_y - min_y];

% Parameters
spacing = 20 * dx;  % horizontal pixel spacing
row_shift = spacing / 2;
vertical_spacing = spacing * sqrt(3) / 2;

% Generate grid
x_coords = 0:spacing:frame_size(1);
y_coords = 0:vertical_spacing:frame_size(2);
[X, Y] = meshgrid(x_coords, y_coords);

% Apply hexagonal shift: shift every other row
for i = 1:size(Y,1)
    if mod(i,2)==0
        X(i,:) = X(i,:) + row_shift;
    end
end

% Flatten and trim to ~450 particles
x = X(:);
y = Y(:);
coords = [x, y];
if size(coords,1) > 450
    coords = coords(1:450,:);
end
num_particles = size(coords,1);
Mdl = KDTreeSearcher(coords);

% -------- Density-Based Disorder --------
window_radius = 30*dx;
samples_per_particle = 60;
density_disorder = zeros(num_particles, 1);

for i = 1:num_particles
    pt = coords(i,:);
    N_counts = zeros(samples_per_particle, 1);
    for s = 1:samples_per_particle
        angle = 2 * pi * rand;
        r = window_radius * sqrt(rand);
        cx = pt(1) + r * cos(angle);
        cy = pt(2) + r * sin(angle);

        if cx < 0 || cx > frame_size(1) || cy < 0 || cy > frame_size(2)
            N_counts(s) = NaN;
            continue;
        end

        dists = sqrt((coords(:,1) - cx).^2 + (coords(:,2) - cy).^2);
        N_counts(s) = sum(dists <= window_radius) - 1;
    end

    valid = ~isnan(N_counts);
    if sum(valid) > 1
        density_disorder(i) = var(N_counts(valid));
    else
        density_disorder(i) = NaN;
    end
end

% -------- G6 Structure-Based Disorder --------
structure_disorder = zeros(num_particles, 1);  % 1 - G6

for i = 1:num_particles
    pt = coords(i,:);
    [idxs, ~] = knnsearch(Mdl, pt, 'K', 7);  % self + 6 neighbors
    neighbor_coords = coords(idxs(2:end), :);

    angles = atan2(neighbor_coords(:,2) - pt(2), neighbor_coords(:,1) - pt(1));
    g6 = mean(exp(1i * 6 * angles));
    structure_disorder(i) = 1 - abs(g6);  % invert to get disorder
end

% -------- Combine Metrics --------
alpha = 0.6;
combined_disorder = alpha * density_disorder + (1 - alpha) * structure_disorder;
min_disorder = min(combined_disorder);

% -------- Plot --------
figure('Units', 'normalized', 'Position', [0.2 0.2 0.5 0.6]);
scatter(coords(:,1), coords(:,2), 15, combined_disorder, 'filled');
colormap(parula); 
colorbar;
colorbar_handle = colorbar;
colorbar_handle.Label.String = 'Disorder';
colorbar_handle.Label.FontSize = 18;
set(gca, 'FontSize', 18);
title(sprintf('Ideal Hexagonal Structure Disorder (Mean = %.4f)', min_disorder));
xlabel('X'); ylabel('Y');
axis equal; grid on;


%% Calculate disorder for each frame of the PK-4 system
function disorder_by_frame = compute_combined_disorder(positions_by_frame, alpha, dx, min_disorder, max_disorder)

num_frames = length(positions_by_frame);
disorder_by_frame = cell(num_frames, 1);

window_radius = 30 * dx;         % Physical scale for sampling
samples_per_particle = 60;       % For local density estimation
max_density_clip = 1e4;          % Upper bound to detect outliers

parfor f = 1:num_frames
    frame_data = positions_by_frame{f};
    if isempty(frame_data)
        continue;
    end

    % Convert pixel positions to µm
    x = frame_data(:,2) * dx;
    y = frame_data(:,3) * dx;
    coords = [x, y];
    N = size(coords, 1);

    if N < 10
        warning('Skipping frame %d (too few particles)', f);
        continue;
    end

    frame_width = max(x) - min(x);
    frame_height = max(y) - min(y);
    frame_size = [frame_width, frame_height];

    Mdl = KDTreeSearcher(coords);

    density_disorder = zeros(N, 1);
    structure_disorder = zeros(N, 1);

    for i = 1:N
        pt = coords(i,:);

        % -- Density Disorder --
        N_counts = zeros(samples_per_particle, 1);
        for s = 1:samples_per_particle
            angle = 2 * pi * rand;
            r = window_radius * sqrt(rand);
            cx = pt(1) + r * cos(angle);
            cy = pt(2) + r * sin(angle);

            if cx < 0 || cx > frame_size(1) || cy < 0 || cy > frame_size(2)
                N_counts(s) = NaN;
                continue;
            end

            dists = sqrt((coords(:,1) - cx).^2 + (coords(:,2) - cy).^2);
            N_counts(s) = sum(dists <= window_radius) - 1;
        end

        valid = ~isnan(N_counts);
        if sum(valid) > 1
            var_val = var(N_counts(valid));
            if var_val > max_density_clip
                fprintf('Frame %d, particle %d: High density disorder = %.2f\n', f, i, var_val);
                density_disorder(i) = NaN;
            else
                density_disorder(i) = var_val;
            end
        else
            density_disorder(i) = NaN;
        end

        % -- Structure Disorder (1 - |G6|) --
        [idxs, ~] = knnsearch(Mdl, pt, 'K', 7);  % 6 neighbors + self
        neighbor_coords = coords(idxs(2:end), :);

        if size(neighbor_coords, 1) < 6
            structure_disorder(i) = NaN;
        else
            angles = atan2(neighbor_coords(:,2) - pt(2), neighbor_coords(:,1) - pt(1));
            g6 = mean(exp(1i * 6 * angles));
            structure_disorder(i) = 1 - abs(g6);
        end
    end

    % Combine and clean
    combined_disorder = alpha * density_disorder + (1 - alpha) * structure_disorder;
    combined_disorder(~isfinite(combined_disorder)) = NaN;  % Remove Inf or NaN

    % Optional clamp for rare spikes
    combined_disorder(combined_disorder > 1e4) = NaN;

    % Normalize using global bounds
    normalized_disorder = (combined_disorder - min_disorder) / (max_disorder - min_disorder);
    normalized_disorder = max(0, min(1, normalized_disorder));  % Clamp to [0, 1]

    disorder_by_frame{f} = [x, y, normalized_disorder];
end
end


disorder_by_frame = compute_combined_disorder(positions_by_frame, alpha, dx, min_disorder, avg_combined_disorder);

%%
% Parameters
threshold = 0.5;  % Define blow-up disorder threshold
disorder_by_frame_fixed = disorder_by_frame;  % Initialize copy
bad_points = [];  % [frame_idx, point_idx, bad_value]

for f = 1:length(disorder_by_frame)
    frame_data = disorder_by_frame{f};
    if isempty(frame_data)
        continue;
    end

    x = frame_data(:,1);
    y = frame_data(:,2);
    d = frame_data(:,3);

    % Find outliers above threshold
    bad_idx = find(d > threshold);

    if ~isempty(bad_idx)
        % Find second-highest non-bad value
        valid_d = d(d <= threshold);
        if isempty(valid_d)
            replacement_value = threshold;  % fallback if frame is all bad
        else
            sorted = sort(valid_d, 'descend');
            if length(sorted) >= 2
                replacement_value = 1.25 * sorted(2);  % 25% more than 2nd max
            else
                replacement_value = 1.25 * sorted(1);  % only 1 value
            end
        end

        % Replace and log
        for i = 1:length(bad_idx)
            idx = bad_idx(i);
            bad_points = [bad_points; f, idx, d(idx)];
            d(idx) = replacement_value;
        end
    end

    % Store corrected data
    disorder_by_frame_fixed{f} = [x, y, d];
end

% Position & disorder by frame

% -------- SETTINGS --------
animate = false;        % Set to true to animate
frame_to_show = 11;    % If animate is false, choose the frame here
pause_time = 0.8;       % Pause between frames if animating

% -------- LOOP OVER FRAMES --------
num_frames = length(disorder_by_frame_fixed);

if animate
    figure('Units', 'normalized', 'Position', [0.1, 0.2, 0.9, 0.6]);

    for f = 1:num_frames
        if isempty(disorder_by_frame_fixed{f})
            continue;
        end

        clf;
        frame_data = disorder_by_frame_fixed{f};
        x = frame_data(:,1);
        y = frame_data(:,2);
        d = frame_data(:,3);

        scatter(x, y, 30, d, 'filled');  
        colormap("turbo"); 
        %caxis([0 max(d)]);  % scale to current frame
        colorbar;
        colorbar_handle = colorbar;
        axis equal;
        title(sprintf('Disorder — Frame %d', f), 'FontSize', 24);
        xlabel('X (µm)', 'FontSize', 18); 
        ylabel('Y (µm)', 'FontSize', 18); 
        grid on;
        colorbar_handle.Label.String = 'Disorder';
        colorbar_handle.Label.FontSize = 18;
        set(gca, 'FontSize', 18);
        drawnow;
        pause(pause_time);
    end


else
    f = frame_to_show;
    frame_data = disorder_by_frame_fixed{f};
    x = frame_data(:,1);
    y = frame_data(:,2);
    d = frame_data(:,3);

    figure('Units', 'normalized', 'Position', [0.1, 0.2, 0.9, 0.6]);
    
    max_y = max(y);
    y_lim = [0, 1.1 * max_y];
    
    
    scatter(x, y, 30, d, 'filled');  
    colormap("turbo"); 
    caxis([0 max(disorder_by_frame_fixed{f}(:,3))]);
    colorbar;
    colorbar_handle = colorbar;
    axis equal;
    title(sprintf('Disorder — Frame %d', frame_to_show), 'FontSize', 24);
    xlabel('X (µm)', 'FontSize', 18); 
    ylabel('Y (µm)', 'FontSize', 18);
    set(gca, 'FontSize', 18);
    colorbar_handle.Label.String = 'Disorder';
    colorbar_handle.Label.FontSize = 18;
    grid on;
end


%% Analysis over all frames

% -------- INITIALIZE --------
num_frames = length(disorder_by_frame_fixed);
min_vals = zeros(num_frames,1);
max_vals = zeros(num_frames,1);
mean_vals = zeros(num_frames,1);

% -------- LOOP OVER FRAMES --------
for f = 1:num_frames
    if isempty(disorder_by_frame_fixed{f})
        min_vals(f) = NaN;
        max_vals(f) = NaN;
        mean_vals(f) = NaN;
        continue;
    end

    d = disorder_by_frame_fixed{f}(:,3);
    min_vals(f) = min(d);
    max_vals(f) = max(d);
    mean_vals(f) = mean(d);
end

% -------- Compute Overall Stats --------
mean_of_means = mean(mean_vals, 'omitnan');
mean_of_max = mean(max_vals, 'omitnan');

fprintf('Mean of Frame Means: %.4f\n', mean_of_means);
fprintf('Mean of Frame Maxima: %.4f\n', mean_of_max);

% -------- PLOT --------
figure('Units', 'normalized', 'Position', [0.1 0.2 0.8 0.5]);
plot(min_vals, 'k.-', 'DisplayName', 'Min');
hold on;
plot(max_vals, 'r.-', 'DisplayName', 'Max');
plot(mean_vals, 'b.-', 'DisplayName', 'Mean');

xlabel('Frame', 'FontSize', 18);
ylabel('Disorder Value', 'FontSize', 18);
title('Disorder Stats Over Time', 'FontSize', 24);
legend('show', 'FontSize', 18);
set(gca, 'FontSize', 18);
grid on;


close all