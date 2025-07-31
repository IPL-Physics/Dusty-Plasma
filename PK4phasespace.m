% Extract columns
data = C1705Pa035mAFov14x2mm;

N0_PARTICLES = max(data{end,1}); % Total number of trajectories detected with ImageJ
N_TIME_STEPS = max(data{end,2}); % Total number of frames analyzed with ImageJ

% Note that your columns may be named X and Y, instead of x and y. Modify the next line accordingly.
alltraj = [data{:,1} data{:,2} data{:,3} data{:,4}]; % Create a matrix containing the relevant vecor columns.

% Data from image analysis is often 2D.
% However, if simulation data is imported, there will be a third coordinate as well and the dimension will be 3D.
N_DIM = 2;

% The timestep is equal to 1/frame rate. For example, for a frame rate 71.4 fr/s, dT = 1/71.4 = 0.014.
% Note that if dT includes more digits, the program gives an error!
framerate = 71.4;
dT = 1/framerate; % Timestep in units of [s]

% Approximate pixel resolution in micrometers from the PK-4 characterization paper: Pustylnik, M. Y., et al. "Plasmakristall-4: New complex (dusty) plasma laboratory on board the International Space Station." Review of Scientific Instruments 87.9 (2016): 093505.
% From table IV of this paper, we see that te pixel sizes vary slightly in
% the x and y directions and are slightly different for camera 1 and camera 2. Here we assume an average pixel size.
dx = 14.20; % Pixel size in units of micrometers [um]

% Cut out short tracks
tracks = cell(N0_PARTICLES, 1); % This is the right format for data input in @msdanalyzer.
short = zeros(N0_PARTICLES, 1); % Array that will count the # of trajectories that are shorter than 10 timesteps
% Note that since we store all tracks except the ones shorter than 10 timesteps, the number of nonempty matrices
% in the cell "tracks" is equal to the total number of particles detected minus the number of particles stored in "short10".
frame_cutoff = 10; % tracks that only have this many or fewer frames will be removed

for i = 1:N0_PARTICLES
    A = alltraj(:,1)==i; % Extract logical data for the ith trajectory.
    B = A.*alltraj; % Extract the actual data for the ith trajectory.
    B(~any(B,2),:)=[]; % Remove all rolls that consist of zeros.
    B = B(:,2:4); % Extract the time frames and the x- and y- positions.
    if size(B,1) < frame_cutoff % Ignoring tracks shorter than frame cutoff timesteps.
        short(i,:) = 1; % Basically, nothing happens in this iteration and the loop will restart with i=i+1
    else
        D = B(:,1); % Extract the time steps for the ith trajectory (as recorded from ImageJ)
        C = B(:,2:3); % Extract x and y positions
        R = dx.*C; % Convert units from pixels to microns
        T = B(:,1);
        time = T * dT; % Create a time sequence corresponding to the image data frame rate in [s].
        tracks{i} = [time R];
    end
end
tracks = tracks(~any(cellfun('isempty', tracks), 2), :);
short(~any(short,2),:)=[]; % Remove all rows that consist of zeros.

N_tracks = length(tracks); % Get the number of tracks >= 10 timesteps
N_shorts = length(short); % Get the number of tracks <10 timesteps

%%
% Initialize arrays to store velocity components
vx = cell(size(tracks));  % Cell array for x-velocity
vy = cell(size(tracks));  % Cell array for y-velocity

% Iterate through each trajectory in tracks
for i = 1:length(tracks)
    % Extract the current trajectory (timesteps, x, and y positions)
    current_track = tracks{i};  % Current trajectory as a double array
    timesteps = current_track(:, 1);  % First column: timesteps
    x_positions = current_track(:, 2);  % Second column: x positions
    y_positions = current_track(:, 3);  % Third column: y positions
    
    % Compute velocities for x and y using differences in positions and timesteps
    vx{i} = [0; diff(x_positions) ./ diff(timesteps)];  % Set initial velocity to 0
    vy{i} = [0; diff(y_positions) ./ diff(timesteps)];  % Set initial velocity to 0
end

% Add vx and vy to tracks
for i = 1:length(tracks)
    tracks{i}(:, 4) = vx{i};  % Add vx as a new column
    tracks{i}(:, 5) = vy{i};  % Add vy as another new column
end
%% Animation of phase space
figure;

% Determine global max and min for fixed axis limits
all_x = cell2mat(cellfun(@(c) c(:, 2), tracks, 'UniformOutput', false));
all_y = cell2mat(cellfun(@(c) c(:, 3), tracks, 'UniformOutput', false));
all_vx = cell2mat(cellfun(@(c) c(:, 4), tracks, 'UniformOutput', false));
all_vy = cell2mat(cellfun(@(c) c(:, 5), tracks, 'UniformOutput', false));

x_limits = [min(all_x), max(all_x)];
vx_limits = [min(all_vx), max(all_vx)];
y_limits = [min(all_y), max(all_y)];
vy_limits = [min(all_vy), max(all_vy)];

% Subplot 1: vx vs x
subplot(1, 2, 1);
xlabel('x');
ylabel('v_x');
title('');
grid on;
xlim(x_limits);
ylim(vx_limits);

% Subplot 2: vy vs y
subplot(1, 2, 2);
xlabel('y');
ylabel('v_y');
title('');
grid on;
xlim(y_limits);
ylim(vy_limits);

% Determine the global time range across all trajectories
all_timesteps = cell2mat(cellfun(@(c) c(:, 1), tracks, 'UniformOutput', false));
min_time = min(all_timesteps);
max_time = max(all_timesteps);
time_steps = min_time:dT:max_time;  % Generate time steps based on dT

% Animation loop for phase space
for t = time_steps
    % Clear plots for current timestep
    subplot(1, 2, 1);
    cla;
    hold on;
    subplot(1, 2, 2);
    cla;
    hold on;
    
    % Initialize arrays for current timestep
    x_all = [];
    vx_all = [];
    y_all = [];
    vy_all = [];
    
    % Loop through each trajectory and collect points for current timestep
    for i = 1:length(tracks)
        current_track = tracks{i};
        % Find the index closest to the current time t
        [~, idx] = min(abs(current_track(:, 1) - t));
        
        if abs(current_track(idx, 1) - t) < dT/2  % Include only if within half dT
            x_all = [x_all; current_track(idx, 2)];
            vx_all = [vx_all; current_track(idx, 4)];
            y_all = [y_all; current_track(idx, 3)];
            vy_all = [vy_all; current_track(idx, 5)];
        end
    end
    
    % Plot phase space points for the current timestep
    subplot(1, 2, 1);
    plot(x_all, vx_all, 'b.', 'MarkerSize', 8);  % vx vs x in blue
    
    subplot(1, 2, 2);
    plot(y_all, vy_all, 'r.', 'MarkerSize', 8);  % vy vs y in red
    
    % Update titles with current time
    sgtitle(sprintf('Phase Space at Time: %.3f', t));
    
    drawnow;  % Update the plots
    pause(0.014);  % Control the animation speed
end

%%
% Select a single particle (e.g., the first particle in tracks)
particle_id = 50;  % Change this to select a different particle
selected_track = tracks{particle_id};

% Extract data for the selected particle
time = selected_track(:, 1);  % Time data
x = selected_track(:, 2);     % Position in x
y = selected_track(:, 3);     % Position in y
vx = selected_track(:, 4);    % Velocity in x
vy = selected_track(:, 5);    % Velocity in y

% Create a figure for the static plot
figure;

% Subplot 1: vx vs x
subplot(1, 2, 1);
hold on;
scatter(x, vx, 30, time, 'filled');  % Scatter plot of phase points, color-coded by time
plot(x, vx, 'b-', 'LineWidth', 1.5);  % Connect points with a line
xlabel('x');
ylabel('v_x');
title('Phase Space: v_x vs x');
grid on;
cb1 = colorbar;  % Add colorbar
colormap(jet);  % Use colormap for time
clim([min(time), max(time)]);  % Ensure color scale reflects time
cb1.Label.String = 'Timesteps (s)';  % Label the colorbar
xlim([min(x), max(x)]);
ylim([min(vx), max(vx)]);

% Subplot 2: vy vs y
subplot(1, 2, 2);
hold on;
scatter(y, vy, 30, time, 'filled');  % Scatter plot of phase points, color-coded by time
plot(y, vy, 'r-', 'LineWidth', 1.5);  % Connect points with a line
xlabel('y');
ylabel('v_y');
title('Phase Space: v_y vs y');
grid on;
cb2 = colorbar;  % Add colorbar
colormap(jet);  % Use the same colormap for consistency
clim([min(time), max(time)]);  % Ensure color scale reflects time
cb2.Label.String = 'Timesteps (s)';  % Label the colorbar
xlim([min(y), max(y)]);
ylim([min(vy), max(vy)]);

% Add an overall title
sgtitle(sprintf('Phase Space Trajectory for Particle %d', particle_id));

%% Contour smoothed and average over all trajectories and timesteps
% Aggregate vx-x and vy-y data across all particles and timesteps
all_x = cell2mat(cellfun(@(c) c(:, 2), tracks, 'UniformOutput', false));  % x positions
all_vx = cell2mat(cellfun(@(c) c(:, 4), tracks, 'UniformOutput', false)); % vx velocities
all_y = cell2mat(cellfun(@(c) c(:, 3), tracks, 'UniformOutput', false));  % y positions
all_vy = cell2mat(cellfun(@(c) c(:, 5), tracks, 'UniformOutput', false)); % vy velocities

% Define bin edges for phase space dimensions
x_edges = linspace(min(all_x), max(all_x), 100);     % 100 bins for x
vx_edges = linspace(min(all_vx), max(all_vx), 100); % 100 bins for vx
y_edges = linspace(min(all_y), max(all_y), 100);     % 100 bins for y
vy_edges = linspace(min(all_vy), max(all_vy), 100); % 100 bins for vy

% Create 2D histograms for phase space density
[N_vx_x, x_bins, vx_bins] = histcounts2(all_x, all_vx, x_edges, vx_edges);  % vx vs x
[N_vy_y, y_bins, vy_bins] = histcounts2(all_y, all_vy, y_edges, vy_edges);  % vy vs y

% Normalize histograms to create density distributions
density_vx_x = N_vx_x / sum(N_vx_x(:));  % Normalize vx vs x
density_vy_y = N_vy_y / sum(N_vy_y(:));  % Normalize vy vs y

% Smooth the density using Gaussian filters
density_vx_x_smooth = imgaussfilt(density_vx_x, 2);  % Smooth vx vs x
density_vy_y_smooth = imgaussfilt(density_vy_y, 2);  % Smooth vy vs y

% Plot phase contour for vx vs x
figure;

subplot(1, 2, 1);
contourf(x_bins(1:end-1), vx_bins(1:end-1), density_vx_x_smooth', 20, 'LineColor', 'none');  % 20 contour levels
cb1 = colorbar;
colormap(jet);
cb1.Label.String = 'Particle-Time Density';  % Label the colorbar
xlabel('x');
ylabel('v_x');
title('Phase Contour: v_x vs x');

% Plot phase contour for vy vs y
subplot(1, 2, 2);
contourf(y_bins(1:end-1), vy_bins(1:end-1), density_vy_y_smooth', 20, 'LineColor', 'none');  % 20 contour levels
cb2 = colorbar;
colormap(jet);
cb2.Label.String = 'Particle-Time Density';  % Label the colorbar
xlabel('y');
ylabel('v_y');
title('Phase Contour: v_y vs y');

%% Smooth contour for all trajectories animated over timesteps
% Determine the unique timesteps across all particles
all_time = cell2mat(cellfun(@(c) c(:, 1), tracks, 'UniformOutput', false));
unique_timesteps = unique(all_time);

% Create figure
figure;

% Iterate through each timestep to generate the animation
for t_idx = 1:length(unique_timesteps)
    t = unique_timesteps(t_idx);
    
    % Initialize arrays for the current timestep
    x_t = [];
    vx_t = [];
    y_t = [];
    vy_t = [];
    
    % Extract data for all particles at the current timestep
    for i = 1:length(tracks)
        track = tracks{i};
        % Find the closest data point to the current timestep
        [~, idx] = min(abs(track(:, 1) - t));
        if abs(track(idx, 1) - t) < 1e-5  % Include only if close to the current timestep
            x_t = [x_t; track(idx, 2)];
            vx_t = [vx_t; track(idx, 4)];
            y_t = [y_t; track(idx, 3)];
            vy_t = [vy_t; track(idx, 5)];
        end
    end
    
    % Define bin edges for histograms
    x_edges = linspace(min(all_x), max(all_x), 100);
    vx_edges = linspace(min(all_vx), max(all_vx), 100);
    y_edges = linspace(min(all_y), max(all_y), 100);
    vy_edges = linspace(min(all_vy), max(all_vy), 100);
    
    % Create 2D histograms for the current timestep
    [N_vx_x, x_bins, vx_bins] = histcounts2(x_t, vx_t, x_edges, vx_edges);
    [N_vy_y, y_bins, vy_bins] = histcounts2(y_t, vy_t, y_edges, vy_edges);
    
    % Normalize the histograms to create density distributions
    density_vx_x = N_vx_x / sum(N_vx_x(:));
    density_vy_y = N_vy_y / sum(N_vy_y(:));
    
    % Smooth the density using a Gaussian filter
    density_vx_x_smooth = imgaussfilt(density_vx_x, 2);
    density_vy_y_smooth = imgaussfilt(density_vy_y, 2);
    
    % Plot phase contour for vx vs x
    subplot(1, 2, 1);
    contourf(x_bins(1:end-1), vx_bins(1:end-1), density_vx_x_smooth', 20, 'LineColor', 'none');
    cb1 = colorbar;
    colormap(jet);
    cb1.Label.String = 'Particle Density';  % Label the colorbar
    xlabel('x');
    ylabel('v_x');
    title(sprintf('Phase Contour: v_x vs x at t = %.2f', t));
    
    % Plot phase contour for vy vs y
    subplot(1, 2, 2);
    contourf(y_bins(1:end-1), vy_bins(1:end-1), density_vy_y_smooth', 20, 'LineColor', 'none');
    cb2 = colorbar;
    colormap(jet);
    cb2.Label.String = 'Particle Density';  % Label the colorbar
    xlabel('y');
    ylabel('v_y');
    title(sprintf('Phase Contour: v_y vs y at t = %.2f', t));
    
    % Update the plots
    sgtitle('Animated Phase Contour Plot');
    drawnow;
    
    % Pause for animation speed
    pause(0.1);
end

%% Smooth contour of timesteps animated for all trajectories
% Create figure
figure;

% Loop through each trajectory to generate the animation
for particle_idx = 1:length(tracks)
    % Extract data for the current trajectory
    track = tracks{particle_idx};
    x_t = track(:, 2);     % x positions
    vx_t = track(:, 4);    % vx velocities
    y_t = track(:, 3);     % y positions
    vy_t = track(:, 5);    % vy velocities
    
    % Define bin edges for histograms
    x_edges = linspace(min(x_t), max(x_t), 100);
    vx_edges = linspace(min(vx_t), max(vx_t), 100);
    y_edges = linspace(min(y_t), max(y_t), 100);
    vy_edges = linspace(min(vy_t), max(vy_t), 100);
    
    % Create 2D histograms for phase space density
    [N_vx_x, x_bins, vx_bins] = histcounts2(x_t, vx_t, x_edges, vx_edges);
    [N_vy_y, y_bins, vy_bins] = histcounts2(y_t, vy_t, y_edges, vy_edges);
    
    % Normalize the histograms to create density distributions
    density_vx_x = N_vx_x / sum(N_vx_x(:));
    density_vy_y = N_vy_y / sum(N_vy_y(:));
    
    % Smooth the density using a Gaussian filter
    density_vx_x_smooth = imgaussfilt(density_vx_x, 2);
    density_vy_y_smooth = imgaussfilt(density_vy_y, 2);
    
    % Plot phase contour for vx vs x
    subplot(1, 2, 1);
    contourf(x_bins(1:end-1), vx_bins(1:end-1), density_vx_x_smooth', 20, 'LineColor', 'none');
    colorbar;
    colormap(jet);
    xlabel('x');
    ylabel('v_x');
    title(sprintf('Phase Contour: v_x vs x for Particle %d', particle_idx));
    
    % Plot phase contour for vy vs y
    subplot(1, 2, 2);
    contourf(y_bins(1:end-1), vy_bins(1:end-1), density_vy_y_smooth', 20, 'LineColor', 'none');
    colorbar;
    colormap(jet);
    xlabel('y');
    ylabel('v_y');
    title(sprintf('Phase Contour: v_y vs y for Particle %d', particle_idx));
    
    % Update the plots
    sgtitle('Animated Phase Contour for Individual Trajectories');
    drawnow;
    
    % Pause for animation speed
    pause(0.5);
end

%% Entropy

% Define bins for coarse-graining
x_edges = linspace(min(all_x), max(all_x), 30);
vx_edges = linspace(min(all_vx), max(all_vx), 30);
entropy_over_time = [];

% Get all unique time steps
all_times = unique(cell2mat(cellfun(@(c) c(:,1), tracks, 'UniformOutput', false)));

for t = all_times'
    x_vals = [];
    vx_vals = [];
    for i = 1:length(tracks)
        [~, idx] = min(abs(tracks{i}(:,1) - t));
        if abs(tracks{i}(idx,1) - t) < dT/2
            x_vals = [x_vals; tracks{i}(idx,2)];
            vx_vals = [vx_vals; tracks{i}(idx,4)];
        end
    end
    if isempty(x_vals)
        continue;
    end
    [N, ~, ~] = histcounts2(x_vals, vx_vals, x_edges, vx_edges);
    P = N / sum(N(:));
    P = P(P > 0);
    S = -sum(P .* log(P));
    entropy_over_time(end+1,:) = [t, S];
end

% Plot entropy vs time
figure;
plot(entropy_over_time(:,1), entropy_over_time(:,2), 'k-o');
xlabel('Time [s]');
ylabel('Phase Space Entropy');
title('Entropy Evolution in Phase Space');
grid on;

%%

% Combine all [x, vx] points
points = cell2mat(cellfun(@(c) c(:,[2 4]), tracks, 'UniformOutput', false));
epsilons = [50, 40, 30, 20, 10];  % Smaller = finer resolution
N_eps = zeros(size(epsilons));

x_range = [min(points(:,1)), max(points(:,1))];
vx_range = [min(points(:,2)), max(points(:,2))];

for i = 1:length(epsilons)
    eps = epsilons(i);
    x_bins = linspace(x_range(1), x_range(2), eps);
    vx_bins = linspace(vx_range(1), vx_range(2), eps);
    [counts,~,~] = histcounts2(points(:,1), points(:,2), x_bins, vx_bins);
    N_eps(i) = sum(counts(:) > 0);
end

% Estimate slope in log-log
log_eps = log(1 ./ epsilons);
log_N = log(N_eps);
p = polyfit(log_eps, log_N, 1);
D = p(1);

fprintf('Estimated Fractal Dimension: %.3f\n', D);

% Plot
figure;
plot(log_eps, log_N, 'o-');
xlabel('log(1/ε)');
ylabel('log(N(ε))');
title('Box-Counting for Fractal Dimension');
grid on;


%%

x_edges = linspace(min(all_x), max(all_x), 20);
vx_edges = linspace(min(all_vx), max(all_vx), 20);

[X_centers, VX_centers] = meshgrid( ...
    0.5*(x_edges(1:end-1) + x_edges(2:end)), ...
    0.5*(vx_edges(1:end-1) + vx_edges(2:end)) );

Jx = zeros(size(X_centers));
Jvx = zeros(size(X_centers));
counts = zeros(size(X_centers));

for i = 1:length(tracks)
    traj = tracks{i};
    for j = 2:size(traj,1)-1
        x = traj(j,2);
        vx = traj(j,4);
        ax = (traj(j+1,4) - traj(j-1,4)) / (traj(j+1,1) - traj(j-1,1)); % central difference

        x_idx = find(histcounts(x, x_edges));
        vx_idx = find(histcounts(vx, vx_edges));
        if ~isempty(x_idx) && ~isempty(vx_idx)
            Jx(vx_idx, x_idx) = Jx(vx_idx, x_idx) + vx;
            Jvx(vx_idx, x_idx) = Jvx(vx_idx, x_idx) + ax;
            counts(vx_idx, x_idx) = counts(vx_idx, x_idx) + 1;
        end
    end
end

% Normalize
Jx = Jx ./ max(counts,1);
Jvx = Jvx ./ max(counts,1);

% Plot current field
figure;
quiver(X_centers, VX_centers, Jx, Jvx, 'b');
xlabel('x [μm]');
ylabel('v_x [μm/s]');
title('Phase Space Current (x vs v_x)');
axis tight;
grid on;



%% Lyaponov Exponent 
% Initialize arrays to store Lyapunov exponents for x and y
lyapunov_exponents_x = zeros(length(tracks), 1);
lyapunov_exponents_y = zeros(length(tracks), 1);

% Loop through each trajectory to compute Lyapunov exponents
for i = 1:length(tracks)
    current_track = tracks{i};
    timesteps = current_track(:, 1);  % Time column
    x_positions = current_track(:, 2);
    y_positions = current_track(:, 3);
    vx = current_track(:, 4);
    vy = current_track(:, 5);
    
    % Compute displacement in phase space for x and y
    dx = diff(x_positions);
    dvx = diff(vx);
    dy = diff(y_positions);
    dvy = diff(vy);
    
    % Phase space displacement (x, vx) and (y, vy)
    displacement_x = sqrt(dx.^2 + dvx.^2);
    displacement_y = sqrt(dy.^2 + dvy.^2);
    
    % Logarithm of displacement (remove zeros to avoid -Inf)
    log_disp_x = log(displacement_x(displacement_x > 0));
    log_disp_y = log(displacement_y(displacement_y > 0));
    
    % Time differences corresponding to displacements
    time_differences = diff(timesteps);
    valid_times = timesteps(2:end);  % Remove the first time point
    
    % Fit linear model to log displacement for x phase space
    if ~isempty(log_disp_x)
        fit_x = polyfit(valid_times(1:length(log_disp_x)), log_disp_x, 1);
        lyapunov_exponents_x(i) = fit_x(1);  % Slope gives Lyapunov exponent for x
    end
    
    % Fit linear model to log displacement for y phase space
    if ~isempty(log_disp_y)
        fit_y = polyfit(valid_times(1:length(log_disp_y)), log_disp_y, 1);
        lyapunov_exponents_y(i) = fit_y(1);  % Slope gives Lyapunov exponent for y
    end
end

% Plot Lyapunov exponents
figure;
scatter(1:length(tracks), lyapunov_exponents_x, 'b.', 'DisplayName', '\lambda_x');  % Blue for x
hold on;
scatter(1:length(tracks), lyapunov_exponents_y, 'r.', 'DisplayName', '\lambda_y');  % Red for y
xlabel('Trajectory Index');
ylabel('Lyapunov Exponent');
title('Lyapunov Exponents for All Trajectories');
legend('Location', 'Best');
grid on;

%%
% Wavelet Analysis
% Initialize arrays to store all displacements
all_x_displacements = [];
all_y_displacements = [];

% Loop through each particle to compute displacements
for i = 1:length(tracks)
    % Extract data for each particle
    data = tracks{i};
    x_position = data(:, 2);
    y_position = data(:, 3);
    
    % Compute displacements for x and y
    x_displacement = diff(x_position);  % Δx
    y_displacement = diff(y_position);  % Δy
    
    % Append to the aggregate arrays
    all_x_displacements = [all_x_displacements; x_displacement];
    all_y_displacements = [all_y_displacements; y_displacement];
end

% Number of decomposition levels
nLevels = 5;

% Perform Discrete Wavelet Transform (Haar)
[coeff_x, l_x] = wavedec(all_x_displacements, nLevels, 'haar');
[coeff_y, l_y] = wavedec(all_y_displacements, nLevels, 'haar');

% Compute energy at each scale for x-displacement
energy_x = zeros(1, nLevels);
for j = 1:nLevels
    details_x = detcoef(coeff_x, l_x, j); % Detail coefficients for scale j
    energy_x(j) = sum(details_x .^ 2);   % Energy (sum of squared coefficients)
end

% Compute energy at each scale for y-displacement
energy_y = zeros(1, nLevels);
for j = 1:nLevels
    details_y = detcoef(coeff_y, l_y, j); % Detail coefficients for scale j
    energy_y(j) = sum(details_y .^ 2);   % Energy (sum of squared coefficients)
end

% Plot energy at each scale
figure;
bar(1:nLevels, [energy_x; energy_y]', 'grouped');
xlabel('Scale (Level)');
ylabel('Energy (Sum of Squared Coefficients)');
legend({'x-displacement', 'y-displacement'});
title('Wavelet Energy at Each Scale for Aggregated Displacements');

%%
% Parameters for temporal analysis
window_size = 10000; % Number of displacement points per window (adjust as needed)
num_windows = floor(length(all_x_displacements) / window_size);

% Initialize arrays to store energy over time
energy_x_time = zeros(num_windows, nLevels);
energy_y_time = zeros(num_windows, nLevels);

% Loop through each time window
for w = 1:num_windows
    % Extract data for current window
    start_idx = (w - 1) * window_size + 1;
    end_idx = w * window_size;
    x_window = all_x_displacements(start_idx:end_idx);
    y_window = all_y_displacements(start_idx:end_idx);
    
    % Perform Discrete Wavelet Transform
    [coeff_x, l_x] = wavedec(x_window, nLevels, 'haar');
    [coeff_y, l_y] = wavedec(y_window, nLevels, 'haar');
    
    % Compute energy at each scale
    for j = 1:nLevels
        details_x = detcoef(coeff_x, l_x, j);
        details_y = detcoef(coeff_y, l_y, j);
        energy_x_time(w, j) = sum(details_x .^ 2);
        energy_y_time(w, j) = sum(details_y .^ 2);
    end
end

% Plot energy evolution over time for x-displacement
figure;
imagesc(1:num_windows, 1:nLevels, energy_x_time');
colorbar;
xlabel('Time Window');
ylabel('Scale (Level)');
title('Wavelet Energy Evolution (x-displacement)');
set(gca, 'YDir', 'normal'); % Ensure scale is displayed top-down

% Plot energy evolution over time for y-displacement
figure;
imagesc(1:num_windows, 1:nLevels, energy_y_time');
colorbar;
xlabel('Time Window');
ylabel('Scale (Level)');
title('Wavelet Energy Evolution (y-displacement)');
set(gca, 'YDir', 'normal');
%%
% Parameters
delta_t = 0.014; % Time step between displacements in seconds
window_size = 5000; % Number of data points per time window (adjust as needed)
T_window = window_size * delta_t; % Time per window in seconds

% Total number of windows
num_windows = floor(length(all_x_displacements) / window_size);

% Define time edges using actual time values
time_edges = linspace(0, num_windows * T_window, num_windows + 1);

% Compute correlations for each time window
correlation_time = zeros(1, num_windows);
for w = 1:num_windows
    % Extract data for the current window
    start_idx = (w - 1) * window_size + 1;
    end_idx = w * window_size;
    x_window = all_x_displacements(start_idx:end_idx);
    y_window = all_y_displacements(start_idx:end_idx);
    
    % Compute Pearson correlation
    correlation_time(w) = corr(x_window, y_window);
end

% Compute time midpoints for each time window
time_midpoints = (time_edges(1:end-1) + time_edges(2:end)) / 2;

% Plot correlation over time
figure;
plot(time_midpoints, correlation_time, '-o');
xlabel('Time (seconds)');
ylabel('Correlation Coefficient');
title('Correlation between x and y Displacements Over Time');
grid on;

% Wavelet energy plots
figure;
imagesc(time_midpoints, 1:nLevels, energy_x_time');
colorbar;
xlabel('Time (seconds)');
ylabel('Scale (Level)');
title('Wavelet Energy Evolution (x-displacement)');
set(gca, 'YDir', 'normal');

figure;
imagesc(time_midpoints, 1:nLevels, energy_y_time');
colorbar;
xlabel('Time (seconds)');
ylabel('Scale (Level)');
title('Wavelet Energy Evolution (y-displacement)');
set(gca, 'YDir', 'normal');
