% Step 1: Open the figure
figPath = '/Users/bradleyandrew/Desktop/Research/PK-4 Experiment/Newest PK4 Data analysis and Ideas/track plots/Camera 1/';
openfig(figPath, 'invisible');

% Step 2: Extract trajectory lines
lines = findall(gcf, 'Type', 'Line');

% Step 3: Initialize containers
all_steps_cell = cell(length(lines), 1);
all_steps_flat = [];  % Will collect all step sizes for overall stats

for i = 1:length(lines)
    x = get(lines(i), 'XData')';
    y = get(lines(i), 'YData')';
    
    if length(x) < 10
        continue;
    end
    
    steps = sqrt(diff(x).^2 + diff(y).^2);  % Step sizes for this trajectory
    all_steps_cell{i} = steps;              % Store in cell
    all_steps_flat = [all_steps_flat; steps];  % Append to global list
end

% Step 4: Statistics from flattened list
min_step = min(all_steps_flat);
mean_step = mean(all_steps_flat);
max_step = max(all_steps_flat);

% Step 5: Detect large jumps (≥ 10x the mean)
threshold = 10 * mean_step;
num_large_jumps = sum(all_steps_flat >= threshold);
percent_large_jumps = 100 * num_large_jumps / length(all_steps_flat);

% Step 6: Print results
fprintf('--- Step Size Statistics (µm) ---\n');
fprintf('Minimum step size: %.4f µm\n', min_step);
fprintf('Average step size: %.4f µm\n', mean_step);
fprintf('Maximum step size: %.4f µm\n', max_step);

fprintf('\n--- Large Jump Detection ---\n');
fprintf('Threshold for large jumps: %.2f µm (10× avg)\n', threshold);
fprintf('Number of large jumps: %d\n', num_large_jumps);
fprintf('Percentage of large jumps: %.2f%%\n', percent_large_jumps);


% Step 7: Check for clustering of large jumps
clustered_large_jumps = 0;

% Loop through each trajectory
for i = 1:length(all_steps_cell)
    steps = all_steps_cell{i};
    large_jump_indices = find(steps >= threshold);
    
    for j = 1:length(large_jump_indices)
        idx = large_jump_indices(j);

        % Look ±5 indices around current index (excluding self)
        for k = max(1, idx-5):min(length(steps), idx+5)
            if k ~= idx && steps(k) >= threshold
                clustered_large_jumps = clustered_large_jumps + 1;
                break;  % Only count once per large jump
            end
        end
    end
end



% Step 8: Compute percentage
percent_clustered_large_jumps = 100 * clustered_large_jumps / num_large_jumps;

% Step 9: Print clustering result
fprintf('\n--- Large Jump Clustering ---\n');
fprintf('Large jumps with another within 10 steps: %d\n', clustered_large_jumps);
fprintf('Percentage of clustered large jumps: %.2f%%\n', percent_clustered_large_jumps);

% Step 10: Calculate average angle of jumps (in degrees)
angles = [];

for i = 1:length(all_steps_cell)
    x = get(lines(i), 'XData')';
    y = get(lines(i), 'YData')';

    if length(x) < 10
        continue;
    end

    dx = diff(x);
    dy = diff(y);

    % Compute angles in radians, then convert to degrees
    theta = atan2d(dy, dx);  % atan2 handles correct quadrant
    angles = [angles; theta];  % Collect all angles
end

% Mean angle (circular mean can also be considered for wrapping)
mean_angle = mean(angles);

% Print result
fprintf('\n--- Average Jump Angle ---\n');
fprintf('Mean jump direction: %.2f degrees\n', mean_angle);

% Step 11: Plot histogram of jump angles
figure;
histogram(angles, 'BinWidth', 1);  % 10° bins
xlabel('Jump Angle (degrees)');
ylabel('Count');
title('Histogram of Jump Directions');
grid on;

% Collect angles only for large jumps
large_jump_angles = [];

for i = 1:length(all_steps_cell)
    x = get(lines(i), 'XData')';
    y = get(lines(i), 'YData')';
    
    if length(x) < 10
        continue;
    end
    
    dx = diff(x);
    dy = diff(y);
    step_sizes = sqrt(dx.^2 + dy.^2);
    
    % Logical index for large jumps
    is_large = step_sizes >= threshold;
    
    % Get angles for large jumps
    theta_large = atan2d(dy(is_large), dx(is_large));
    large_jump_angles = [large_jump_angles; theta_large];
end

% Compute and print average angle of large jumps
mean_large_angle = mean(large_jump_angles);
fprintf('\n--- Average Angle of Large Jumps ---\n');
fprintf('Mean direction of large jumps: %.2f degrees\n', mean_large_angle);

% Plot histogram for large jump angles
figure;
histogram(large_jump_angles, 'BinWidth', 10);
xlabel('Jump Angle (degrees)');
ylabel('Count');
title('Histogram of Large Jump Directions');
grid on;

% Initialize accumulators
total_dx_all = 0; total_dy_all = 0; count_all = 0;
total_dx_large = 0; total_dy_large = 0; count_large = 0;

for i = 1:length(all_steps_cell)
    x = get(lines(i), 'XData')';
    y = get(lines(i), 'YData')';

    if length(x) < 2
        continue;
    end

    dx = diff(x);
    dy = diff(y);
    step_sizes = sqrt(dx.^2 + dy.^2);
    
    % For all jumps
    total_dx_all = total_dx_all + sum(abs(dx));
    total_dy_all = total_dy_all + sum(abs(dy));
    count_all = count_all + length(dx);

    % For large jumps
    is_large = step_sizes >= threshold;
    total_dx_large = total_dx_large + sum(abs(dx(is_large)));
    total_dy_large = total_dy_large + sum(abs(dy(is_large)));
    count_large = count_large + sum(is_large);
end

% Compute means
mean_dx_all = total_dx_all / count_all;
mean_dy_all = total_dy_all / count_all;
mean_dx_large = total_dx_large / count_large;
mean_dy_large = total_dy_large / count_large;

% Compute relative boost
x_ratio = mean_dx_large / mean_dx_all;
y_ratio = mean_dy_large / mean_dy_all;

fprintf('\n--- Directional Preference Comparison ---\n');
fprintf('Mean |Δx| (all): %.2f µm | large: %.2f µm | x ratio: %.2f\n', mean_dx_all, mean_dx_large, x_ratio);
fprintf('Mean |Δy| (all): %.2f µm | large: %.2f µm | y ratio: %.2f\n', mean_dy_all, mean_dy_large, y_ratio);


%%

% Initialize
angles_all = [];
angles_large = [];

for f = 1:length(files)
    figPath = fullfile(figDir, files(f).name);
    openfig(figPath, 'invisible');
    lines = findall(gcf, 'Type', 'Line');

    all_steps = [];

    for i = 1:length(lines)
        x = get(lines(i), 'XData')';
        y = get(lines(i), 'YData')';

        if length(x) < 2
            continue;
        end

        dx = diff(x);
        dy = diff(y);
        steps = sqrt(dx.^2 + dy.^2);

        angles = atan2d(dy, dx);  % angles in degrees
        angles_all = [angles_all; angles];
        all_steps = [all_steps; steps];
    end

    threshold = 10 * mean(all_steps);

    for i = 1:length(lines)
        x = get(lines(i), 'XData')';
        y = get(lines(i), 'YData')';

        if length(x) < 2
            continue;
        end

        dx = diff(x);
        dy = diff(y);
        steps = sqrt(dx.^2 + dy.^2);
        angles = atan2d(dy, dx);
        is_large = steps >= threshold;

        angles_large = [angles_large; angles(is_large)];
    end

    close(gcf);
end

% Clean up
angles_all = angles_all(~isnan(angles_all));
angles_large = angles_large(~isnan(angles_large));

% T-test: compare angle distributions
[h, p, ci, stats] = ttest2(angles_large, angles_all, 'Vartype', 'unequal');

fprintf('\n--- Angle-Based Directional Preference Test ---\n');
fprintf('Mean angle (all): %.2f degrees\n', mean(angles_all));
fprintf('Mean angle (large): %.2f degrees\n', mean(angles_large));
fprintf('t-statistic: %.4f\n', stats.tstat);
fprintf('p-value: %.4f\n', p);
if h
    fprintf('Result: Large jumps have a significantly different directional angle.\n');
else
    fprintf('Result: No significant difference in angle between large and typical jumps.\n');
end

[p_wilcoxon, h_wilcoxon, stats_wilcoxon] = ranksum(angles_large, angles_all);

fprintf('\n--- Mann–Whitney U Test (Angle) ---\n');
fprintf('p-value: %.4f\n', p_wilcoxon);
if h_wilcoxon
    fprintf('Result: Statistically significant difference in angular distributions.\n');
else
    fprintf('Result: No significant difference detected.\n');
end


%%
% analyze_jump_components

% Set directory with .fig files
figDir = '/Users/bradleyandrew/Desktop/Research/PK-4 Experiment/Newest PK4 Data analysis and Ideas/track plots/Camera 1/';
files = dir(fullfile(figDir, '*1000above10frs.fig'));

% Initialize metadata
pressure_all = [];
current_all = [];

% Initialize jump component arrays
mean_dx_all = [];
mean_dy_all = [];
mean_dx_large = [];
mean_dy_large = [];
mean_dx_clustered = [];
mean_dy_clustered = [];

for f = 1:length(files)
    % Extract pressure and current from filename
    name = files(f).name;
    tokens = regexp(name, '(\d+)Pa (\d+\.?\d*)mA', 'tokens');
    if isempty(tokens)
        continue;
    end
    pressure = str2double(tokens{1}{1});
    current = str2double(tokens{1}{2});

    % Open the figure invisibly
    figPath = fullfile(figDir, name);
    fig = openfig(figPath, 'invisible');
    lines = findall(fig, 'Type', 'Line');

    % Containers for this file's data
    dx_all = [];
    dy_all = [];
    dx_large = [];
    dy_large = [];
    dx_clustered = [];
    dy_clustered = [];

    % Process each trajectory line
    for i = 1:length(lines)
        x = get(lines(i), 'XData')';
        y = get(lines(i), 'YData')';

        if length(x) < 2
            continue;
        end

        dx = diff(x);
        dy = diff(y);
        steps = sqrt(dx.^2 + dy.^2);
        step_mean = mean(steps);
        threshold = 10 * step_mean;

        dx_all = [dx_all; dx];
        dy_all = [dy_all; dy];

        % Identify large jumps
        large_jump_idx = find(steps >= threshold);
        dx_large = [dx_large; dx(large_jump_idx)];
        dy_large = [dy_large; dy(large_jump_idx)];

        % Identify clustered large jumps as grouped displacements
        clustered_dx_vals = [];
        clustered_dy_vals = [];

        if ~isempty(large_jump_idx)
            large_jump_idx = sort(large_jump_idx);  % Ensure increasing order

            current_cluster = large_jump_idx(1);

            for j = 2:length(large_jump_idx)
                if large_jump_idx(j) - current_cluster(end) <= 10
                    current_cluster(end+1) = large_jump_idx(j);
                else
                    clustered_dx_vals(end+1) = sum(dx(current_cluster));
                    clustered_dy_vals(end+1) = sum(dy(current_cluster));
                    current_cluster = large_jump_idx(j);
                end
            end

            % Final cluster
            if ~isempty(current_cluster)
                clustered_dx_vals(end+1) = sum(dx(current_cluster));
                clustered_dy_vals(end+1) = sum(dy(current_cluster));
            end
        end

        % Store cluster-level results
        dx_clustered = [dx_clustered; clustered_dx_vals(:)];
        dy_clustered = [dy_clustered; clustered_dy_vals(:)];
    end

    % Save values
    pressure_all(end+1) = pressure;
    current_all(end+1) = current;

    mean_dx_all(end+1) = mean(abs(dx_all));
    mean_dy_all(end+1) = mean(abs(dy_all));
    mean_dx_large(end+1) = mean(abs(dx_large));
    mean_dy_large(end+1) = mean(abs(dy_large));
    mean_dx_clustered(end+1) = mean(abs(dx_clustered));
    mean_dy_clustered(end+1) = mean(abs(dy_clustered));

    close(fig);
end

% --- Force all to column vectors ---
pressure_all = pressure_all(:);
current_all = current_all(:);
mean_dx_all = mean_dx_all(:);
mean_dy_all = mean_dy_all(:);
mean_dx_large = mean_dx_large(:);
mean_dy_large = mean_dy_large(:);
mean_dx_clustered = mean_dx_clustered(:);
mean_dy_clustered = mean_dy_clustered(:);

% --- Metrics for correlation analysis ---
metrics = {
    mean_dx_all, 'Mean dx';
    mean_dy_all, 'Mean dy';
    mean_dx_large, 'Large dx';
    mean_dy_large, 'Large dy';
    mean_dx_clustered, 'Clustered dx';
    mean_dy_clustered, 'Clustered dy';
};

fprintf('\n--- Pearson Correlation with Pressure & Current ---\n');
for i = 1:size(metrics,1)
    vec = metrics{i,1};
    label = metrics{i,2};

    % Debug print
    fprintf('\n[%s] Size check: vec = %s, pressure = %s\n', ...
        label, mat2str(size(vec)), mat2str(size(pressure_all)));

    if length(vec) ~= length(pressure_all)
        warning('Length mismatch for %s. Skipping.\n', label);
        continue;
    end

    [r_p, p_p] = corr(pressure_all, vec, 'Type', 'Pearson');
    [r_c, p_c] = corr(current_all, vec, 'Type', 'Pearson');

    fprintf('%s vs Pressure: r = %.3f, p = %.3f\n', label, r_p, p_p);
    fprintf('%s vs Current:  r = %.3f, p = %.3f\n', label, r_c, p_c);
end
%% Cluster Directionality Check with Condition Tracking (angle std < 45°)

% Initialize trackers
angle_spreads_per_cluster = [];
cluster_counts_per_file = struct();
detailed_cluster_list = [];
all_cluster_sizes = [];
all_cluster_lengths = [];

for f = 1:length(files)
    name = files(f).name;
    figPath = fullfile(figDir, name);
    fig = openfig(figPath, 'invisible');
    lines = findall(fig, 'Type', 'Line');

    % Extract pressure and current from filename
    tokens = regexp(name, '(\d+)Pa (\d+\.?\d*)mA', 'tokens');
    if isempty(tokens)
        continue;
    end
    pressure = str2double(tokens{1}{1});
    current = str2double(tokens{1}{2});

    count_under_45 = 0;

    for i = 1:length(lines)
        x = get(lines(i), 'XData')';
        y = get(lines(i), 'YData')';

        if length(x) < 3
            continue;
        end

        dx = diff(x);
        dy = diff(y);
        steps = sqrt(dx.^2 + dy.^2);
        angles = atan2d(dy, dx);
        step_mean = mean(steps);
        threshold = 10 * step_mean;

        large_jump_idx = find(steps >= threshold);
        if isempty(large_jump_idx)
            continue;
        end

        large_jump_idx = sort(large_jump_idx);
        current_cluster = large_jump_idx(1);

        for j = 2:length(large_jump_idx)
            if large_jump_idx(j) - current_cluster(end) <= 10
                current_cluster(end+1) = large_jump_idx(j);
            else
                if length(current_cluster) >= 2
                    cluster_angles = angles(current_cluster);
                    cluster_std = std(cluster_angles);
                    angle_spreads_per_cluster(end+1) = cluster_std;

                    cluster_dx = dx(current_cluster);
                    cluster_dy = dy(current_cluster);
                    cluster_lengths = sqrt(cluster_dx.^2 + cluster_dy.^2);
                    cluster_total_length = sum(cluster_lengths);
                    all_cluster_lengths(end+1) = cluster_total_length;

                    if cluster_std < 45
                        count_under_45 = count_under_45 + 1;
                        detailed_cluster_list(end+1).pressure = pressure;
                        detailed_cluster_list(end).current = current;
                        detailed_cluster_list(end).std = cluster_std;
                        detailed_cluster_list(end).file = name;
                        detailed_cluster_list(end).sum_abs_dx = sum(abs(cluster_dx));
                        detailed_cluster_list(end).sum_abs_dy = sum(abs(cluster_dy));
                        detailed_cluster_list(end).mean_angle = mean(atan2d(cluster_dy, cluster_dx));
                        detailed_cluster_list(end).cluster_size = length(current_cluster);
                        detailed_cluster_list(end).cluster_span = current_cluster(end) - current_cluster(1) + 1;
                        detailed_cluster_list(end).total_length = cluster_total_length;
                    end
                end
                current_cluster = large_jump_idx(j);
            end
        end

        % Final cluster
        if length(current_cluster) >= 2
            cluster_angles = angles(current_cluster);
            cluster_std = std(cluster_angles);
            angle_spreads_per_cluster(end+1) = cluster_std;

            cluster_dx = dx(current_cluster);
            cluster_dy = dy(current_cluster);
            cluster_lengths = sqrt(cluster_dx.^2 + cluster_dy.^2);
            cluster_total_length = sum(cluster_lengths);
            all_cluster_lengths(end+1) = cluster_total_length;

            if cluster_std < 45
                count_under_45 = count_under_45 + 1;
                detailed_cluster_list(end+1).pressure = pressure;
                detailed_cluster_list(end).current = current;
                detailed_cluster_list(end).std = cluster_std;
                detailed_cluster_list(end).file = name;
                detailed_cluster_list(end).sum_abs_dx = sum(abs(cluster_dx));
                detailed_cluster_list(end).sum_abs_dy = sum(abs(cluster_dy));
                detailed_cluster_list(end).mean_angle = mean(atan2d(cluster_dy, cluster_dx));
                detailed_cluster_list(end).cluster_size = length(current_cluster);
                detailed_cluster_list(end).cluster_span = current_cluster(end) - current_cluster(1) + 1;
                detailed_cluster_list(end).total_length = cluster_total_length;
            end
        end
    end

    close(fig);

    cluster_counts_per_file(end+1).pressure = pressure;
    cluster_counts_per_file(end).current = current;
    cluster_counts_per_file(end).file = name;
    cluster_counts_per_file(end).num_under_45 = count_under_45;
end

% --- Plot histogram of all angle spreads ---
figure;
histogram(angle_spreads_per_cluster, 'BinWidth', 5);
xlabel('Angular Std Dev within Cluster (°)');
ylabel('Count');
title('Directional Spread of Clustered Large Jumps');

% --- Summary Stats ---
num_total = length(angle_spreads_per_cluster);
num_under_45 = sum(angle_spreads_per_cluster < 45);
fprintf('\n--- Cluster Directionality Summary ---\n');
fprintf('Total clusters analyzed: %d\n', num_total);
fprintf('Clusters with std < 45°: %d\n', num_under_45);
fprintf('Mean angular std deviation: %.2f°\n', mean(angle_spreads_per_cluster));
fprintf('Median angular std deviation: %.2f°\n', median(angle_spreads_per_cluster));

% --- Aligned Cluster Breakdown ---
fprintf('\nClusters with angle std dev < 45°:\n');
for i = 1:length(detailed_cluster_list)
    d = detailed_cluster_list(i);
    dir_str = "mixed";
    if d.sum_abs_dx > 2 * d.sum_abs_dy
        dir_str = "mostly x";
    elseif d.sum_abs_dy > 2 * d.sum_abs_dx
        dir_str = "mostly y";
    end
    fprintf('→ %s | %.2f Pa | %.2f mA | Std Dev: %.2f° | Avg Angle: %.1f° | %s\n', ...
        d.file, d.pressure, d.current, d.std, d.mean_angle, dir_str);
end

% --- Final Averages ---
cluster_sizes = [detailed_cluster_list.cluster_size];
avg_cluster_size = mean(cluster_sizes);
fprintf('\nAverage number of large jumps per aligned cluster (std < 45°): %.2f\n', avg_cluster_size);

cluster_lengths = [detailed_cluster_list.total_length];
fprintf('Average total displacement length of aligned clusters (std < 45°): %.3f\n', mean(cluster_lengths));
fprintf('Average total displacement length of all clusters: %.3f\n', mean(all_cluster_lengths));


%% Analyze Behavior Before and After Large Jumps (Across All Files)

% Initialize containers
angle_diff_before_all = [];
angle_diff_after_all = [];
mag_ratio_before_all = [];
mag_ratio_after_all = [];

% Re-process each figure for before/after jump analysis
for f = 1:length(files)
    name = files(f).name;
    figPath = fullfile(figDir, name);
    fig = openfig(figPath, 'invisible');
    lines = findall(fig, 'Type', 'Line');

    for i = 1:length(lines)
        x = get(lines(i), 'XData')';
        y = get(lines(i), 'YData')';

        if length(x) < 3
            continue;
        end

        dx = diff(x);
        dy = diff(y);
        steps = sqrt(dx.^2 + dy.^2);
        angles = atan2d(dy, dx);
        step_mean = mean(steps);
        threshold = 10 * step_mean;

        large_jump_idx = find(steps >= threshold);

        for idx = large_jump_idx'
            if idx > 1
                angle_b = abs(angles(idx-1) - angles(idx));
                angle_b = mod(angle_b, 360);
                if angle_b > 180
                    angle_b = 360 - angle_b;
                end
                angle_diff_before_all(end+1) = angle_b;
                mag_ratio_before_all(end+1) = steps(idx-1) / step_mean;
            end

            if idx < length(angles)
                angle_a = abs(angles(idx+1) - angles(idx));
                angle_a = mod(angle_a, 360);
                if angle_a > 180
                    angle_a = 360 - angle_a;
                end
                angle_diff_after_all(end+1) = angle_a;
                mag_ratio_after_all(end+1) = steps(idx+1) / step_mean;
            end
        end
    end

    close(fig);
end

% --- Histograms ---
figure;
histogram(angle_diff_before_all, 'BinWidth', 10);
xlabel('Angle Difference (°) Before Large Jump');
ylabel('Count');
title('Direction Change Before Large Jumps');

figure;
histogram(angle_diff_after_all, 'BinWidth', 10);
xlabel('Angle Difference (°) After Large Jump');
ylabel('Count');
title('Direction Change After Large Jumps');

figure;
histogram(mag_ratio_before_all, 'BinWidth', 0.5);
xlabel('Step Size / Mean (Before)');
ylabel('Count');
title('Relative Step Magnitude Before Large Jumps');

figure;
histogram(mag_ratio_after_all, 'BinWidth', 0.5);
xlabel('Step Size / Mean (After)');
ylabel('Count');
title('Relative Step Magnitude After Large Jumps');

% --- Summary Statistics ---
fprintf('\n--- Jump Neighborhood Summary Across All Files ---\n');
fprintf('Mean angle change BEFORE: %.2f°\n', mean(angle_diff_before_all));
fprintf('Mean angle change AFTER:  %.2f°\n', mean(angle_diff_after_all));
fprintf('Mean magnitude ratio BEFORE: %.2f×\n', mean(mag_ratio_before_all));
fprintf('Mean magnitude ratio AFTER:  %.2f×\n', mean(mag_ratio_after_all));

%% Step Size Heatmap Across Trajectories (Time-Based Burst Visualization)

% Load one example figure (or loop through files if desired)
exampleFile = files(1).name;  % or pick any index
figPath = fullfile(figDir, exampleFile);
fig = openfig(figPath, 'invisible');
lines = findall(fig, 'Type', 'Line');

max_steps = 0;
step_matrix = [];

% Build step size matrix: rows = particles, columns = time steps
for i = 1:length(lines)
    x = get(lines(i), 'XData')';
    y = get(lines(i), 'YData')';
    
    if length(x) < 3
        continue;
    end

    dx = diff(x);
    dy = diff(y);
    steps = sqrt(dx.^2 + dy.^2);

    max_steps = max(max_steps, length(steps));
    step_matrix{i} = steps;  % store variable-length rows
end

% Pad all rows to equal length with NaN
step_array = nan(length(step_matrix), max_steps);
for i = 1:length(step_matrix)
    s = step_matrix{i};
    step_array(i, 1:length(s)) = s;
end

% Normalize rows (optional)
normalize = false;
if normalize
    step_array = step_array ./ mean(step_array, 2, 'omitnan');
end

% Plot heatmap
figure;
imagesc(step_array);
colormap('hot');  % or try 'parula', 'jet', 'turbo'
colorbar;
xlabel('Step Index (Time)');
ylabel('Trajectory Index');
title('Step Size Heatmap (Time-Based Bursts)');



%% Mean and Max Step Size for Each Pressure/Current Case

% Set figure directory
figDir = '/Users/bradleyandrew/Desktop/Research/PK-4 Experiment/Newest PK4 Data analysis and Ideas/track plots/Camera 1/';
files = dir(fullfile(figDir, '*1000above10frs.fig'));

% Initialize result table
results = [];

for f = 1:length(files)
    name = files(f).name;

    % Extract pressure and current from filename
    tokens = regexp(name, '(\d+)Pa (\d+\.?\d*)mA', 'tokens');
    if isempty(tokens)
        continue;
    end
    pressure = str2double(tokens{1}{1});
    current = str2double(tokens{1}{2});

    % Open figure invisibly
    figPath = fullfile(figDir, name);
    fig = openfig(figPath, 'invisible');
    lines = findall(fig, 'Type', 'Line');

    % Collect step sizes
    all_steps = [];
    for i = 1:length(lines)
        x = get(lines(i), 'XData')';
        y = get(lines(i), 'YData')';

        if length(x) < 2
            continue;
        end

        steps = sqrt(diff(x).^2 + diff(y).^2);
        all_steps = [all_steps; steps];
    end

    % Compute stats
    mean_step = mean(all_steps);
    max_step = max(all_steps);

    % Store
    results(end+1).pressure = pressure;
    results(end).current = current;
    results(end).mean = mean_step;
    results(end).max = max_step;

    close(fig);
end

% --- Print results for Python copy ---
fprintf('\n--- Mean and Max Step Sizes by Pressure/Current ---\n');
for r = 1:length(results)
    fprintf('{"pressure": %d, "current": %.2f, "mean": %.4f, "max": %.4f},\n', ...
        results(r).pressure, results(r).current, results(r).mean, results(r).max);
end
