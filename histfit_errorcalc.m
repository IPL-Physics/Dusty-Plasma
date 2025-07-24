%%
% -------------------------------
% Load FIG file and extract data
% -------------------------------
fileNames = {'70Pa 07mA hist Vx.fig','70Pa 07mA hist Vy.fig'};

% Store stats
R_all = []; RMSE_all = []; chi2_all = []; label_all = []; bin_counts = [];

% -----------------------------------------------
% Declare Stats function
% -----------------------------------------------
function [R2, RMSE, chi2] = fit_stats(y_actual, y_pred, numParams)
    valid = y_pred > 0;
    SS_res = sum((y_actual - y_pred).^2);
    SS_tot = sum((y_actual - mean(y_actual)).^2);
    R2 = 1 - SS_res / SS_tot;
    RMSE = sqrt(mean((y_actual - y_pred).^2));
    chi2_val = sum(((y_actual(valid) - y_pred(valid)).^2) ./ y_pred(valid));
    chi2 = chi2_val / (sum(valid) - numParams);
end


for i = 1:length(fileNames)
    fig = openfig(fileNames{i}, 'invisible');  % open silently
    ax = gca;
    allObjs = findall(ax, 'Type', 'line');     % all curves
    barObj = findall(ax, 'Type', 'bar');       % histogram

    % Histogram data
    x_data = barObj.XData;
    y_data = barObj.YData;

    % Fit curves
    blackIdx   = find(arrayfun(@(h) isequal(h.Color, [0 0 0]), allObjs));
    magentaIdx = find(arrayfun(@(h) isequal(h.Color, [1 0 1]), allObjs));
    greenIdx   = find(arrayfun(@(h) isequal(h.Color, [0 1 0]), allObjs));

    if isempty(blackIdx)
        warning('Black (Maxwellian) curve not found in %s', fileNames{i});
    end
    if isempty(magentaIdx)
        warning('Magenta (single-q) curve not found in %s', fileNames{i});
    end
    if isempty(greenIdx)
        warning('Green (two-q) curve not found in %s', fileNames{i});
    end

    % Define x overlap region for interpolation
    all_x_ranges = [x_data, ...
                    allObjs(blackIdx).XData, ...
                    allObjs(magentaIdx).XData, ...
                    allObjs(greenIdx).XData];
    x_min = max([min(all_x_ranges)]);
    x_max = min([max(all_x_ranges)]);

    validIdx = (x_data >= x_min) & (x_data <= x_max);
    x_clipped = x_data(validIdx);
    y_actual_clipped = y_data(validIdx);

    % Interpolate and store stats
    if ~isempty(blackIdx)
        y_fit_black = interp1(allObjs(blackIdx).XData, allObjs(blackIdx).YData, x_clipped, 'linear', 'extrap');
        [R, RMSE, chi2] = fit_stats(y_actual_clipped, y_fit_black, 2);
        R_all(end+1) = R; RMSE_all(end+1) = RMSE; chi2_all(end+1) = chi2;
        label_all{end+1} = 'Black';
        bin_counts(end+1) = length(x_data);
        fprintf('\n-- Black (Maxwellian) fit -- %s\n', fileNames{i});
        fprintf('R^2 = %.4f, RMSE = %.4f, Chi^2/dof = %.4f\n', R, RMSE, chi2);
    end

    if ~isempty(magentaIdx)
        y_fit_magenta = interp1(allObjs(magentaIdx).XData, allObjs(magentaIdx).YData, x_clipped, 'linear', 'extrap');
        [R, RMSE, chi2] = fit_stats(y_actual_clipped, y_fit_magenta, 2);
        R_all(end+1) = R; RMSE_all(end+1) = RMSE; chi2_all(end+1) = chi2;
        label_all{end+1} = 'Magenta';
        bin_counts(end+1) = length(x_data);
        fprintf('\n-- Magenta (Single-q) fit -- %s\n', fileNames{i});
        fprintf('R^2 = %.4f, RMSE = %.4f, Chi^2/dof = %.4f\n', R, RMSE, chi2);
    end

    if ~isempty(greenIdx)
        y_fit_green = interp1(allObjs(greenIdx).XData, allObjs(greenIdx).YData, x_clipped, 'linear', 'extrap');
        [R, RMSE, chi2] = fit_stats(y_actual_clipped, y_fit_green, 4);
        R_all(end+1) = R; RMSE_all(end+1) = RMSE; chi2_all(end+1) = chi2;
        label_all{end+1} = 'Green';
        bin_counts(end+1) = length(x_data);
        fprintf('\n-- Green (Two-q) fit -- %s\n', fileNames{i});
        fprintf('R^2 = %.4f, RMSE = %.4f, Chi^2/dof = %.4f\n', R, RMSE, chi2);
    end

    close(fig);
end

%%
fileNames = {'70Pa 07mA hist Vx.fig','70Pa 07mA hist Vy.fig'};
bin_range = 10:10:2000;

for f = 1:length(fileNames)
    clear barObj allLines x_data y_data raw_data fit_curve y_fit  % clear figure-specific vars
    % ---------------------------
    % Open the original .fig file
    % ---------------------------
    fig = openfig(fileNames{f}, 'invisible');  % open silently
    ax = gca;
    
    % ---------------------------
    % Extract histogram bar object robustly
    % ---------------------------
    allGraphics = findall(fig, '-property', 'Type'); % look in full figure
    barObj = [];
    for g = 1:length(allGraphics)
        if strcmp(get(allGraphics(g), 'Type'), 'bar') && ...
           isprop(allGraphics(g), 'XData') && isprop(allGraphics(g), 'YData')
            barObj = allGraphics(g);
            break;
        end
    end
    
    if isempty(barObj)
        warning('No valid histogram bar object found in %s — skipping.', fileNames{f});
        continue;
    end
    
    x_data = barObj.XData;
    y_data = barObj.YData;
    bin_width = mode(diff(x_data));
    fprintf('%s had ~%d bins\n', fileNames{f}, length(x_data));

    % ---------------------------
    % Reconstruct raw data
    % ---------------------------
    raw_data = repelem(x_data, round(y_data));
    
    % ---------------------------
    % Detect fit lines by color
    % ---------------------------
    allLines = findall(ax, 'Type', 'line');
    
    blackIdx   = find(arrayfun(@(h) isequal(h.Color, [0 0 0]), allLines));
    magentaIdx = find(arrayfun(@(h) isequal(h.Color, [1 0 1]), allLines));
    greenIdx   = find(arrayfun(@(h) isequal(h.Color, [0 1 0]), allLines));

    fits = {};
    fitColors = {};
    fitLabels = {};
    markerSizes = [];

    if ~isempty(blackIdx)
        fits{end+1} = blackIdx;
        fitColors{end+1} = [0 0 0];
        fitLabels{end+1} = 'Maxwellian';
        markerSizes(end+1) = 30;
    end
    if ~isempty(magentaIdx)
        fits{end+1} = magentaIdx;
        fitColors{end+1} = [1 0 1];
        fitLabels{end+1} = 'Single-q';
        markerSizes(end+1) = 50;
    end
    if ~isempty(greenIdx)
        fits{end+1} = greenIdx;
        fitColors{end+1} = [0 1 0];
        fitLabels{end+1} = 'Two-q';
        markerSizes(end+1) = 70;
    end

    numFits = length(fits);
    if numFits == 0
        warning('No fits found in %s — skipping.', fileNames{f});
        continue;
    end

    % Initialize results
    R_all = zeros(length(bin_range), numFits);
    RMSE_all = zeros(length(bin_range), numFits);
    chi2_all = zeros(length(bin_range), numFits);

    for b = 1:length(bin_range)
        nbins = bin_range(b);
        [y_rehist, edges] = histcounts(raw_data, nbins);
        x_centers = (edges(1:end-1) + edges(2:end)) / 2;

        for j = 1:numFits
            idx = fits{j};
            fit_curve = allLines(idx);
            if isempty(fit_curve.XData)
                continue;
            end

            y_fit = interp1(fit_curve.XData, fit_curve.YData, x_centers, 'linear', 'extrap');
            [R, RMSE, chi2] = fit_stats(y_rehist, y_fit, 3);
            R_all(b,j) = R;
            RMSE_all(b,j) = RMSE;
            chi2_all(b,j) = chi2;
        end
    end

    % ---------------------------------------
    % Plotting Fit Stats for This .fig File
    % ---------------------------------------
    figure('Name', ['Fit Stats for ' fileNames{f} ' vs. Bin Count'], ...
           'Color', 'w', 'Position', [100 100 1400 500]);

    stats = {R_all, RMSE_all, chi2_all};
    titles = {'R^2', 'RMSE', '\chi^2/dof'};
    ylabels = {'R^2', 'RMSE', '\chi^2/dof'};
    ylimits = {[0.7 1], [0 10], [0 3]};  % Custom Y axis limits

    for s = 1:3
        subplot(1,3,s); hold on;
        for j = 1:numFits
            scatter(bin_range, stats{s}(:,j), markerSizes(j), ...
                'MarkerFaceColor', fitColors{j}, ...
                'MarkerEdgeColor', 'k');
        end
        xlabel('# Bins'); 
        ylabel(ylabels{s}); 
        title(titles{s});

        xlim([400 1000]);              % Set consistent X-axis range
        ylim(ylimits{s});              % Custom Y-axis range per stat

        if s == 3
            legend(fitLabels, 'Location', 'northeast');
        end
        grid on;
    end


    close(fig);
end


