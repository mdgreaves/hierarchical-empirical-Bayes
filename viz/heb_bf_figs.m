function heb_bf_figs()

% =========================================================================
% heb_bf_figs.m
% =========================================================================
% Reproduces manuscript data figures based on compact derived inputs in
% fig_data.mat:
%   - Fig. 3  : Bayesian model-averaged prior-variance transformation
%   - Fig. S3 : Test-retest reliability
%   - Fig. 4  : Out-of-sample validation (session 1)
%   - Fig. S4 : Out-of-sample validation (session 2)
%
% Inputs
% ------
% fig_data.mat with variables:
%   explore, retest_holdout, schaefer_networks
%
% =========================================================================
% Load data and figure metadata
% =========================================================================
data_file   = './data/fig_data.mat';
load(data_file, "explore", "retest_holdout", "schaefer_networks");
fig_names   = {['Fig. 3: Bayesian model-averaged prior-variance', ...
    ' transformation for 17 brain networks'];
    ['Fig. S3: Test-retest reliability of hierarchical empirical', ...
    ' Bayes models'];
    ['Fig. 4: Out-of-sample validation of hierarchical empirical', ...
    ' Bayes models (session-1 results)'];
    ['Fig. S4: Out-of-sample validation of hierarchical empirical', ...
    ' Bayes models (session-2 results)']};
font_size   = 20;

% =========================================================================
% Exploration and face validation (Fig. 3)
% =========================================================================
fig = figure;
tiledlayout(fig, 4, 5);

set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], ...
    'Name', fig_names{1}, 'NumberTitle', 'off');

% Add common axis labels for Fig. 3
annotation('textbox', [0.093, 0.3, 0.25, 0.30], ...
    'String', 'Prior variance of effective connectivity', ...
    'Rotation', 90, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', font_size, ...
    'EdgeColor', 'none');
annotation('textbox', [0.325, 0.06, 0.2, 0.01], ...
    'String', 'Structural connectivity', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', font_size, ...
    'EdgeColor', 'none');

% Continuous functions with confidence envelope
fun_plot_tiles  = [1:5, 6:9, 11:14, 16:19];
x_mesh          = 0:0.01:1;
log_BFs         = zeros(length(schaefer_networks), 1);

% A subplots
for i_network = 1:length(schaefer_networks)
    nexttile(fun_plot_tiles(i_network), [1, 1]);
    rgb_code = network_colour(i_network);
    [~, se_upper, se_lower] = plot_bma_se(x_mesh, ...
        explore.bf_ndgrid{i_network},...
        explore.hp_ndgrid.alphas, explore.hp_ndgrid.betas,...
        explore.hp_bma.alphas{i_network},...
        explore.hp_bma.betas{i_network}, 1/2, rgb_code);

    % Adjust plot for readability
    y_upper_at_x1 = interp1(x_mesh, se_upper, 1, 'linear', 'extrap');
    y_lower_at_x1 = interp1(x_mesh, se_lower, 1, 'linear', 'extrap');
    y_range = y_upper_at_x1-y_lower_at_x1;
    if ((explore.hp_bma.alphas{i_network}+ ...
            explore.hp_bma.betas{i_network})- ...
            explore.hp_bma.alphas{i_network}) > 1e-3*5
    ylim([min(explore.hp_bma.alphas{i_network}, ...
        y_lower_at_x1-(y_range/3)),...
        y_upper_at_x1+(y_range/3)]);
    end
    hold on;
    xticks([0, 1])
    yticks([explore.hp_bma.alphas{i_network},...
        (explore.hp_bma.alphas{i_network}+ ...
        explore.hp_bma.betas{i_network})]);
    title(schaefer_networks{i_network, 2}, 'FontSize', font_size-1,...
        'HorizontalAlignment', 'center', 'FontWeight','normal')
    axesObjects = findall(gcf, 'Type', 'axes');
    set(axesObjects, 'FontSize', font_size);
    
    if i_network == 1
    text(axesObjects, min(get(axesObjects, 'xlim')),...
        max(get(axesObjects, 'ylim')), 'A', ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
     'FontSize', font_size+2, 'Units', 'normalized', ...
     'Position', [0.06 1.40], 'FontWeight', 'bold');
    end

    % Normalised log BFs for (final) bar plot
    log_BFs(i_network) = explore.BF{i_network}/100;
end

% B subplot
ax      = nexttile(10, [3, 1]);
barObj  = bar(ax, flip(log_BFs), 'FaceColor', 'flat', 'Horizontal', 'on');
set(barObj, 'EdgeColor', 'none');

% Adjust x-ticks and labels
set(ax, 'ytick', 1:length(schaefer_networks), 'yticklabel', ...
    flip(schaefer_networks(:,2)));
xline(ax, 3, '--r', 'LineWidth', 1.5);
ytickangle(ax, 45);

% Colouring the bars using the network colours
for i_network = 1:length(schaefer_networks)
    i_network_fliped = flip(1:length(schaefer_networks));
    barObj.CData(i_network, :) = network_colour(...
        i_network_fliped(i_network));
end

% Additional adjustments
xlabel(ax, {'Group-level log-Bayes factor', '(scaled by sample size)'},...
    'FontSize', font_size-1);
set(ax, 'FontSize', font_size-1);
box(ax, 'off');
    text(ax, min(get(ax, 'xlim')), max(get(ax, 'ylim')), 'B', ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
     'FontSize', font_size+2, 'Units', 'normalized', ...
     'Position', [0.06 1.07], 'FontWeight', 'bold');

% =========================================================================
% Test-retest and out-of-sample validation (Figs. 4, S3–S4)
% =========================================================================
for i_validation = 1:size(retest_holdout.vPEB,1)

    % Create a new figure and tiled layout
    fig = figure;
    t = tiledlayout(fig, 3, 6);
    set(fig, 'Units', 'normalized', 'OuterPosition', [0 1/2 1 1/1.3], ...
    'Name', fig_names{i_validation+1}, 'NumberTitle', 'off');

    % Adding common X and Y labels
    ylabel(t,'Sample-size-scaled log-Bayes factor', 'FontSize', font_size);
    annotation('textbox', [0.325, 0.06, 0.2, 0.01], ...
    'String', 'Subject- and group-level result', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', font_size, ...
    'EdgeColor', 'none');

    for i_network = 1:length(schaefer_networks)
        rgb_code = network_colour(i_network);
        nexttile;
        norm_logBF = retest_holdout.Fdiff{...
            i_validation, i_network}/...
            (length(retest_holdout.iFs{...
            i_validation, i_network}));
        bar(length(retest_holdout.iFs{...
            i_validation, i_network})/2,...
            norm_logBF,...
            'FaceColor', rgb_code,...
            'FaceAlpha', 0.3, 'LineStyle', 'none', 'BarWidth',...
            length(retest_holdout.iFs{i_validation, i_network}));
        hold on;

        % Split data into positive and negative components
        subject_log_BF =...
            retest_holdout.viFs{i_validation, i_network} -...
            retest_holdout.iFs{i_validation, i_network};

        % Plot positive values
        bar(subject_log_BF, 'FaceColor', rgb_code, 'LineStyle', 'none');

        % Additional plotting (e.g., yline, xticks, titles)
        yline(3, 'r--', 'LineWidth', 1.5);
        padding = length(retest_holdout.iFs{...
            i_validation, i_network})/10;
        xlim([-padding, length(retest_holdout.iFs{...
            i_validation, i_network})+padding])
        xticks([1, length(retest_holdout.iFs{...
            i_validation, i_network})]);
        maxVal = max(subject_log_BF);
        if round(norm_logBF,3) < maxVal*0.15
            yticks([round(norm_logBF,3), maxVal]);
        else
            yticks([0, round(norm_logBF,3), maxVal]);
        end

        % Add title
        title(schaefer_networks{i_network, 2}, 'FontSize', font_size-1, ...
            'FontWeight', 'normal');
        axesObjects = findall(gcf, 'Type', 'axes');
        set(axesObjects, 'FontSize', font_size-1);

        if i_network == 1
            text(axesObjects, min(get(axesObjects, 'xlim')),...
                max(get(axesObjects, 'ylim')), 'A', ...
                'VerticalAlignment', 'top', 'HorizontalAlignment',...
                'right', 'FontSize', font_size+2, 'Units',...
                'normalized', 'Position', [0.06 1.40],...
                'FontWeight', 'bold');
        end
    end

    % Adding histogram of effective connectivity (in Hz)
    ax = nexttile;
    for i_network = 1:length(schaefer_networks)
        A = full(retest_holdout.vPEB{i_validation,...
            i_network}.Ep);
        n = sqrt(numel(A));
        A = reshape(A, n, n);
        A(logical(eye(n))) = -0.5 * exp(A(logical(eye(n))));
        Ep = A(:);

        histogram(ax, Ep, 'FaceColor', [0.5, 0.5, 0.5],...
            'LineStyle', 'none');
        hold on;
    end

    % Additional adjustments
    xlabel(ax, {'MAP effective', 'connectivity (Hz)'}, 'FontSize', ...
        font_size-1);
    ylabel(ax, 'Frequency', 'FontSize', font_size-1);
    set(ax, 'FontSize', font_size-1);
    box(ax, 'off');
    text(ax, min(get(ax, 'xlim')), max(get(ax, 'ylim')), 'B', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
        'FontSize', font_size+2, 'Units', 'normalized', ...
        'Position', [0.06 1.3], 'FontWeight', 'bold');
end

end

% =========================================================================
% Helper Functions
% =========================================================================

function rgb_code = network_colour(network_number)
% Colour map taken from Schaefer labels
cmap = [234, 147, 35;...
    140, 51, 76;...
    123, 141, 177;...
    251, 254, 2;...
    210, 62, 80;...
    4, 1, 132;...
    74, 155, 61;...
    0, 118, 16;...
    220, 248, 165;...
    122, 135, 51;...
    196, 58, 251;...
    255, 152, 241;...
    74, 131, 177;...
    47, 205, 162;...
    16, 48, 255;...
    124, 19, 134;...
    25, 0, 0];

% Darken selected for greater accessibility
cmap(4, :) = cmap(4, :) * 0.85;
cmap(9, :) = cmap(9, :) * 0.85;

cmap = cmap./255;
rgb_code = cmap(network_number, :);
end

function [f_bma, se_upper, se_lower] = plot_bma_se(x_values, Fs, ...
    alphas, betas, alpha_mean,...
    beta_mean, max_var, rgb_code)

% Create ND grid of parameters governing data-to-variance mapping
[Alphas, Betas] = ndgrid(alphas, betas);

% Identify valid combinations
epsilon = 1e-5;
valid_indices = (Alphas >= epsilon) & (Alphas <= max_var)...
    & (Alphas + Betas >= epsilon) & (Alphas + Betas...
    <= max_var);

% Extract valid combinations
valid_alphas = Alphas(valid_indices);
valid_betas  = Betas(valid_indices);

% Obtain weights
Fs_matrix = cell2mat(Fs);

% Compute evidence with numerical stability
Fs_max = max(Fs_matrix(:));
evidence = exp(Fs_matrix - Fs_max);
total_evidence = sum(evidence(:));
if total_evidence == 0
    warning('Total evidence is zero. Using uniform weights.');
    weights = ones(size(Fs_matrix)) / numel(Fs_matrix);
else
    weights = evidence / total_evidence;
end

% Calculate variances of alpha and beta
alpha_var = sum(weights .* (valid_alphas - alpha_mean).^2);
beta_var = sum(weights .* (valid_betas - beta_mean).^2);

% Initialize BMA function and variance
f_bma = alpha_mean + x_values .* beta_mean;
f_var = alpha_var + (x_values.^2) .* beta_var;

% Compute confidence envelope
se_upper = f_bma + sqrt(f_var);
se_lower = f_bma - sqrt(f_var);

% Plot the BMA function with confidence envelope
fill([x_values, fliplr(x_values)], [se_upper,...
    fliplr(se_lower)],...
    rgb_code, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;
plot(x_values, f_bma, 'Color', rgb_code, 'LineWidth', 3.5);

end