
% =========================================================================
% heb_brain_fig.m
% =========================================================================
% Reproduces manuscript Fig. 5 using compact derived inputs in ./data.
%
% Panel A: cortical surface projection of network-specific scale
%          hyperparameters (beta).
% Panel B: ranked bar plot of network-specific beta values.
% Panel C: relationship between network-specific beta values and mean
%          principal-gradient position.
%
% Required inputs in ./data
% --------------------------
% fig_data.mat containing variables: explore, useful, schaefer_networks
% Schaefer annotation files: lh/rh *.annot
% Cortical surfaces: lh/rh *.gii
%
% External function requirements
% ------------------------------
% The script uses `gifti` and `read_annotation`, which must be available
% on the MATLAB path. `gifti` is commonly provided by SPM12 (or other
% GIfTI toolbox distributions); `read_annotation` is typically provided
% via FreeSurfer MATLAB utilities.
% =========================================================================
% Load data and figure metadata
% =========================================================================
data_dir    = './data';
data_file   = fullfile(data_dir, 'fig_data.mat');
parc_info   = 'Schaefer2018_200Parcels_17Networks_order_info.txt';
load(data_file, "explore", "useful", "schaefer_networks");
fig_name    = ['Fig. 5: Structural modulation of ', ...
    'effective connectivity across large-scale cortical networks'];
font_size   = 20;

% Validate required external functions (toolbox-agnostic check)
required_fns = {'gifti', 'read_annotation'};
missing_fns = required_fns(~cellfun(@(f) ...
    (exist(f, 'file') == 2) || (exist(f, 'class') == 8), required_fns));
if ~isempty(missing_fns)
    error(['Missing required function(s): %s. Add them to your MATLAB', ...
        ' path. `gifti` is commonly provided by SPM12 (or another', ...
        ' GIfTI toolbox), and `read_annotation` by FreeSurfer', ...
        ' MATLAB utilities.'], strjoin(missing_fns, ', '));
end

% Set up basic figure
fig = figure;
set(fig, 'Color', 'white', 'Units', 'normalized', 'Position', ...
    [0, 1, 1, 2/3], 'Name', fig_name, 'NumberTitle', 'off')

% Load network order and hyperparameters
info    = readcell(fullfile(data_dir, parc_info));
info    = info(1:2:length(info), 1);
betas   = cell2mat(explore.hp_bma.betas);

% =========================================================================
% SUBPLOT C: Scatterplot showing the relationship between network-specific 
% scale hyperparameters and the mean position of each network along the 
% principal gradient of functional connectivity 
% =========================================================================
ax_c    = subplot(2,4,[7,8]);
x       = zscore(useful.mean_pgfc);
y       = betas;

% Correlation and confidence interval
[RHO, PVAL] = corr(x', y');
n           = length(x);
z           = atanh(RHO);                  
SE          = 1 / sqrt(n - 3);           
z_crit      = 1.96;                   

% Confidence interval in z-space
z_lower = z - z_crit * SE;
z_upper = z + z_crit * SE;

% Convert back to r-space
r_lower = tanh(z_lower);
r_upper = tanh(z_upper);

% Plotting the correlation
for i_network   = 1:length(schaefer_networks)
    rgb         = network_colour(i_network);
    scatter(ax_c, x(i_network), y(i_network),...
        'filled', 'MarkerFaceColor', rgb, 'SizeData', 150);
    hold on;
end

% Add least-squares regression line in gray dashed style
coeffs  = polyfit(x, y, 1);
x_fit   = linspace(min(x), max(x), 100);
y_fit   = polyval(coeffs, x_fit);
plot(ax_c, x_fit, y_fit, '--', 'Color', [0.4 0.4 0.4], ...
    'LineWidth', 1.5);

xlabel({'Mean position (principal gradient)'});
ylabel('Hyperparameter: ');
ylabel('Hyperparameter: ', 'FontSize', font_size);

% Add β as via annotation (position is [x y w h] in normalized units)
annotation('textbox', [0.15105, 0.518, 0.05, 0.05], ...
    'String', 'β', ...
    'FontName', 'Times New Roman', ...
    'FontSize', font_size, ...
    'FontAngle', 'italic', ...
    'LineStyle', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'Rotation', 90);
ci_str = sprintf('95%% CI [%.2f, %.2f]', r_lower, r_upper);

% Place the correlation information in the upper right-hand corner
xlims = xlim; 
ylims = ylim; 
text(xlims(2), ylims(2), ['{\it r}', sprintf(' = %.2f, ', RHO), ...
    ci_str], 'VerticalAlignment', 'top', 'HorizontalAlignment', ...
    'right', 'FontSize', font_size, 'Units', 'normalized', ...
     'Position', [1, 1.10]);
set(ax_c, 'FontSize', font_size);
set(ax_c, 'Position', [0.7 0.1100 0.25 0.4])

% =========================================================================
% SUBPLOT B: Bar plot of scale hyperparameters ranked by magnitude
% =========================================================================
ax_b = subplot(2,4,[5,6]);

[betas_sorted, sortIdx] = sort(betas, 'ascend');
networks_sorted = schaefer_networks(sortIdx, 2);

% Use the predefined axis 'ax' for plotting
barObj = bar(ax_b, betas_sorted, 'FaceColor', 'flat');
set(barObj, 'EdgeColor', 'none');

% Adjust x-ticks and labels
set(ax_b, 'xtick', 1:length(networks_sorted), 'xticklabel', ...
    networks_sorted);
xtickangle(ax_b, 45);

% Coloring the bars using the cool colormap
colormap(ax_b, 'cool');  
nColors     = length(betas_sorted);
cMap        = colormap(ax_b);

% Normalize betas_sorted to the range [0, 1]
betas_min           = min(betas_sorted);
betas_max           = max(betas_sorted);
betas_normalized    = (betas_sorted - betas_min) / (betas_max - betas_min);

% Map normalized beta values to colormap indices
colorIndices = round(1 + betas_normalized * (size(cMap, 1) - 1));

% Assign colors to each bar individually based on normalized beta values
for i = 1:nColors
    barObj.CData(i, :) = cMap(colorIndices(i), :);
end

% Optional: Additional adjustments
box(ax_b, 'off');
axis(ax_b, 'tight');

% Adjust fontsize for subplot C
ylabel('Hyperparameter: ', 'FontSize', font_size);
ylim([0, 0.15]);
set(ax_b, 'xticklabel', networks_sorted, 'FontSize', font_size); 
xtickangle(ax_b, 45);

% Add β as an annotation textbox (position is [x y w h] in normalized units)
annotation('textbox', [0.6108, 0.47, 0.05, 0.05], ...
    'String', 'β', ...
    'FontName', 'Times New Roman', ...
    'FontSize', font_size, ...
    'FontAngle', 'italic', ...
    'LineStyle', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'Rotation', 90);  % Rotate to match y-axis label orientation

% =========================================================================
% SUBPLOT A: Cortical projection of network-specific scale hyperparameters
% =========================================================================
ax_a1 = subplot(2,4,1);
heb_brain_viz(betas, schaefer_networks(:,1), ax_a1, 'p2');
ax_b1 = subplot(2,4,2);
heb_brain_viz(betas, schaefer_networks(:,1), ax_b1, 'p1');
ax_c1 = subplot(2,4,3);
heb_brain_viz(betas, schaefer_networks(:,1), ax_c1, 'p3');
ax_d1 = subplot(2,4,4);
heb_brain_viz(betas, schaefer_networks(:,1), ax_d1, 'p4');
set(ax_a1, 'Position', [0.0 0.5838 0.4 0.4])
set(ax_b1, 'Position', [0.2 0.5838 0.4 0.4])
set(ax_c1, 'Position', [0.4 0.5838 0.4 0.4])
set(ax_d1, 'Position', [0.6 0.5838 0.4 0.4])
c = colorbar(ax_d1, 'eastoutside');
colormap(ax_d1, 'cool');
c.Limits = [0, 1]; 
c.Ticks = [0, 1];
c.FontSize = font_size;
c.Label.String = 'arb.';
c.Label.FontSize = font_size;
c.Label.Position = [1.25, 0.5, 0];
set(c, 'Position', [0.9 0.689 0.009 0.25])
set(ax_b, 'Position', [0.18 0.27 0.3628 0.33])
set(ax_c, 'Position', [0.64 0.19 0.2500 0.4000])

% Add annotations
annotation('textbox', [0.12, 0.93, 0.05, 0.05], ...
    'String', 'A', ...
    'FontWeight', 'bold', ...
    'FontSize', font_size, ...
    'LineStyle', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle');
annotation('textbox', [0.12, 0.63, 0.05, 0.05], ...
    'String', 'B', ...
    'FontWeight', 'bold', ...
    'FontSize', font_size, ...
    'LineStyle', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle');
annotation('textbox', [0.58, 0.615, 0.05, 0.05], ...
    'String', 'C', ...
    'FontWeight', 'bold', ...
    'FontSize', font_size, ...
    'LineStyle', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle');

% =========================================================================
% Helper Functions
% =========================================================================
function heb_brain_viz(params, networks, ax_ob, ax_select)

% Sort parameters and corresponding networks
[params_sorted, sortIdx] = sort(params, 'ascend');
networks_sorted = networks(sortIdx);

% Coloring the bars using the cool colormap
cMap = colormap('cool');

% Normalize betas_sorted to the range [0, 1]
params_min = min(params_sorted);
params_max = max(params_sorted);
params_normalized = (params_sorted - params_min) / ...
    (params_max - params_min);

% Map normalized beta values to colormap indices
colorIndices = round(1 + params_normalized * (size(cMap, 1) - 1));

% Load atlas labels and cortical meshes from local data directory
addpath('./data');
% Paths to annotation and pial surface files
lh_annot_file = 'lh.Schaefer2018_200Parcels_17Networks_order.annot';
rh_annot_file = 'rh.Schaefer2018_200Parcels_17Networks_order.annot';

lh_pial_file = 'lh.pial.gii';
rh_pial_file = 'rh.pial.gii';

% Read the left hemisphere annotation
[vertices_lh, label_lh, colortable_lh] = read_annotation(lh_annot_file);

% Read the right hemisphere annotation
[vertices_rh, label_rh, colortable_rh] = read_annotation(rh_annot_file);

% Load the left hemisphere pial surface
lh_pial = gifti(lh_pial_file);
vertices_lh_pial = lh_pial.vertices;
faces_lh = lh_pial.faces;

% Load the right hemisphere pial surface
rh_pial = gifti(rh_pial_file);
vertices_rh_pial = rh_pial.vertices;
faces_rh = rh_pial.faces;

% Create color maps for vertices (initialize to grey)
vertex_colors_lh = repmat([0.7, 0.7, 0.7], size(vertices_lh_pial, 1), 1);
vertex_colors_rh = repmat([0.7, 0.7, 0.7], size(vertices_rh_pial, 1), 1);

% Assign network-specific colors to parcel vertices (both hemispheres)
for i_network = 1:length(networks_sorted)
    regions = networks_sorted{i_network};
    % Find the indices of the regions in the colortables that belong to
    % the specified network
    selected_region_indices_lh = ...
        find(contains(colortable_lh.struct_names, regions));
    selected_region_indices_rh = ...
        find(contains(colortable_rh.struct_names, regions));

    % Get the corresponding labels for the selected regions
    selected_labels_lh = ...
        colortable_lh.table(selected_region_indices_lh, 5);
    selected_labels_rh = ...
        colortable_rh.table(selected_region_indices_rh, 5);

    % Map each vertex to the corresponding color in the colortable for
    % the left hemisphere
    for i = 1:length(vertices_lh)
        if ismember(label_lh(i), selected_labels_lh)
            vertex_colors_lh(i, :) = cMap(colorIndices(i_network), :);
        end
    end

    % Map each vertex to the corresponding color in the colortable for
    % the right hemisphere
    for i = 1:length(vertices_rh)
        if ismember(label_rh(i), selected_labels_rh)
            vertex_colors_rh(i, :) = cMap(colorIndices(i_network), :);
        end
    end
end

% Render the requested hemisphere/view
% Plot the selected panel
set(ax_ob);

switch ax_select
    case 'p2'
        % Left Hemisphere - Medial View
        patch('Vertices', vertices_lh_pial, 'Faces', faces_lh, ...
            'FaceVertexCData', vertex_colors_lh, ...
            'FaceColor', 'interp', 'EdgeColor', 'none');
        lighting gouraud;
        camlight('headlight');
        view([-90, 0]);
        lighting flat;
        camlight('headlight');
        axis equal off;
        xlim([-100, 100]); ylim([-110, 90]); zlim([-100, 100]);
        material dull;

    case 'p4'
        % Right Hemisphere - Medial View
        patch('Vertices', vertices_rh_pial, 'Faces', faces_rh, ...
            'FaceVertexCData', vertex_colors_rh, ...
            'FaceColor', 'interp', 'EdgeColor', 'none');
        lighting gouraud;
        camlight('headlight');
        view([90, 0]);
        lighting flat;
        camlight('headlight');
        axis equal off;
        xlim([-100, 100]); ylim([-110, 90]); zlim([-100, 100]);
        material dull;

    case 'p1'
        % Left Hemisphere - Lateral View
        patch('Vertices', vertices_lh_pial, 'Faces', faces_lh, ...
            'FaceVertexCData', vertex_colors_lh, ...
            'FaceColor', 'interp', 'EdgeColor', 'none');
        lighting gouraud;
        camlight('headlight');
        view([90, 0]);
        lighting flat;
        camlight('headlight');
        axis equal off;
        xlim([-100, 100]); ylim([-120, 80]); zlim([-100, 100]);
        material dull;

    case 'p3'
        % Right Hemisphere - Lateral View
        patch('Vertices', vertices_rh_pial, 'Faces', faces_rh, ...
            'FaceVertexCData', vertex_colors_rh, ...
            'FaceColor', 'interp', 'EdgeColor', 'none');
        lighting gouraud;
        camlight('headlight');
        view([-90, 0]);
        lighting flat;
        camlight('headlight');
        axis equal off;
        xlim([-100, 100]); ylim([-120, 80]); zlim([-100, 100]);
        material dull;
end

rmpath('./data');
end
