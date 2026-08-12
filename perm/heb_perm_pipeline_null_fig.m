function fig = heb_perm_pipeline_null_fig(data_file, out_dir, panel_b_kind)
% =========================================================================
% heb_perm_pipeline_null_fig.m
% =========================================================================
% Figure script for the pipeline-level permutation null under scrambled
% structural connectivity (Supporting Information Fig. S5).
%
% Context:
% This script visualizes network-wise outputs from the full-pipeline
% permutation analysis in which each network-specific structural
% connectivity matrix was scrambled, the hyperparameter grid/BMA procedure
% was rerun, and each hierarchical model was re-inverted under the resulting
% permuted-structural-connectivity-implied prior covariance.
%
% Loaded quantity:
%   ./plperm/perm_si_figure_data_B1000.mat
% containing a structure named perm_si with B=1000 permutation results for
% each of the 17 networks.
%
% What is shown:
% Panel A (17 subplots):
%   - Gray histograms = sample-size-scaled group-level log-Bayes factors
%     under full-pipeline scrambled-connectivity permutations.
%   - Red vertical line = empirical structurally informed value.
%   - Text annotation = empirical percentile rank within the permutation
%     distribution.
%
% Panel B (single scatter):
%   - y-axis = empirical prior-variance modulation.
%   - x-axis = prior-variance modulation implied by the modal permuted
%     alpha,beta mapping.
%   - Dashed line: equality (x = y).
%
% Optional panel_b_kind values:
%   'joint_param_mode'    PVMI at the joint modal permuted alpha/beta pair
%   'pvmi_mode'           mode of the permuted PVMI distribution
%   'marginal_param_mode' PVMI at separate marginal alpha and beta modes
%   'beta_joint_mode'     empirical beta versus joint modal permuted beta
%   'beta_mode'           empirical beta versus marginal modal permuted beta
%   'logbf_mode'          PVMI for the permutation nearest modal logBF
%   'median'              median of the permuted PVMI distribution
%
% Usage:
%   cd perm
%   heb_perm_pipeline_null_fig
%
% =========================================================================

if nargin < 1 || isempty(data_file)
    data_file = fullfile(pwd, 'plperm', 'perm_si_figure_data_B1000.mat');
end
if nargin < 2 || isempty(out_dir)
    out_dir = fullfile(pwd, 'plperm', 'figures');
end
if nargin < 3 || isempty(panel_b_kind)
    panel_b_kind = 'joint_param_mode';
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

S = load(data_file, 'perm_si');
perm_si = S.perm_si;
rows = perm_si.by_network;

font_size = 18;
title_size = 14;
line_width = 1.4;
red = [0.8 0.1 0.1];
grey = 0.75 * [1 1 1];
fig_visible = 'on';
if ~usejava('desktop')
    fig_visible = 'off';
end

fig = figure('Color', 'w', ...
    'Name', ['Pipeline permutation null: ', panel_b_kind], ...
    'Units', 'normalized', ...
    'Visible', fig_visible, ...
    'Position', [0 1/2 1 1/1.3]);

hist_ylim = histogram_ylim(rows);
ax = cell(numel(rows) + 1, 1);

for i = 1:numel(rows)
    ax{i} = subplot(3, 6, i);
    plot_network_histogram(ax{i}, rows(i), hist_ylim, grey, red, ...
        font_size, title_size, line_width);
end

ax{18} = subplot(3, 6, 18);
plot_panel_b(ax{18}, rows, perm_si.colors, font_size, ...
    line_width, panel_b_kind);

annotation('textbox', [0.13, 0.94, 0.05, 0.05], ...
    'String', '\bf A', ...
    'FontSize', 20, ...
    'LineStyle', 'none', ...
    'Interpreter', 'tex');

posB = get(ax{18}, 'Position');
annotation('textbox', [posB(1)-0.02, posB(2)+posB(4), 0.05, 0.05], ...
    'String', '\bf B', ...
    'FontSize', 20, ...
    'LineStyle', 'none', ...
    'Interpreter', 'tex');

annotation('textbox', [0.12 0.02 0.60 0.04], ...
    'String', 'Sample-size-scaled log-Bayes factor', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 20, ...
    'LineStyle', 'none', ...
    'FitBoxToText', 'off');

annotation('textbox', [0.24 0.24 0.25 0.58], ...
    'String', 'Count', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 20, ...
    'Rotation', 90, ...
    'LineStyle', 'none', ...
    'FitBoxToText', 'off');

suffix = clean_suffix(panel_b_kind);
fig_file = fullfile(out_dir, ['heb_perm_pipeline_null_fig_', suffix, '.fig']);
png_file = fullfile(out_dir, ['heb_perm_pipeline_null_fig_', suffix, '.png']);
pdf_file = fullfile(out_dir, ['heb_perm_pipeline_null_fig_', suffix, '.pdf']);

savefig(fig, fig_file);
exportgraphics(fig, png_file, 'Resolution', 300);
exportgraphics(fig, pdf_file, 'ContentType', 'vector');

fprintf('\nSaved:\n');
fprintf('  %s\n', fig_file);
fprintf('  %s\n', png_file);
fprintf('  %s\n', pdf_file);

end

function plot_network_histogram(ax, row, hist_ylim, grey, red, ...
    font_size, title_size, line_width)
x = row.perm_scaled_logBF;
x = x(isfinite(x));

if isempty(x)
    title(ax, row.network_name, 'FontSize', title_size, ...
        'FontWeight', 'normal');
    axis(ax, 'off');
    return
end

nbins = max(10, min(28, round(sqrt(numel(x)))));
histogram(ax, x, nbins, 'FaceColor', grey, 'FaceAlpha', 1, ...
    'EdgeColor', 'none');
hold(ax, 'on');
xline(ax, row.empirical_scaled_logBF, '-', 'Color', red, ...
    'LineWidth', line_width);

xl = [min([x(:); row.empirical_scaled_logBF]), ...
    max([x(:); row.empirical_scaled_logBF])];
pad = max(0.02, diff(xl) * 0.12);
if diff(xl) == 0
    xl = xl + [-0.05 0.05];
else
    xl = xl + [-pad pad];
end
xlim(ax, xl);
ylim(ax, hist_ylim);

title(ax, row.network_name, 'FontSize', title_size, ...
    'FontWeight', 'normal');
add_permutation_summary_text(ax, row, font_size, red);
box(ax, 'off');
set(ax, 'FontSize', font_size, 'TickDir', 'out');
xlabel(ax, '');
ylabel(ax, '');
end

function add_permutation_summary_text(ax, row, font_size, red)
if ~isfield(row, 'perm_scaled_logBF_percentile') || ...
        ~isfield(row, 'perm_scaled_logBF_bayes_p')
    return
end

pct = row.perm_scaled_logBF_percentile;
if ~isfinite(pct)
    return
end

txt = sprintf('pct. = %.1f', pct);
xl = xlim(ax);
offset = 0.025 * diff(xl);
yl = ylim(ax);
y = yl(1) + 0.92 * diff(yl);

if strcmp(row.net, 'SomMotA')
    x = row.empirical_scaled_logBF - offset;
    h_align = 'right';
else
    x = row.empirical_scaled_logBF + offset;
    h_align = 'left';
end

text(ax, x, y, txt, ...
    'HorizontalAlignment', h_align, ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 18, ...
    'Color', [0.15 0.15 0.15], ...
    'Interpreter', 'tex');
end

function plot_panel_b(ax, rows, colors, font_size, line_width, mode_kind)
[x, y, x_label_lines, y_label_lines] = panel_b_values(rows, mode_kind);

ok = isfinite(x) & isfinite(y);
lims = [min([x(ok); y(ok)]), max([x(ok); y(ok)])];
if isempty(lims) || any(~isfinite(lims))
    lims = [0 1];
end
pad = max(0.03, diff(lims) * 0.12);
if diff(lims) == 0
    lims = lims + [-0.05 0.05];
else
    lims = lims + [-pad pad];
end

hold(ax, 'on');
plot(ax, lims, lims, 'k--', 'LineWidth', line_width);
for i = 1:numel(rows)
    if ok(i)
        scatter(ax, x(i), y(i), 60, colors(i, :), 'filled', ...
            'MarkerEdgeColor', 'none');
    end
end
xlim(ax, lims);
ylim(ax, lims);
axis(ax, 'square');

xlabel(ax, x_label_lines, 'FontSize', font_size);
ylabel(ax, y_label_lines, 'FontSize', font_size);
box(ax, 'off');
set(ax, 'FontSize', font_size, 'TickDir', 'out');
end

function [x, y, x_label_lines, y_label_lines] = panel_b_values(rows, mode_kind)
switch lower(mode_kind)
    case {'beta_joint_mode', 'beta_joint', 'beta'}
        x = [rows.perm_beta_joint_mode]';
        y = [rows.empirical_beta]';
        x_label_lines = {'Modal permuted', '\beta'};
        y_label_lines = {'Empirical', '\beta'};
    case {'beta_mode', 'beta_marginal'}
        x = [rows.perm_beta_mode]';
        y = [rows.empirical_beta]';
        x_label_lines = {'Marginal modal', 'permuted \beta'};
        y_label_lines = {'Empirical', '\beta'};
    case {'joint_param_mode', 'joint', 'param_mode'}
        x = [rows.perm_pvmi_at_joint_param_mode]';
        y = [rows.empirical_pvmi]';
        x_label_lines = {'Permuted prior-variance', 'modulation at modal \alpha,\beta'};
        y_label_lines = {'Empirical prior-', 'variance modulation'};
    case {'marginal_param_mode', 'marginal'}
        x = [rows.perm_pvmi_at_marginal_param_mode]';
        y = [rows.empirical_pvmi]';
        x_label_lines = {'Permuted prior-variance', 'modulation at marginal modes'};
        y_label_lines = {'Empirical prior-', 'variance modulation'};
    case {'pvmi_mode', 'mode'}
        x = [rows.perm_pvmi_mode]';
        y = [rows.empirical_pvmi]';
        x_label_lines = {'Mode of permuted', 'prior-variance modulation'};
        y_label_lines = {'Empirical prior-', 'variance modulation'};
    case {'logbf_mode', 'evidence_mode'}
        x = [rows.perm_pvmi_at_logBF_mode]';
        y = [rows.empirical_pvmi]';
        x_label_lines = {'Permuted prior-variance', 'modulation at modal evidence'};
        y_label_lines = {'Empirical prior-', 'variance modulation'};
    case {'median', 'pvmi_median'}
        x = [rows.perm_pvmi_median]';
        y = [rows.empirical_pvmi]';
        x_label_lines = {'Median permuted', 'prior-variance modulation'};
        y_label_lines = {'Empirical prior-', 'variance modulation'};
    otherwise
        error('Unknown panel_b_kind: %s', mode_kind);
end
end

function yl = histogram_ylim(rows)
mx = 0;
for i = 1:numel(rows)
    x = rows(i).perm_scaled_logBF;
    x = x(isfinite(x));
    if isempty(x)
        continue
    end
    nbins = max(10, min(28, round(sqrt(numel(x)))));
    counts = histcounts(x, nbins);
    mx = max(mx, max(counts));
end
if mx <= 0
    yl = [0 1];
else
    yl = [0 ceil(mx * 1.08)];
end
end

function suffix = clean_suffix(s)
suffix = regexprep(lower(s), '[^a-z0-9]+', '_');
suffix = regexprep(suffix, '^_+|_+$', '');
end
