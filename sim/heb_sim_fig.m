function heb_sim_fig()
% HEB_SIM_FIG - Generates the main text figure for structurally informed  
%               estimation of effective connectivity (in silico analyses).  
%  
% Syntax: heb_sim_fig()  
%  
% Description: 
%   This figure corresponds to the **main text results** of the  
%   associated publication. For additional analyses and supporting  
%   information figure, see `heb_sim_sup_fig`.  
%  
% Example:  
%   heb_sim_fig(); % Generates the main text figure  
%  
% See also: heb_sim, heb_sim_run, heb_sim_sup_fig  

% Load required data
file = dir(fullfile(pwd, 'output', sprintf('*snr%d_*.mat', 5)));
load(fullfile(file(1).folder, file(1).name), 'SC', 'Ag', 'As', 'W',...
    'HEB', 'EXP', 'DCM', 'n', 's', 'a', 'b', 'pv');

% Results
%--------------------------------------------------------------------------
% Results for (structurally informed) MVAR model
W_group = mean(W, 3);
mvar_group_pearsons_r = corr(Ag(:), W_group(:));
mvar_rmse_group = sqrt(mean((W_group(:) - Ag(:)).^2));

% Results for (structurally informed) effective connectivity
heb_Ag = full(HEB.Ep);

% Correct for log-normal scaling (convert to Hz)
A = heb_Ag;
A = reshape(A, sqrt(numel(A)), sqrt(numel(A)));
A(logical(eye(sqrt(numel(A))))) =...
    -0.5 * exp(A(logical(eye(sqrt(numel(A))))));
heb_Ag = A(:);

% Grop-level statistics
heb_group_pearsons_r = corr(Ag(:), heb_Ag(:));
heb_rmse_group = sqrt(mean((heb_Ag(:) - Ag(:)).^2));

% Subject-level Pearson's r
mvar_pearsons_r = zeros(s, 1);
heb_pearsons_r = zeros(s, 1);

% Subject-level RMSE
mvar_sub_rmse = zeros(s, 1);
heb_sub_rmse = zeros(s, 1);

% Subject-level polarity error (%)
mvar_sign_pe = zeros(s, 1);
heb_sign_pe = zeros(s, 1);

for i = 1:s
    % Ground truth for subject
    A_sub = As{i};

    % Correct for log-normal scaling (convert to Hz)
    A = heb_Ag;
    A(logical(eye(sqrt(numel(A))))) =...
        -0.5 * exp(A(logical(eye(sqrt(numel(A))))));
    heb_Ag = A;

    % MVAR Pearson's r, RMSE, and PE
    W_v = W(:, :, i);
    mvar_pearsons_r(i) = corr(A_sub(:), W_v(:));
    mvar_sub_rmse(i) = sqrt(mean((W_v(:) - A_sub(:)).^2));
    mvar_sign_pe(i) = (sum(sign(W_v(:)) ~= sign(A_sub(:))))/...
        n^2*100;

    % HEB Pearson's r, RMSE, and PE
    heb_As = full(DCM{i}.Ep);
    heb_pearsons_r(i) = corr(A_sub(:), heb_As(:));
    heb_sub_rmse(i) = sqrt(mean((heb_As(:) - A_sub(:)).^2));
    heb_sign_pe(i) = (sum(sign(heb_As(:)) ~= sign(A_sub(:))))/...
        n^2*100;
end

% Figure
%--------------------------------------------------------------------------

fig = figure;
set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
font_size = 20;
xy_ax = max(1, ceil(max(abs([W_group(:); heb_Ag(:)]))));

%--------------------------------------------------------------------------

% Subplot A: Normalized structural connectivity
ax_a = subplot(2, 4, 1);
set(ax_a, 'OuterPosition', [0.058 0.57 0.205 0.35]); % 0.04 (x)
imagesc(SC);
xticks(1:n)
yticks(1:n)
xlabel('Region', 'FontSize', font_size);
ylabel('Region', 'FontSize', font_size);
c = colorbar(ax_a, 'eastoutside');
colormap(ax_a, 'cool');
c.Limits = [0, 1]; 
c.Ticks = [0, 1];
c.FontSize = font_size;
c.Label.String = 'arb.';
c.Label.FontSize = font_size;
c.Label.Position = [1.25, 0.5, 0];
set(gca, 'FontSize', font_size);

% Add subplot label
xlims = xlim; 
ylims = ylim; 
text(xlims(1), ylims(2), 'A',...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
    'FontSize', font_size+1, 'FontWeight', 'bold', 'Units',...
    'normalized', 'Position', [0.06, 1.20]);
hold off;

%--------------------------------------------------------------------------

% Subplot B: Effective connectivity
ax_b = subplot(2, 4, 2);
set(ax_b, 'OuterPosition', [0.3 0.57 0.205 0.35]); % 0.28 (x)
imagesc(Ag);
xticks(1:n)
yticks(1:n)
xlabel('Target region', 'FontSize', font_size);
ylabel('Source region', 'FontSize', font_size);
c = colorbar(ax_b, 'eastoutside');
C = redwhiteblue(-1, 1, 100);
colormap(ax_b, C);
clim([-1 1]);
c.Ticks = round([-1, 1],2);
c.FontSize = font_size; 
c.Label.String = 'Hz';
c.Label.FontSize = font_size;
c.Label.Position = [1.25, 0, 0];
set(gca, 'FontSize', font_size);

% Add subplot label
text(xlims(1), ylims(2), 'B',...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
    'FontSize', font_size+1, 'FontWeight', 'bold', 'Units',...
    'normalized', 'Position', [0.06, 1.20]);
hold off;

%--------------------------------------------------------------------------

% Subplot C: True vs. Estimated Parameters (HEB)
ax = subplot(2, 4, [3, 4]);
set(ax, 'OuterPosition', [0.52 0.56 0.46 0.40]);
scatter(Ag(:), heb_Ag(:), 150, 'g', 'filled');
hold on;
plot([-xy_ax, xy_ax], [-xy_ax, xy_ax], 'k--', 'LineWidth', 1.2);
ylim([-xy_ax, xy_ax]);
xlim([-xy_ax, xy_ax]);
xlabel('True group-level effective connectivity (Hz)', 'FontSize',...
    font_size);
ylabel({'MAP group-level', 'effective connectivity (Hz)'},...
    'FontSize', font_size);

% Place error metrics in the upper right-hand corner
xlims = xlim; 
ylims = ylim; 
text(xlims(2), ylims(2), ['\it{r}\rm = ', num2str(...
    heb_group_pearsons_r, '%.2f'),...
    '; RMSE: ', num2str(heb_rmse_group, '%.3f')],...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
    'FontSize', font_size, 'Units', 'normalized', 'Position', [1, 1.10]);
set(gca, 'FontSize', font_size);
grid on;

% Add subplot label
text(xlims(1), ylims(2), 'C',...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
    'FontSize', font_size+1, 'FontWeight', 'bold', 'Units',...
    'normalized', 'Position', [0.06, 1.15]);
hold off;

%--------------------------------------------------------------------------

% Subplot D: Ground-truth variance transformation versus BMA transformation
ax = subplot(2, 4, 5);
set(ax, 'OuterPosition', [0.02 0.1100 0.23 0.35]);
x = 0:0.01:1;
heb_bma_se(x, EXP.Fs, EXP.params.alphas, EXP.params.betas,...
    EXP.winning.alpha_bma, EXP.winning.beta_bma, pv, [0, 1, 0])
hold on;
plot(x, (a + x .* b), 'k--', 'LineWidth', 1.2)
xticks([0, 1])
yticks(round([EXP.winning.alpha_bma,...
    EXP.winning.alpha_bma+EXP.winning.beta_bma],2));
xlabel({'Normalized structural',' connectivity'}, 'FontSize', font_size);
ylabel({'Variance of group-level',' effective connectivity'},...
    'FontSize', font_size);
set(gca, 'FontSize', font_size);

% Add subplot label
text(xlims(1), ylims(2), 'D',...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
    'FontSize', font_size+1, 'FontWeight', 'bold', 'Units',...
    'normalized', 'Position', [0.06, 1.20]);
hold off;

%--------------------------------------------------------------------------

% Subplot E: True vs. Estimated Parameters (MVAR)
ax = subplot(2, 4, 6);
set(ax, 'OuterPosition', [0.265 0.1100 0.26 0.35])
scatter(Ag(:), W_group(:), 150, 'b', 'filled');
hold on;
plot([-xy_ax, xy_ax], [-xy_ax, xy_ax], 'k--', 'LineWidth', 1.2);
ylim([-xy_ax, xy_ax]);
xlim([-xy_ax, xy_ax]);
xlabel({'True group-level',' effective connectivity (Hz)'},...
    'FontSize', font_size);
ylabel({'Mean directed functional',...
    'connectivity (unitless)'}, 'FontSize', font_size);

% Place error metrics in the upper right-hand corner
xlims = xlim; 
ylims = ylim; 
text(xlims(2), ylims(2), ['\it{r}\rm = ', num2str(...
    mvar_group_pearsons_r, '%.2f'),...
    '; RMSE: ', num2str(mvar_rmse_group, '%.3f')],...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
    'FontSize', font_size, 'Units', 'normalized', 'Position', [1, 1.12]);
set(gca, 'FontSize', font_size);
grid on;

% Add subplot label
text(xlims(1), ylims(2), 'E',...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
    'FontSize', font_size+1, 'FontWeight', 'bold', 'Units',...
    'normalized', 'Position', [0.06, 1.2]);
hold off;

%--------------------------------------------------------------------------

% Subplot F: Pearson's r at the First, Subject Level
ax = subplot(2, 4, [7, 8]);
set(ax, 'OuterPosition', [0.54 0.1 0.44 0.40]);

% Create grouped bar plot
h1 = bar([heb_pearsons_r, mvar_pearsons_r], 'grouped');

% Adjust GroupWidth to reduce space between grouped bars
set(h1, 'GroupWidth', 0.88);

% Adjust individual bar positions using XOffset
h1(1).XOffset = -0.07;  
h1(2).XOffset =  0.07;  

% Set colors
h1(1).FaceColor = 'g'; % Green (HEB)
h1(2).FaceColor = 'b'; % Blue (MVAR)
h1(1).EdgeAlpha = 0; 
h1(2).EdgeAlpha = 0;
xline(1.5:1:49.5, "LineStyle","--", "Color",'k')
hold off;

% Set limits and labels
ylim([0, 1]);
xticks([1, s]);
xlabel('Instantiation (simulated subject)', 'FontSize', font_size);
ylabel('Pearson''s \it{r}', 'FontSize', font_size);
set(gca, 'FontSize', font_size);

% Create the legend
lgd = legend(h1, 'HEB', 'MVAR', 'FontSize', font_size,...
    'Location', 'northoutside');

% Set legend entries side-by-side
lgd.NumColumns = 2;

% Adjust legend position if needed
lgd.Position = [lgd.Position(1) + 0.117, lgd.Position(2) + 0.035,...
    lgd.Position(3), lgd.Position(4)];

% Add subplot label
xlims = xlim; 
ylims = ylim; 
text(xlims(1), ylims(2), 'F',...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
    'FontSize', font_size+1, 'FontWeight', 'bold', 'Units',...
    'normalized', 'Position', [0.06, 1.15]);

%--------------------------------------------------------------------------

% Supporting Information Figure
heb_sim_sup_fig(n, s, heb_sub_rmse, mvar_sub_rmse,...
    heb_sign_pe, mvar_sign_pe);

end
