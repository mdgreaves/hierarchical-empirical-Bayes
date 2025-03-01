function heb_sim_sup_fig(n, s, heb_sub_rmse, mvar_sub_rmse,...
    heb_sign_pe, mvar_sign_pe)
% HEB_SIM_SUP_FIG - Generates the supporting information figure  
%                   for structurally informed effective connectivity.  
%  
% Syntax: heb_sim_sup_fig(n, s, heb_sub_rmse, mvar_sub_rmse,  
%                         heb_sign_pe, mvar_sign_pe)  
%  
% Description:  
%   This function generates additional analyses that supplement the main  
%   text figure created by `heb_sim_fig`. The supporting information  
%   figure evaluates model performance at both levels, comparing  
%   Multivariate Autoregressive (MVAR) modeling and Hierarchical  
%   Empirical Bayes (HEB) applied to Dynamic Causal Models (DCMs).  
%  
%   The figure includes:  
%     - RMSE across simulated subjects for HEB and MVAR (Subplot A).  
%     - Polarity error (% of incorrect sign assignments) across subjects  
%       for both methods (Subplot B).  
%     - Macro-F1 scores computed at different signal-to-noise ratios (SNR)  
%       for both methods (Subplot C).  
%  
% Input Arguments:  
%   n              - Number of brain regions  
%   s              - Number of simulated subjects  
%   heb_sub_rmse   - RMSE values for HEB at the subject level  
%   mvar_sub_rmse  - RMSE values for MVAR at the subject level  
%   heb_sign_pe    - Polarity error (%) for HEB at the subject level  
%   mvar_sign_pe   - Polarity error (%) for MVAR at the subject level  
%  
% Example:  
%   heb_sim_sup_fig(6, 50, heb_sub_rmse, mvar_sub_rmse, heb_sign_pe,  
%                   mvar_sign_pe);  
%  
% See also: heb_sim, heb_sim_fig  

% Subplot 6: Macro F-score
snr = [1, 5, 10, 50, 100, 200];
[heb_f_scores, mvar_f_scores] = deal(nan(1, length(snr)));

% Loop through each SNR value
for i = 1:length(snr)
    % Find the file for the current SNR
    file = dir(fullfile(pwd, 'output', sprintf('*_snr%s_*.mat',...
        num2str(snr(i)))));

    % Load required data
    load(fullfile(file(1).folder, file(1).name), ...
        'Ag', 'W', 'HEB');
    W_group = mean(W, 3);
    heb_Ag = full(HEB.Ep);
    heb_Ag = reshape(heb_Ag, n, n);

    % Compute F-scores
    heb_f_scores(i) = heb_fscore(Ag, heb_Ag);
    mvar_f_scores(i) = heb_fscore(Ag, W_group);
end

% Figure
%--------------------------------------------------------------------------

fig = figure;
set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
font_size = 20;

%--------------------------------------------------------------------------

% Subplot A: RMSE at the First, Subject Level
ax_1 = subplot(1,3,1);
set(ax_1, 'OuterPosition', [0.0005 0.5 0.35 .4]);

% Create grouped bar plot
h1 = bar([heb_sub_rmse, mvar_sub_rmse], 'grouped');

% Adjust GroupWidth to reduce space between grouped bars
set(h1, 'GroupWidth', 0.88); 

% Adjust individual bar positions using XOffset
h1(1).XOffset = -0.07;  
h1(2).XOffset =  0.07;  

% Set colors
h1(1).FaceColor = 'g'; 
h1(2).FaceColor = 'b'; 
h1(1).EdgeAlpha = 0; 
h1(2).EdgeAlpha = 0;
xline(1.5:1:49.5, "LineStyle","--", "Color",'k')

% Set limits and labels
ylim([0, 1]);
xticks([1, s]);
xlabel('Instantiation (simulated subject)', 'FontSize', font_size);
ylabel('RMSE', 'FontSize', font_size);
set(gca, 'FontSize', font_size);

% Create the legend
lgd = legend(h1, 'HEB', 'MVAR', 'FontSize', font_size, 'Location',...
    'northoutside');

% Set legend entries side-by-side
lgd.NumColumns = 2;

% Adjust legend position if needed
lgd.Position = [lgd.Position(1) + 0.091, lgd.Position(2) + 0.04,...
    lgd.Position(3), lgd.Position(4)];

% Add subplot label
xlims = xlim; 
ylims = ylim; 
text(xlims(1), ylims(2), 'A',...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
    'FontSize', font_size+1, 'FontWeight', 'bold', 'Units',...
    'normalized', 'Position', [0.06, 1.15]);
hold off;

%--------------------------------------------------------------------------

% Subplot C: Macro-F1
ax_3 = subplot(1,3,3);
set(ax_3, 'OuterPosition', [0.75, 0.5, 0.25, .4])
plot(1:numel(snr), heb_f_scores, '-o', 'LineWidth', 3.5, 'Color', 'g',...
    'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'HEB');
hold on;
plot(1:numel(snr), mvar_f_scores, '-o', 'LineWidth', 3.5, 'Color', 'b',...
    'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'MVAR');
hold off;

% Customize plot
xticks(1:numel(snr)); 
xticklabels(string(snr)); 
xlabel('SNR');
ylabel('Macro F1-score');
set(gca, 'FontSize', font_size);

% Create the legend and assign to ldg2
ldg2 = legend('Location', 'best', 'FontSize', font_size, 'Location',...
    'northoutside'); 

% Set legend entries side-by-side
ldg2.NumColumns = 2;

% Adjust legend position if needed
ldg2.Position = [ldg2.Position(1) + 0.028, ldg2.Position(2) + 0.051,...
    ldg2.Position(3), ldg2.Position(4)];

xlim([1 numel(snr)]);
grid on;

% Add subplot label
text(xlims(1), ylims(2), 'C',...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
    'FontSize', font_size+1, 'FontWeight', 'bold', 'Units',...
    'normalized', 'Position', [0.06, 1.17]);
hold off;

%--------------------------------------------------------------------------

% Subplot B: PE at the First, Subject Level
ax_2 = subplot(1,3,2);
set(ax_2, 'OuterPosition', [0.375 0.5 .35 .4]);

% Create grouped bar plot
h2 = bar([heb_sign_pe, mvar_sign_pe], 'grouped');

% Adjust GroupWidth to reduce space between grouped bars
set(h2, 'GroupWidth', 0.88); 

% Adjust individual bar positions using XOffset
h2(1).XOffset = -0.07; 
h2(2).XOffset =  0.07;  

% Set colors
h2(1).FaceColor = 'g'; 
h2(2).FaceColor = 'b'; 
h2(1).EdgeAlpha = 0; 
h2(2).EdgeAlpha = 0;
xline(1.5:1:49.5, "LineStyle","--", "Color",'k')

% Set limits and labels
ylim([0, 100]);
xticks([1, s]);
xlabel('Instantiation (simulated subject)', 'FontSize', font_size);
ylabel('Polarity error (%)', 'FontSize', font_size);
set(gca, 'FontSize', font_size);

% Create the legend
lgd = legend(h2, 'HEB', 'MVAR', 'FontSize', font_size, 'Location',...
    'northoutside');

% Set legend entries side-by-side
lgd.NumColumns = 2; % Number of columns for the legend entries

% Adjust legend position if needed
lgd.Position = [lgd.Position(1) + 0.091, lgd.Position(2) + 0.04,...
    lgd.Position(3), lgd.Position(4)];

% Add subplot label
xlims = xlim; 
ylims = ylim; 
text(xlims(1), ylims(2), 'B',...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
    'FontSize', font_size+1, 'FontWeight', 'bold', 'Units',...
    'normalized', 'Position', [0.06, 1.15]);
hold off;

end