function heb_sim_run()
% HEB_SIM_RUN - Runs simulations across multiple SNR levels and generates 
%               figures for the main text and supporting information.
%  
% This function is associated with the study: Greaves et al. (2024).
% DOI: https://doi.org/10.1101/2024.04.03.587831
%
% Description:
%   This function executes the `heb_sim` function across a range of
%   signal-to-noise ratio (SNR) values. The goal is to simulate
%   ground-truth effective connectivity, and then perform parameter 
%   recovery via both a) the inversion of a hierarchical empirical Bayes 
%   model, and b) the fitting of a structurally informed multivariate 
%   autoregressive (MVAR) model. The function then calls `heb_sim_fig` to 
%   generate figures for both the main text and supporting information, 
%   consistent with the analyses presented in the associated publication.
%   Note that running this function as is creates a folder called 'output' 
%   in the currect working directory.%
%  
% Example Usage:  
%   heb_sim_run(); % Runs simulations for all SNR levels
%  
% See also: heb_sim, heb_sim_fig  

% Add the simulation directory to MATLAB's path  
addpath('./sim/');

% Define SNR levels to evaluate  
snr_levels = [1, 5, 10, 50, 100, 200];

% Loop through each SNR level and run the simulation  
for i = 1:length(snr_levels)
    fprintf('Running heb_sim for SNR = %d...\n', snr_levels(i));
    
    % Start timer for this simulation
    snr_tic = tic;
    
    % Run simulation
    heb_sim(snr_levels(i)); 
    
    % Display elapsed time
    elapsed_time = toc(snr_tic);
    fprintf('Completed SNR = %d in %.2f seconds.\n',...
        snr_levels(i), elapsed_time);
end

% Generate figures
fprintf('Generating figures...\n');
heb_sim_fig();

end
