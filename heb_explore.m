function HEB = heb_explore(P, C, params)
% HEB_EXPLORE Function to explore hierarchical empirical Bayes models.
%
% This function is associated with the study: Greaves et al. (2024). 
% DOI: https://doi.org/10.1101/2024.04.03.587831
%
% DESCRIPTION:
% This function explores the hierarchical empirical Bayes models by
% implementing Bayesian model reduction (BMR) to obtain updated posteriors
% and relative model evidence for different reduced models based on 
% data-to-prior-variance mappings.
%
% INPUTS:
%   P        -  Cell array of DCMs (which constitute the first level of the 
%               HEB model). 
%   C        -  Normalized [0,1] structural connectivity data used to 
%               inform third-level empirical priors. In principle, this 
%               could be other relevant data.
%   params   -  Structure containing the following fields:
%                .alphas: array of alpha values (intercept) to be explored.
%                .betas:  array of beta values (slope) to be explored.
%                .field:  string specifying the field(s) to be treated as 
%                         random effects.
%                .name:   string, name of the network used for output.
%
% OUTPUT:
%   HEB      -  Structure containing the results of the hierarchical 
%               empirical Bayes exploration, including posteriors, model 
%               structure, free energy estimates, and the Bayesian model 
%               average (BMA) parameters for the data-to-variance mappings.
%
% REQUIREMENTS:
%   - SPM12 must be installed and added to the MATLAB path.
%
% USAGE:
%   HEB = heb_explore(P, C, params)
%
% Example:
%   P = {...}; % Define subject-level DCM structures
%   C = ...;   % Define normalized structural connectivity data
%   params = struct('alphas', linspace(-0.5, 0.5, 30), ...
%                   'betas', linspace(0, 0.5, 30), ...
%                   'field', 'A', ...
%                   'name', 'default_network');
%   HEB = heb_explore(P, C, params);

% Store params in output structure.
HEB.params = params;

% Invert hierarchical empirical Bayes model without incorporating 
% third-level empirical priors (i.e., SPM's parametric empirical Bayes 
% [PEB] model).
[HEB_null, HEB_subs_null] = heb_capture(P, struct(), params.field);

% Store posteriors, model structure and free energy in output structure.
HEB.HEB_null.Ep = HEB_null.Ep; 
HEB.HEB_null.Cp = HEB_null.Cp; 
HEB.HEB_null.M = HEB_null.M; 
HEB.HEB_null.F = HEB_null.F; 
HEB.HEB_null.OutputVL = HEB_null.OutputVL;
HEB.HEB_null.F_subs_null = arrayfun(@(p1) HEB_subs_null{p1, 1}.F,...
    1:length(HEB_subs_null))';

% Create ND grid of parameters governing data-to-variance mapping.
[Alphas, Betas] = ndgrid(params.alphas, params.betas);

% Identify valid combinations
p = P{1}.M.pC(sub2ind(size(P{1}.M.pC), (1:P{1}.n^2)', (1:P{1}.n^2)'));
epsilon = 1e-5;
valid_indices = (Alphas >= epsilon) & (Alphas <= p(2)) & ...
                (Alphas + Betas >= epsilon) & (Alphas + Betas <= p(2));

% Extract valid combinations
valid_alphas = Alphas(valid_indices);
valid_betas  = Betas(valid_indices);

% Use Bayesian model reduction (BMR) to implement hierarchical empirical 
% Bayes model and obtain the updated posteriors and (relative) model 
% evidence associated with different reduced models resulting from 
% different data-to-prior-variance mappings.
fprintf(['Exploring HEB models with BMR procedure. Model: ',...
    sprintf('%s\n\n', params.name)]);
[HEB.Fs, HEB.sEs, HEB.sCs] = arrayfun(@(p1, p2)...
    heb_bmr(p1, p2, C, HEB_null), valid_alphas, valid_betas,...
    'UniformOutput', false);

% Find the maximum free energy (F) across all parameter regimes, and then 
% locate the corresponding posterior expectations and uncertainty.
Fs_matrix = cell2mat(HEB.Fs);
[HEB.winning.maxF, idx] = max(Fs_matrix(:));
[i, j, k] = ind2sub(size(HEB.Fs), idx);
HEB.winning.Ep = HEB.sEs{i, j, k};
HEB.winning.Cp = HEB.sCs{i, j, k};

% Identify the parameters that yield the maximum F.
HEB.winning.alpha = valid_alphas(i, j, k);
HEB.winning.beta = valid_betas(i, j, k);

% Check Fs_matrix values
if any(isinf(Fs_matrix(:))) || any(isnan(Fs_matrix(:)))
    error('Fs_matrix contains Inf or NaN values.');
end

% Identify parameters of the Bayesian model average (BMA) data-to-variance 
% mapping, by first calculating weights based on the normalized evidence.
evidence = exp(Fs_matrix);
total_evidence = sum(evidence(:));
if total_evidence == 0
    error(['Total evidence is zero. ',...
        'This can cause division by zero in weights calculation.']);
end
weights = evidence / total_evidence;

% Check for NaN values in weights
if any(isnan(weights(:)))
    error(['Weights contain NaN values. ',...
        'Check the computation of Fs_matrix and evidence.']);
end

% Ensure dimensions are consistent
if any(size(weights) ~= size(valid_alphas)) || any(size(weights)...
        ~= size(valid_betas))
    error('Mismatch in matrix dimensions.');
end

% Compute Bayesian model average parameters
HEB.winning.alpha_bma = sum(sum(weights .* valid_alphas));
HEB.winning.beta_bma = sum(sum(weights .* valid_betas));

end
