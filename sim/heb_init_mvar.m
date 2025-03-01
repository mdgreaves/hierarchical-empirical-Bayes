function [W, c] = heb_init_mvar(y, SC)
% heb_init_mvar - Estimates MVAR weights using structural connectivity.
%
% Syntax:  [W, c] = heb_init_mvar(timeSeries, SC)
%
% Inputs:
%   y          - Matrix (n x T) of time-series data for n regions and T
%                time points. 
%   SC         - Structural connectivity matrix (n x n), where SC(i,j) ≠ 0
%                indicates a connection from region j to region i.
%
% Outputs:
%   W          - MVAR weight matrix (n x n).
%   c          - Offset vector (n x 1).
%
% Description:
%   This function estimates the weights (W) and biases (c) of a
%   multivariate autoregressive (MVAR) model for a set of time-series data.
%   Structural connectivity (SC) is used to define potential connections.
%   Ordinary least squares regression is performed for each region using
%   time-shifted data from its connected neighbors as predictors.
%
% Example:
%   [W, c] = heb_init_mvar(timeSeries, SC);
%
% Adapted from:
% Tanner, J., Faskowitz, J., Teixeira, A.S. et al. A multi-modal,
% asymmetric, weighted, and signed description of anatomical connectivity.
% Nat Commun 15, 5865 (2024). https://doi.org/10.1038/s41467-024-50248-6

% Validate inputs
[n, T] = size(y);

% Initialize outputs
W = zeros(n, n); % Weight matrix
c = zeros(n, 1); % Offset terms

% Loop over each region
for i = 1:n
    % Get connected neighbors (excluding self-loops)
    neighbors = find(SC(i, :) ~= 0 & (1:n) ~= i);

    if isempty(neighbors)
        warning('Region %d has no connected neighbors.', i);
        continue;
    end

    % Prepare predictor matrix (time-shifted neighbors)
    predictors = y(neighbors, 1:T-1)';

    % Prepare target (current time point for region i)
    target = y(i, 2:T)';

    % Add a bias term to predictors
    X = [predictors, ones(size(predictors, 1), 1)];

    % Perform ordinary least squares regression
    coeffs = (X' * X) \ (X' * target);

    % Extract weights and bias
    W(i, neighbors) = coeffs(1:end-1);  % Assign weights to neighbors
    c(i) = coeffs(end);                 % Assign bias
end

end
