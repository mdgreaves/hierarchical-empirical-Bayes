function heb_bma_se(x_values, Fs, alphas, betas, alpha_mean,...
    beta_mean, max_var, rgb_code)
% HEB_BMA_SE - Computes and plots the Bayesian Model Averaged  
%              function with standard error envelope.  
%  
% Syntax: heb_bma_se(x_values, Fs, alphas, betas, alpha_mean,  
%                    beta_mean, max_var, rgb_code)  
%  
% Description:  
%   This function computes the Bayesian Model Averaged (BMA) function  
%   and its associated standard error envelope based on sampled  
%   parameters (alphas and betas) and model evidence (Fs). The function  
%   visualizes the expected function along with its confidence bounds.  
%  
% Inputs:  
%   x_values   - Vector of x-axis values  
%   Fs         - Cell array of (relative) log-evidence for different models  
%   alphas     - Vector of alpha parameter samples  
%   betas      - Vector of beta parameter samples  
%   alpha_mean - Mean value of alpha from Bayesian Model Averaging  
%   beta_mean  - Mean value of beta from Bayesian Model Averaging  
%   max_var    - Maximum variance threshold for valid parameter pairs  
%   rgb_code   - RGB color code for the plot  
%  
% Outputs:  
%   This function does not return a value but generates a plot with:  
%     - A shaded confidence envelope around the BMA function.  
%     - A solid line representing the BMA-estimated function.  
%  
% Example:  
%   heb_bma_se(0:0.01:1, Fs, alphas, betas, 0.5, 0.2, 1, [0, 1, 0]);  
%  

% Create ND grid of parameters governing data-to-variance mapping.
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
evidence = exp(Fs_matrix);
total_evidence = sum(evidence(:));
if total_evidence == 0
    error(['Total evidence is zero. ',...
        'This can cause division by zero in weights calculation.']);
end
weights = evidence / total_evidence;

% Calculate variances of alpha and beta
alpha_var = sum(weights .* (valid_alphas - alpha_mean).^2);
beta_var = sum(weights .* (valid_betas - beta_mean).^2);

% Initialize BMA function and variance
f_bma = alpha_mean + x_values .* beta_mean;
f_var = alpha_var + (x_values.^2) .* beta_var;

% Compute confidence envelope
confidence_upper = f_bma + sqrt(f_var); % removing 1.96 here for SE
confidence_lower = f_bma - sqrt(f_var);

% Plot the BMA function with confidence envelope
fill([x_values, fliplr(x_values)], [confidence_upper,...
    fliplr(confidence_lower)],...
    rgb_code, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;
plot(x_values, f_bma, 'Color', rgb_code, 'LineWidth', 3.5);

end