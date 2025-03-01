function [Ag, SC] = heb_connectivity(n, a, b, dn, vg, T, TR, tau)
% heb_connectivity - Generates group-level effective and structural
%                    connectivity.
%
% Syntax:  [Ag, SC] = heb_connectivity(n, a, b, dn, vg, T, TR))
%
% Inputs:
%   n      - Number of regions (size of the square matrix).
%   a      - Scaling parameter for the effective connectivity.
%   b      - Weight parameter for combining structural connectivity.
%   dn     - Density of (non-zero) structural connectivity.
%   vg     - Variance of self connections
%   T      - Number of time points (observations) for BOLD simulation.
%   TR     - Repetition time for BOLD simulation.
%
% Outputs:
%   Ag     - Group-level effective connectivity matrix (n x n).
%   SC     - Structural connectivity matrix (n x n), normalized.
%
% Description:
%   This function generates a group-level effective connectivity matrix
%   (Ag) and a corresponding structural connectivity matrix (SC). The
%   matrices are generated with user-defined sparsity, scaling, and maximum
%   connection strength. Stability is ensured by verifying all eigenvalues
%   have negative real parts.
%
% Example:
%   [Ag, SC] = heb_connectivity(6, 1/10, 1/3, 1/3, 1/64, 1e3, 2);
%
% See also: eig, sprandnsym

Ag = nan(n);
is_valid = false;

% Generate group-level Jacobian and validate with BOLD simulation
while ~is_valid
    % Temporarily suppress and monitor warnings
    warn_state = warning('off', 'MATLAB:singularMatrix');
    lastwarn(''); % Clear previous warnings

    try
        is_stable = false;
        while ~is_stable
            % % Initialize a symetric random sparse matrix
            while true
                C = sprandsym(n, dn);
                if any(~any(C, 1)) || any(~any(C, 2)) 
                    % Exit the loop if at least one region is unconnected
                    break; 
                end
            end

            % Obtain structural connectivity
            C = C .* ~eye(n);           % Adjust diagonal entries
            SC = full(abs(C)...         % Normalize to max absolute value
                / max(abs(C(:))));

            % Obtain diagonal effective connectivity (log-normal variance)
            D = log(-2 * (-1/2 + sqrt(vg) * randn(n, 1)));
            
            % Obtain off-diagonal effective connectivity (normal variance)
            J = sqrt((b * SC + a) .* ~eye(n)) .* randn(n);
            J(logical(eye(n))) = D;
            
            % Verify stability (all eigenvalues should have
            % negative real parts)
            eig_vals = eig(J);
            if ~any(real(eig_vals) >= 0)
                is_stable = true;
                fprintf('Stable Jacobian for group...\n');
                fprintf(['Assessing fitness for numerical ',...
                    'integration...\n']);
            end
        end

        % Validate the Jacobian by attempting to generate BOLD time
        % series data
        heb_bold(T, TR, J, 4, tau);

        % Check if warnings occurred
        [~, warn_id] = lastwarn;
        if isempty(warn_id)
            % If no warnings, Jacobian and time-series are valid
            Ag = J;
            is_valid = true;
        else
            % If a warning occurred, retry
            fprintf('Warning detected: %s\n', warn_id);
            fprintf('Retrying for group...\n');
        end

    catch ME
        % Handle errors during Jacobian validation
        fprintf('Error detected: %s\n', ME.message);
        fprintf('Retrying...\n');
    end

    % Restore warning state
    warning(warn_state);
end
end
