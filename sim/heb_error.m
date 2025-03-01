function [Y, As] = heb_error(Ag, s, ve, vi, T, TR, snr, tau)
% heb_error - Generates subject-level Jacobians and BOLD time series data.
%
% Syntax:  [Y, As] = heb_error(Ag, s, ve, vi, T, TR, snr)
%
% Inputs:
%   Ag  - Group-level Jacobian (n x n) representing effective connectivity.
%   s   - Number of subjects for which to generate data.
%   ve  - Variance of Gaussian noise added to Jacobian off-diagonal.
%   vi  - Variance of Gaussian noise added to Jacobian diagonal.
%   T   - Number of time points (observations) for BOLD simulation.
%   TR  - Repetition time for BOLD simulation.
%   snr    - Desired signal-to-noise ratio.
%
% Outputs:
%   Y   - Cell array (s x 1) of BOLD time-series data for each subject.
%   As  - Cell array (s x 1) of subject-level Jacobians.
%
% Description:
%   This function generates subject-specific Jacobians by adding Gaussian
%   noise to a group-level Jacobian and ensuring their stability. It
%   validates the Jacobians by simulating BOLD time series for each
%   subject, ensuring no warnings or errors occur during the process.
%
% Example:
%   [Y, As] = heb_error(Ag, 10, 1/64, 1/256, 1000, 2, 4);
%
% See also: heb_bold, eig

[n, m] = size(Ag);
assert(n == m, 'Group-level Jacobian must be square.');

% Initialize cell array to store subject-level Jacobians and
% time-series data
As = cell(s, 1);
Y = cell(s, 1);

for i = 1:s
    is_valid = false;

    % Generate subject-level Jacobian and validate with BOLD simulation
    while ~is_valid
        % Temporarily suppress and monitor warnings
        warn_state = warning('off', 'MATLAB:singularMatrix');
        lastwarn(''); % Clear previous warnings

        try
            is_stable = false;
            while ~is_stable
                % Obtain diagonal effective connectivity 
                % (log-normal variance)
                D = log(-2 * (-1/2 + sqrt(vi) * randn(n, 1)));

                % Obtain off-diagonal effective connectivity 
                % (normal variance)
                J = Ag + sqrt(ve) * ~eye(n) .* rand(n);
                J(logical(eye(n))) = D;

                % Verify stability (all eigenvalues should have
                % negative real parts)
                eig_vals = eig(full(J));
                if ~any(real(eig_vals) >= 0)
                    is_stable = true;
                    fprintf('Stable Jacobian for subject %s...\n',...
                        num2str(i));
                    fprintf(['Assessing fitness for numerical ',...
                        'integration...\n']);
                end
            end

            % Validate the Jacobian by attempting to generate BOLD time
            % series data
            y = heb_bold(T, TR, J, snr, tau);

            % Check if warnings occurred
            [~, warn_id] = lastwarn;
            if isempty(warn_id)
                % If no warnings, Jacobian and time-series are valid
                As{i} = J;
                Y{i} = y;
                is_valid = true;
            else
                % If a warning occurred, retry
                fprintf('Warning detected: %s\n', warn_id);
                fprintf('Retrying for subject %s...\n', num2str(i));
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
end
