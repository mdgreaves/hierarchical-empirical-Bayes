function [PEB, DCM] = heb_capture(P, M, field)
% HEB_CAPTURE Function to capture output from spm_dcm_peb and extract 
% specific information.
%
% This function is associated with the study: Greaves et al. (2024). 
% DOI: https://doi.org/10.1101/2024.04.03.587831
%
% DESCRIPTION:
% This function inverts SPM's parametric empirical Bayes (PEB) and captures 
% the output, extracting specific information related to the variational 
% Laplace (VL) iterations. It parses the output to obtain the iteration 
% number, free energy (F), and change in free energy (dF).
%
% INPUTS:
%   P       -   Cell array of DCMs (which constitute the first level of the 
%               PEB model). 
%   M       -   Structure containing model specifications.
%   field   -   String specifying the field(s) to be treated as random 
%               effects.
%
% OUTPUTS:
%   PEB     -   Structure containing the PEB results, including extracted 
%               iteration data.
%   DCM     -   Cell array of (reduced) DCM structures after implementing 
%               the PEB analysis.
%
% REQUIREMENTS:
%   - SPM12 must be installed and added to the MATLAB path.
%
% USAGE:
%   [PEB, DCM] = heb_capture(P, M, field)
%
% Example:
%   P = {...}; % Define subject-level DCM structures
%   M = struct(); % Define model specifications
%   field = 'A'; % Specify the field
%   [PEB, DCM] = heb_capture(P, M, field);

% Create an anonymous function that calls spm_dcm_peb with the inputs
fun = @() spm_dcm_peb(P, M, field); varName = getVarNames(fun);

% Capture the output from the function call
output_str = evalc(sprintf('[PEB, P] = %s();', varName{1}));
disp(output_str);

% Regular expression pattern to extract data
pattern = 'VL Iteration (\d+)\s+:\s+F = (-?\d+\.\d+) dF: (-?\d+\.\d+)\s+\[([-+]\d+\.\d+)\]';

% Using regexp to match the pattern
matches = regexp(output_str, pattern, 'tokens');

% Initialize result
result = struct('Iteration', {}, 'F', {}, 'dF', {}, 'Value', {});

% Parse the data into the result structure
for i = 1:length(matches)
    result(i).Iteration = str2double(matches{i}{1});
    result(i).F = str2double(matches{i}{2});
    result(i).dF = str2double(matches{i}{3});
    result(i).Value = str2double(matches{i}{4});
end

if exist('PEB', 'var') && exist('P', 'var')
    PEB.OutputVL = result;
    DCM = P;
else
    [PEB, DCM] = deal([]);
end

end

function varNames = getVarNames(varargin)
nVars = nargin;
varNames = cell(1, nVars);
for i = 1:nVars
    varNames{i} = inputname(i);
end
end
