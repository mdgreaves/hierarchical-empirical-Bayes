# Hierarchical Empirical Bayes Model

This repository contains MATLAB functions for implementing the hierarchical empirical Bayes model of effective connectivity described in Greaves et al. (2024).

## Prerequisites

The functions in this repository rely on the statistical parametric mapping toolbox (SPM12), and—as this is an empirical Bayes procedure—require access to real data. To use these functions, you will need:
1. A cell array of > 2 DCMs inverted under identical prior assumptions that describe an effective connectivity network of *n* > 2 regions. Note that per the procedures reported in the associated study, it is assumed that the prior variance for all intraregional effective connections is fixed at 1/64, and that the prior variance for all interregional effective connections is fixed at 1/2 (see the *n* diagonal elements of `DCM.M.pC`). It is, however, easy to bypass this requirement—i.e., the assert commands—if one wishes to explore a different procedure.
2. A normalized [0,1] structural connectivity matrix (or matrix containing other relevant data) that has the same dimensions as the `A` (transition) matrices (such that the connection in `C(i,j)` corresponds to the connection in `DCM.Ep.A(i,j)`).

## Example Workflow

1. **Explore hierarchical empirical Bayes models in a test sample**:
   - Store the inverted (test) DCMs in a cell array `P`, and correctly-formatted secondary data—e.g., structural connectivity—into variable `C`.

   Example:
   ```matlab
   % Network name
   network = 'DMN';

   % Test directory
   test_dir = fullfile(pwd, 'test_subjects');

   % Load DCMs
   DCM_filelist = dir(fullfile(test_dir, '**', sprintf('*DCM_%s.mat', network)));

   % Load each DCM file into the cell array P
   P = cellfun(@(f) load(fullfile(f.folder, f.name), 'DCM').DCM, num2cell(DCM_filelist), 'UniformOutput', false);

   % Load structural connectivity data
   C = load(fullfile(test_dir, dir(fullfile(test_dir, '*SC.mat')).name)).SC;

   % Explore hierarchical empirical Bayes models
   heb_study(P, C, network);
   ```

2. **Assess the consistency of the Bayesian model average (BMA) data-to-variance mapping**:
   - Repeat the steps outlined above, loading the inverted (holdout) DCMs in a cell array `Pv`, and secondary data into variable `Cv`.
   - Store the path to the `HEB` file that was saved during the previous step.

   Example:
   ```matlab
   % Store the path to the 'HEB' file that was saved
   HEB = fullfile(pwd, sprintf('HEB_explore_%s.mat', network));

   % Validate Bayesian model average (BMA) data-to-variance mapping
   heb_study(Pv, Cv, network, HEB);
   ```

## Flexibility and Interpretation

The code provided in this repository can be easily modified to consider different data-to-prior-variance mappings, and thus consider different hypotheses regarding the relationship between structural and effective connectivity. The output provided by the `heb_study` function is easy to interpret for those with some knowledge of DCM (or similar models inverted using SPM's nonlinear systems identification function). It can be used to address a number of questions related to the coupling between structural and effective connectivity, for example.

## Contact

Questions? Please feel free to reach out.

## References

Greaves et al. (2024). Structurally informed resting-state effective connectivity recapitulates cortical hierarchy. DOI: [https://doi.org/10.1101/2024.04.03.587831](https://doi.org/10.1101/2024.04.03.587831)
