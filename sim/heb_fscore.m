function F_score = heb_fscore(ground_truth, inferred)
% HEB_FSCORE computes the macro-averaged F1-score for inferred connectivity  
%   
% This function compares inferred connectivity to a ground truth matrix and   
% computes the macro-averaged F1-score across three classes: positive,    
% negative, and zero. It first thresholds values to classify them, then   
% computes precision, recall, and F1-score for each class before taking   
% the mean. The function assumes that the matrices have identical sizes.   
%  
% Inputs:  
%   ground_truth - NxM matrix of ground truth values  
%   inferred     - NxM matrix of inferred values  
%  
% Output:  
%   F_score      - Macro-averaged F1-score across classes  

    % Validate input dimensions
    if ~isequal(size(ground_truth), size(inferred))
        error(['Ground truth and inferred matrices ',...
            'must have the same dimensions.']);
    end
    
    % Vectorize matrices for comparison
    gt = ground_truth(:);
    inf = inferred(:);
    
    % Define thresholds for classifying values
    epsilon = 1e-5;
    gt_classes = classify_values(gt, epsilon);
    inf_classes = classify_values(inf, epsilon);
    
    % Define class labels
    classes = [-1, 0, 1];
    F_scores = zeros(1, numel(classes));
    
    % Calculate F-scores for each class
    for i = 1:numel(classes)

        % Define binary vectors for the current class
        gt_class = (gt_classes == classes(i));
        inf_class = (inf_classes == classes(i));
        
        % Calculate True Positives, False Positives, False Negatives
        TP = sum(gt_class & inf_class);
        FP = sum(~gt_class & inf_class);
        FN = sum(gt_class & ~inf_class);
        
        % Calculate precision and recall
        precision = TP / (TP + FP + eps);
        recall = TP / (TP + FN + eps);
        
        % Calculate F-score for the class
        F_scores(i) = 2 * (precision * recall) /...
            (precision + recall + eps);
    end
    
    % Compute macro-averaged F-score (average across classes)
    F_score = mean(F_scores, 'omitnan');
end

function class = classify_values(values, epsilon)
    % Classify values into 1 (positive), -1 (negative), or 0 (zero)
    class = zeros(size(values)); 
    class(values > epsilon) = 1; 
    class(values < -epsilon) = -1;
end
