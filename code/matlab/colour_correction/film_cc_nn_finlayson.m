clear; clc; close all;

%% Load data
inputFile     = '/home/oem/eliza/data/xyz_lab_rgb/colorchecker/yoda_halogen_fuji_underexp_balanced_colorchecker.mat';
referenceFile = '/home/oem/eliza/data/xyz_lab_rgb/reference/xrite_cc_reference_official.mat';

S_input       = load(inputFile);
S_ref         = load(referenceFile);

% Use gamma-encoded RGB as input and reference
RGB_input = S_input.patchRGB;
RGB_ref   = S_ref.RGB_ref;

N = size(RGB_input, 1);
rng(42);  % for reproducibility
K = 6;
cv = cvpartition(N, 'KFold', K);

mean_de = zeros(K, 1);
max_de  = zeros(K, 1);

for fold = 1:K
    fprintf('Fold %d/%d\n', fold, K);

    % Train/test split
    train_idx = training(cv, fold);
    test_idx  = test(cv, fold);

    % Separate luminance and chromaticity
    sum_input = sum(RGB_input, 2);
    chrom_input = RGB_input ./ max(sum_input, 1e-6);  % normalize
    chrom_input = chrom_input(:,1:2);  % drop 3rd for redundancy
    lum_input   = sum_input;

    sum_ref = sum(RGB_ref, 2);
    chrom_ref = RGB_ref ./ max(sum_ref, 1e-6);
    chrom_ref = chrom_ref(:,1:2);  % same drop
    lum_ref   = sum_ref;

    %% Train chromaticity net (2 → 2)
    layers_chrom = [
        featureInputLayer(2)
        fullyConnectedLayer(16)
        reluLayer
        fullyConnectedLayer(16)
        reluLayer
        fullyConnectedLayer(2)  % Output 2D chromaticity
        regressionLayer];

    opts = trainingOptions('adam', ...
        'MaxEpochs', 300, ...
        'MiniBatchSize', 8, ...
        'Verbose', false, ...
        'Plots', 'none');

    chrom_net = trainNetwork(chrom_input(train_idx, :), chrom_ref(train_idx, :), layers_chrom, opts);

    %% Train luminance net (1 → 1)
    layers_lum = [
        featureInputLayer(1)
        fullyConnectedLayer(1)  % Linear luminance mapping
        regressionLayer];

    lum_net = trainNetwork(lum_input(train_idx), lum_ref(train_idx), layers_lum, opts);

    %% Predict full set
    chrom_pred = predict(chrom_net, chrom_input);
    lum_pred   = predict(lum_net, lum_input);

    % Reconstruct 3-channel chromaticity
    chrom_pred = chrom_pred ./ max(sum(chrom_pred, 2), 1e-6);  % re-normalize
    chrom_pred_full = [chrom_pred, 1 - sum(chrom_pred, 2)];

    % Reconstruct RGB
    RGB_pred = chrom_pred_full .* lum_pred;
    RGB_pred = max(RGB_pred, 0);  % Clamp negatives

    %% Evaluate error in Lab space
    XYZ_pred = prophoto2xyz(RGB_pred, true);
    Lab_pred = xyz2lab_custom(XYZ_pred);

    XYZ_ref = prophoto2xyz(RGB_ref, true);
    Lab_ref = xyz2lab_custom(XYZ_ref);

    deltaE = deltaE2000(Lab_pred, Lab_ref);

    mean_de(fold) = mean(deltaE(test_idx));
    max_de(fold)  = max(deltaE(test_idx));
end

%% Final Results
fprintf('\n--- Exposure-Invariant Neural Network (RGB → RGB) ---\n');
fprintf('ΔE2000 Mean = %.2f ± %.2f | Max = %.2f ± %.2f\n', ...
    mean(mean_de), std(mean_de), mean(max_de), std(max_de));
