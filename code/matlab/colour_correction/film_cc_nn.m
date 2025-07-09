clear; close all;

%% Load input and reference patch data
inputFile     = '/home/oem/eliza/data/xyz_lab_rgb/colorchecker/yoda_halogen_fuji_underexp_balanced_colorchecker.mat';
referenceFile = '/home/oem/eliza/data/xyz_lab_rgb/reference/xrite_cc_reference_official.mat';

input     = load(inputFile);
reference = load(referenceFile);

XYZ_input     = input.patchXYZ;
Lab_input     = input.patchLab;
RGB_input     = input.patchRGB;
RGB_lin_input = input.patchRGB_lin;

XYZ_ref     = reference.XYZ_ref;
Lab_ref     = reference.Lab_ref;
RGB_ref     = reference.RGB_ref;
RGB_lin_ref = reference.RGB_lin_ref;

%% Prepare data
N = size(XYZ_input, 1);
rng(42);
K = 6;
cv = cvpartition(N, 'KFold', K);

fields = {'XYZ', 'Lab', 'RGB', 'RGB_lin'};
ref_struct = struct('XYZ', XYZ_ref, 'Lab', Lab_ref, 'RGB', RGB_ref, 'RGB_lin', RGB_lin_ref);
input_struct = struct('XYZ', XYZ_input, 'Lab', Lab_input, 'RGB', RGB_input, 'RGB_lin', RGB_lin_input);
mean_err = struct(); max_err = struct();

%% Neural network regression
for f = 1:numel(fields)
    name = fields{f};
    fprintf('\nEvaluating correction in %s space (Neural Network)\n', name);
    mean_de = zeros(K, 1);
    max_de  = zeros(K, 1);

    for k = 1:K
        train_idx = training(cv, k);
        test_idx  = test(cv, k);

        X_train_raw = input_struct.(name)(train_idx, :)';
        Y_train_raw = ref_struct.(name)(train_idx, :)';
        X_all_raw   = input_struct.(name)';

        % Normalize
        [X_train, psX] = mapminmax(X_train_raw);
        [Y_train, psY] = mapminmax(Y_train_raw);
        X_all = mapminmax('apply', X_all_raw, psX);

        %% Improved MLP: 3-64-32-16-3
        net = fitnet([64 32 16], 'trainlm');
        net.trainParam.showWindow = false;
        net.trainParam.showCommandLine = false;
        net.trainParam.epochs = 1000;
        net.performFcn = 'mse';
        net.trainParam.max_fail = 20;  % Early stopping
        net.divideParam.trainRatio = 0.8;
        net.divideParam.valRatio   = 0.15;
        net.divideParam.testRatio  = 0.05;
        net.performParam.regularization = 0.02;

        net = train(net, X_train, Y_train);

        % Predict
        Y_pred_norm = net(X_all);
        Y_pred = mapminmax('reverse', Y_pred_norm, psY)';

        %% Evaluate in Lab space
        if strcmp(name, 'Lab')
            predLab = Y_pred;
            refLab  = Lab_ref;

        elseif strcmp(name, 'RGB')
            XYZ_pred = prophoto2xyz(Y_pred, true);
            predLab  = xyz2lab_custom(XYZ_pred);
            refLab   = Lab_ref;

        elseif strcmp(name, 'RGB_lin')
            XYZ_pred = prophoto2xyz(Y_pred, false);
            predLab  = xyz2lab_custom(XYZ_pred);
            refLab   = Lab_ref;

        else  % XYZ
            predLab = xyz2lab_custom(Y_pred);
            refLab  = Lab_ref;
        end

        % Compute deltaE only on test samples
        deltaE = deltaE2000(predLab(test_idx, :), refLab(test_idx, :));
        mean_de(k) = mean(deltaE);
        max_de(k)  = max(deltaE);

        fprintf(' Fold %d | ΔE2000 Mean = %.2f | Max = %.2f\n', ...
            k, mean_de(k), max_de(k));
    end

    % Store results
    mean_err.(name) = mean_de;
    max_err.(name)  = max_de;
end

%% Summary
fprintf('\n---- Cross-Validation Summary (ΔE2000, Neural Net) ----\n');
for f = 1:numel(fields)
    name = fields{f};
    fprintf('%-8s | Mean = %.2f ± %.2f | Max = %.2f ± %.2f\n', ...
        name, ...
        mean(mean_err.(name)), std(mean_err.(name)), ...
        mean(max_err.(name)),  std(max_err.(name)));
end
