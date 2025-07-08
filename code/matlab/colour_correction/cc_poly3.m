clc; clear; close all;

%% Load precomputed patch data
inputFile = '/home/oem/eliza/data/xyz_lab_rgb/colorchecker/yoda_halogen_fuji_underexp_balanced_colorchecker.mat';  % <-- update this
S = load(inputFile);

XYZ_input     = S.patchXYZ;
Lab_input     = S.patchLab;
RGB_input     = S.patchRGB;
RGB_lin_input = S.patchRGB_lin;

% Assume same values are targets, for synthetic testing (normally use measured reference)
XYZ_ref     = XYZ_input;
Lab_ref     = Lab_input;
RGB_ref     = RGB_input;
RGB_lin_ref = RGB_lin_input;

% Reshape to N×3
N = size(XYZ_input, 1);
m = 6; n = 4;  % for plotting

rng(42);  % for reproducibility
K = 6;
cv = cvpartition(N, 'KFold', K);

% Storage
mean_err = struct();
max_err  = struct();

fields = {'XYZ', 'Lab', 'RGB', 'RGB_lin'};
ref_struct = struct('XYZ', XYZ_ref, 'Lab', Lab_ref, 'RGB', RGB_ref, 'RGB_lin', RGB_lin_ref);
input_struct = struct('XYZ', XYZ_input, 'Lab', Lab_input, 'RGB', RGB_input, 'RGB_lin', RGB_lin_input);

% Loop over each representation
for f = 1:numel(fields)
    name = fields{f};
    fprintf('\nEvaluating correction in %s space\n', name);
    mean_de = zeros(K, 1);
    max_de  = zeros(K, 1);

    for k = 1:K
        train_idx = training(cv, k);
        test_idx  = test(cv, k);

        X_train = poly3_features(input_struct.(name)(train_idx, :));
        Y_train = ref_struct.(name)(train_idx, :);

        % Linear regression
        coeffs = pinv(X_train) * Y_train;

        % Apply to all data
        X_all = poly3_features(input_struct.(name));
        Y_pred = X_all * coeffs;

        % Evaluate in Lab space
        if strcmp(name, 'Lab')
            predLab = Y_pred;
            refLab = Lab_ref;
        elseif strcmp(name, 'RGB')
            predLab = rgb2lab(Y_pred, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');
            refLab = rgb2lab(RGB_ref, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');
        elseif strcmp(name, 'RGB_lin')
            predLab = rgb2lab(Y_pred, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');
            refLab = rgb2lab(RGB_lin_ref, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');
        else
            predLab = xyz2lab(Y_pred);
            refLab  = xyz2lab(XYZ_ref);
        end

        deltaE = deltaE2000(predLab, refLab);

        mean_de(k) = mean(deltaE(test_idx));
        max_de(k)  = max(deltaE(test_idx));
        fprintf(' Fold %d | ΔE2000 Mean = %.2f | Max = %.2f\n', k, mean_de(k), max_de(k));
    end

    % Store results
    mean_err.(name) = mean_de;
    max_err.(name)  = max_de;
end

%% Summary
fprintf('\n---- Cross-Validation Summary (ΔE2000) ----\n');
for f = 1:numel(fields)
    name = fields{f};
    fprintf('%-8s | Mean = %.2f ± %.2f | Max = %.2f ± %.2f\n', ...
        name, ...
        mean(mean_err.(name)), std(mean_err.(name)), ...
        mean(max_err.(name)),  std(max_err.(name)));
end

%% --- 3rd-degree polynomial features
function X_poly = poly3_features(input_data)
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);
    X_poly = [ ...
        ones(size(a)), ...
        a, b, c, ...
        a.^2, b.^2, c.^2, ...
        a.*b, a.*c, b.*c, ...
        a.^3, b.^3, c.^3, ...
        a.^2.*b, a.^2.*c, ...
        b.^2.*a, b.^2.*c, ...
        c.^2.*a, c.^2.*b, ...
        a.*b.*c];
end
