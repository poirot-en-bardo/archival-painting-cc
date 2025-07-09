clear; close all;

%% Load patch data
inputFile = '/home/oem/eliza/data/xyz_lab_rgb/colorchecker/yoda_halogen_fuji_underexp_balanced_colorchecker.mat';  
refFile = '/home/oem/eliza/masters-thesis/data/CCpassport_spectra.txt';
ill = importdata('../../../data/CIE_D50.txt');
CMFs = importdata('../../../data/CIE2degCMFs_1931.txt');

S = load(inputFile);
XYZ_input     = S.patchXYZ;
Lab_input     = S.patchLab;
RGB_input     = S.patchRGB;
RGB_lin_input = S.patchRGB_lin;

%% Load reference spectral data
refData = readmatrix(refFile, 'NumHeaderLines', 1);  
ref_wl = refData(:,1);
ref_spectra = refData(:,2:end);

% Limit to 380–780 nm
valid_idx = ref_wl >= 380 & ref_wl <= 780;
ref_wl = ref_wl(valid_idx);
ref_spectra = ref_spectra(valid_idx, :);

illIP = interp1(ill(:,1), ill(:,2), ref_wl, 'spline');
CMFsIP = interp1(CMFs(:,1), CMFs(:,2:4), ref_wl, 'spline');
sp_tristREF = CMFsIP .* illIP;
k = 100 / sum(sp_tristREF(:,2));
patchXYZ_ref = k * (ref_spectra' * sp_tristREF);

patchLab_ref      = xyz2lab_custom(patchXYZ_ref);
patchRGB_ref      = xyz2prophoto(patchXYZ_ref / 100, true);
patchRGB_lin_ref  = xyz2prophoto(patchXYZ_ref / 100, false);

XYZ_ref     = patchXYZ_ref;
Lab_ref     = patchLab_ref;
RGB_ref     = patchRGB_ref;
RGB_lin_ref = patchRGB_lin_ref;

%% Setup
N = size(XYZ_input, 1);
rng(42);
K = 6;
cv = cvpartition(N, 'KFold', K);

fields = {'XYZ', 'Lab', 'RGB', 'RGB_lin'};
ref_struct = struct('XYZ', XYZ_ref, 'Lab', Lab_ref, 'RGB', RGB_ref, 'RGB_lin', RGB_lin_ref);
input_struct = struct('XYZ', XYZ_input, 'Lab', Lab_input, 'RGB', RGB_input, 'RGB_lin', RGB_lin_input);

mean_err = struct(); max_err = struct();

%% Neural network training
for f = 1:numel(fields)
    name = fields{f};
    fprintf('\nEvaluating in %s space (Improved NN)\n', name);
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

        % Ensemble training
        best_net = [];
        best_perf = inf;

        for trial = 1:10
            net = fitnet([64 32], 'trainlm');  % Bayesian regularization
            net.trainParam.showWindow = false;
            net.trainParam.showCommandLine = false;

            % Early stopping via internal split
            net.divideParam.trainRatio = 0.7;
            net.divideParam.valRatio = 0.2;
            net.divideParam.testRatio = 0.1;

            [net_trained, tr] = train(net, X_train, Y_train);
            Y_val = net_trained(X_train);
            perf = perform(net_trained, Y_train, Y_val);

            if perf < best_perf
                best_perf = perf;
                best_net = net_trained;
            end
        end

        % Predict and reverse normalization
        Y_pred_norm = best_net(X_all);
        Y_pred = mapminmax('reverse', Y_pred_norm, psY)';

        % Convert to Lab for evaluation
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
            predLab = xyz2lab_custom(Y_pred);
            refLab  = xyz2lab_custom(XYZ_ref);
        end

        deltaE = deltaE2000(predLab, refLab);

        mean_de(k) = mean(deltaE(test_idx));
        max_de(k)  = max(deltaE(test_idx));
        fprintf(' Fold %d | ΔE2000 Mean = %.2f | Max = %.2f\n', k, mean_de(k), max_de(k));
    end

    mean_err.(name) = mean_de;
    max_err.(name)  = max_de;
end

%% Summary
fprintf('\n---- Improved NN Summary (ΔE2000) ----\n');
for f = 1:numel(fields)
    name = fields{f};
    fprintf('%-8s | Mean = %.2f ± %.2f | Max = %.2f ± %.2f\n', ...
        name, ...
        mean(mean_err.(name)), std(mean_err.(name)), ...
        mean(max_err.(name)),  std(max_err.(name)));
end
