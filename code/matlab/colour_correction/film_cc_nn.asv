clear; close all;

%% Load precomputed patch data
inputFile = '/home/oem/eliza/data/xyz_lab_rgb/colorchecker/yoda_halogen_fuji_underexp_balanced_colorchecker.mat';  
refFile = '/home/oem/eliza/masters-thesis/data/CCpassport_spectra.txt';
ill = importdata('../../../data/CIE_D50.txt');
CMFs = importdata('../../../data/CIE2degCMFs_1931.txt');

S = load(inputFile);

XYZ_input     = S.patchXYZ;
Lab_input     = S.patchLab;
RGB_input     = S.patchRGB;
RGB_lin_input = S.patchRGB_lin;

%% Load spectral reflectance reference file
refData = readmatrix(refFile, 'NumHeaderLines', 1);  
ref_wl = refData(:,1);   % Wavelengths (1st column)
ref_spectra = refData(:,2:end);    % 24 columns of reflectance data

% Trim wavelengths to 380–780 nm
valid_idx = ref_wl >= 380 & ref_wl <= 780;
ref_wl = ref_wl(valid_idx);
ref_spectra = ref_spectra(valid_idx, :);

% Interpolate illuminant and CMFs to match wavelengths
illIP = interp1(ill(:,1), ill(:,2), ref_wl, 'spline');
CMFsIP = interp1(CMFs(:,1), CMFs(:,2:4), ref_wl, 'spline');  % x̄ ȳ z̄

% Compute reference XYZ
sp_tristREF = CMFsIP .* illIP;
k = 100 / sum(sp_tristREF(:,2));
patchXYZ_ref = k * (ref_spectra' * sp_tristREF);  % [24×3]

% Convert to Lab and RGB spaces
patchLab_ref      = xyz2lab_custom(patchXYZ_ref);
patchRGB_ref      = xyz2prophoto(patchXYZ_ref / 100, true);
patchRGB_lin_ref  = xyz2prophoto(patchXYZ_ref / 100, false);

% Set reference variables
XYZ_ref     = patchXYZ_ref;
Lab_ref     = patchLab_ref;
RGB_ref     = patchRGB_ref;
RGB_lin_ref = patchRGB_lin_ref;

%% Prepare data
N = size(XYZ_input, 1);
rng(42);  % Reproducibility
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

        % Train neural network
        % net = fitnet([10 10]);  % Two hidden layers
        % net.trainParam.showWindow = false;
        % net.trainParam.showCommandLine = false;
        % net = train(net, X_train, Y_train);

        % % Predict all
        % Y_pred_norm = net(X_all);
        % Y_pred = mapminmax('reverse', Y_pred_norm, psY)';

        % Improved Neural Network
        best_net = [];
        best_perf = inf;
    
        for trial = 1:5  % Multiple initializations
            net = fitnet([32 16], 'trainlm');  % Larger, deeper network
            net.performFcn = 'mse';  % Mean Squared Error
            net.trainParam.showWindow = false;
            net.trainParam.showCommandLine = false;
            
            % Regularization to prevent overfitting
            net.performParam.regularization = 0.01;
            
            % Early stopping
            net.divideParam.trainRatio = 0.7;
            net.divideParam.valRatio = 0.15;
            net.divideParam.testRatio = 0.15;
        
            [net_trained, tr] = train(net, X_train, Y_train);
            
            % Validate on test set (built-in)
            Y_val_pred = net_trained(X_train);
            perf = perform(net_trained, Y_train, Y_val_pred);
        
            if perf < best_perf
                best_perf = perf;
                best_net = net_trained;
            end
        end
        
        % Predict using best network
        Y_pred_norm = best_net(X_all);
        Y_pred = mapminmax('reverse', Y_pred_norm, psY)';


        
        % Evaluate in Lab space
        if strcmp(name, 'Lab')
            predLab = Y_pred;
            refLab = Lab_ref;
        
        elseif strcmp(name, 'RGB')
            % Convert predicted RGB to Lab via XYZ
            XYZ_pred = prophoto2xyz(Y_pred, true);
            predLab = xyz2lab_custom(XYZ_pred);
        
            refLab = Lab_ref;
        
        elseif strcmp(name, 'RGB_lin')
            XYZ_pred = prophoto2xyz(Y_pred, false);
            predLab = xyz2lab_custom(XYZ_pred);
        
            % Use spectral Lab ref
            refLab = Lab_ref;
        
        else  % XYZ
            predLab = xyz2lab_custom(Y_pred);
            refLab  = Lab_ref;
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
fprintf('\n---- Cross-Validation Summary (ΔE2000, Neural Network) ----\n');
for f = 1:numel(fields)
    name = fields{f};
    fprintf('%-8s | Mean = %.2f ± %.2f | Max = %.2f ± %.2f\n', ...
        name, ...
        mean(mean_err.(name)), std(mean_err.(name)), ...
        mean(max_err.(name)),  std(max_err.(name)));
end
