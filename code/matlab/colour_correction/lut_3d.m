clc; clear; close all;

%% Set random seed for reproducibility
rng(10);

%% Load spectral data
disp('Loading spectral data...');
ill = importdata('../../data/CIE_D65.txt'); % Illuminant
CMFs = importdata('../../data/CIE2degCMFs_1931.txt'); % Color matching functions

%% Load the reference and test cubes
cubeFile = "../../data/colorChecker_SG/cubes/cubeCC_fuji-frame4.hdr";
refFile = "../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd_in] = size(inCUBE);
lincube = reshape(inCUBE,[],bd_in);

hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
lincube_ref = reshape(refCUBE,[],size(refCUBE,3));

%% Convert Spectral Data to RGB
disp('Processing RGB values...');
rgb_input = process_rgb(lincube, bands, ill, CMFs);
rgb_ref = process_rgb(lincube_ref, bands_ref, ill, CMFs);

%% Normalize Data
disp('Normalizing data...');
X_min = min(rgb_input, [], 1);
X_max = max(rgb_input, [], 1);
Y_min = min(rgb_ref, [], 1);
Y_max = max(rgb_ref, [], 1);

rgb_input_norm = (rgb_input - X_min) ./ (X_max - X_min);
rgb_ref_norm = (rgb_ref - Y_min) ./ (Y_max - Y_min);

%% Train-Test Split
rng(10);
num_samples = size(rgb_input_norm, 1);
train_ratio = 0.8;
train_idx = randperm(num_samples, round(train_ratio * num_samples));
test_idx = setdiff(1:num_samples, train_idx);

X_train = rgb_input_norm(train_idx, :);
Y_train = rgb_ref_norm(train_idx, :);
X_test = rgb_input_norm(test_idx, :);
Y_test = rgb_ref_norm(test_idx, :);

%% Generate a Regular 3D Grid for the LUT
disp('Generating True 3D LUT...');

% Define LUT grid resolution
LUT_size = 17; % 17x17x17 grid
[xq, yq, zq] = ndgrid(linspace(0, 1, LUT_size), ...
                       linspace(0, 1, LUT_size), ...
                       linspace(0, 1, LUT_size));

% Use scattered interpolants with nearest extrapolation
LUT_r = scatteredInterpolant(X_train(:,1), X_train(:,2), X_train(:,3), Y_train(:,1), 'linear', 'nearest');
LUT_g = scatteredInterpolant(X_train(:,1), X_train(:,2), X_train(:,3), Y_train(:,2), 'linear', 'nearest');
LUT_b = scatteredInterpolant(X_train(:,1), X_train(:,2), X_train(:,3), Y_train(:,3), 'linear', 'nearest');

% Compute LUT values
LUT_x = LUT_r(xq, yq, zq);
LUT_y = LUT_g(xq, yq, zq);
LUT_z = LUT_b(xq, yq, zq);

% Store as gridded interpolant
LUT_interp = griddedInterpolant(xq, yq, zq, cat(4, LUT_x, LUT_y, LUT_z), 'linear', 'nearest');

%% Apply LUT to Test Data
disp('Applying LUT to test set...');
Y_pred_test = LUT_interp(X_test(:,1), X_test(:,2), X_test(:,3));

% Ensure proper shape for rgb2lab()
Y_pred_test = squeeze(Y_pred_test);
if size(Y_pred_test,2) ~= 3
    Y_pred_test = reshape(Y_pred_test, [], 3);
end

%% Compute ΔE2000 for Test Data
disp('Evaluating test performance...');
lab_test_ref = rgb2lab(Y_test);
lab_test_pred = rgb2lab(Y_pred_test);
test_error_map = deltaE2000(lab_test_pred, lab_test_ref);

mean_test_error = mean(test_error_map(:));
max_test_error = max(test_error_map(:));
fprintf('Mean ΔE2000 Error on Test Set: %.4f\n', mean_test_error);
fprintf('Max ΔE2000 Error on Test Set: %.4f\n', max_test_error);

%% Apply LUT to Full Image
disp('Applying LUT to full image...');
corrected_rgb_norm = LUT_interp(rgb_input_norm(:,1), rgb_input_norm(:,2), rgb_input_norm(:,3));

% Ensure the output is in the correct format
corrected_rgb_norm = squeeze(corrected_rgb_norm); % Remove extra dimensions if any
if size(corrected_rgb_norm,2) ~= 3
    corrected_rgb_norm = reshape(corrected_rgb_norm, [], 3);
end

corrected_rgb = corrected_rgb_norm .* (Y_max - Y_min) + Y_min;

lab_ref = rgb2lab(rgb_ref);
lab_corr = rgb2lab(corrected_rgb);
error_map = deltaE2000(lab_corr, lab_ref);

mean_full_error = mean(error_map(:));
max_full_error = max(error_map(:));
fprintf('Mean ΔE2000 Error for Full Image: %.4f\n', mean_full_error);
fprintf('Max ΔE2000 Error for Full Image: %.4f\n', max_full_error);

%% Display Error Map
figure;
imagesc(reshape(error_map, m, n));
colormap jet;
colorbar;
clim([0 10]); 
title('\DeltaE2000 Error Map');
xlabel('Width');
ylabel('Height');

%% Function to Process RGB Data
function rgb = process_rgb(lincube, bands, ill, CMFs)
    illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
    CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
    sp_tristREF = CMFsIP .* illIP;
    trist = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

    % Convert XYZ to RGB (Adobe RGB)
    rgb = xyz2rgb(trist, 'ColorSpace', 'adobe-rgb-1998');
end
