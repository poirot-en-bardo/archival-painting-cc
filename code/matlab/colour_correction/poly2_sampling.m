%% Define parameters
clc; clear; close all;
roof = double(intmax('uint16')); % Histogram normalization factor

% Load spectral data
ill = importdata('../..//../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../../../data/CIE2degCMFs_1931.txt');
CMFs = CMFs_1931;
rng(10);

%% Select the reference and the cube to be corrected
cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_fuji-frame4.hdr";
refFile = "../../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

%% Load and process the cubes
hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd] = size(inCUBE);
lincube = reshape(inCUBE, [], bd);

hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
% Assuming that the spatial dimensions match for input and reference:
[m, n, bd] = size(refCUBE);
lincube_ref = reshape(refCUBE, [], bd);

% Interpolate illuminant and CMFs for the input cube
illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
sp_tristREF = CMFsIP .* illIP;

% Compute XYZ values for the input cube
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

% Interpolate illuminant and CMFs for the reference cube
illIP = interp1(ill(:,1), ill(:,2), bands_ref, 'spline');
CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands_ref, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,3), bands_ref, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,4), bands_ref, 'spline')];
sp_tristREF = CMFsIP .* illIP;

% Compute XYZ values for the reference cube
xyz_ref = (lincube_ref * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

%% Stratified Train-Test Split by Separating the Border Patches
% Determine the spatial indices of the patches in the color checker.
% We assume the color checker is arranged in an m x n grid.
[patch_rows, patch_cols] = ind2sub([m, n], (1:(m*n))');

% Define border indices as those in the first or last row or column
border_idx = find(patch_rows == 1 | patch_rows == m | patch_cols == 1 | patch_cols == n);
non_border_idx = setdiff((1:(m*n))', border_idx);

% Stratified sampling: 80% train, 20% test for border patches
perm_border = randperm(length(border_idx));
train_idx_border = border_idx(perm_border(1:round(0.8 * length(border_idx))));
test_idx_border  = border_idx(perm_border(round(0.8 * length(border_idx)) + 1:end));

% Stratified sampling: 80% train, 20% test for non-border patches
perm_non_border = randperm(length(non_border_idx));
train_idx_non_border = non_border_idx(perm_non_border(1:round(0.8 * length(non_border_idx))));
test_idx_non_border  = non_border_idx(perm_non_border(round(0.8 * length(non_border_idx)) + 1:end));

% Combine the indices from both groups
train_idx = sort([train_idx_border; train_idx_non_border]);
test_idx  = sort([test_idx_border; test_idx_non_border]);

% Define training data for XYZ space
xyz_input_train = xyz_input(train_idx, :);
xyz_ref_train   = xyz_ref(train_idx, :);


%% Polynomial Regression Model for XYZ (Second Order)
% Use a second-order expansion: [X, X.^2]
X_poly = [xyz_input_train, xyz_input_train.^2];
coeffs_xyz = (X_poly' * X_poly) \ (X_poly' * xyz_ref_train);

% Apply the model to the full dataset
X_poly_all = [xyz_input, xyz_input.^2];
corrected_xyz = X_poly_all * coeffs_xyz;

%% Convert to Lab and Apply Second Order Polynomial Correction
lab_input = xyz2lab(xyz_input);
lab_ref   = xyz2lab(xyz_ref);

lab_input_train = lab_input(train_idx, :);
lab_ref_train   = lab_ref(train_idx, :);

X_poly_lab = [lab_input_train, lab_input_train.^2];
coeffs_lab = (X_poly_lab' * X_poly_lab) \ (X_poly_lab' * lab_ref_train);

X_poly_lab_all = [lab_input, lab_input.^2];
corrected_lab = X_poly_lab_all * coeffs_lab;

%% Polynomial Correction in RGB Space
rgb_input = xyz2rgb(xyz_input);
rgb_ref   = xyz2rgb(xyz_ref);

rgb_input_train = rgb_input(train_idx, :);
rgb_ref_train   = rgb_ref(train_idx, :);

X_poly_rgb = [rgb_input_train, rgb_input_train.^2];
coeffs_rgb = (X_poly_rgb' * X_poly_rgb) \ (X_poly_rgb' * rgb_ref_train);

X_poly_rgb_all = [rgb_input, rgb_input.^2];
corrected_rgb = X_poly_rgb_all * coeffs_rgb;

%% Display Results for XYZ, Lab, and RGB
output_folder = '../results/error_maps';  % Output folder to save the error maps
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

lab_from_xyz_ref = xyz2lab(xyz_ref);
lab_from_xyz_cor = xyz2lab(corrected_xyz);

lab_from_rgb_ref = rgb2lab(rgb_ref);
lab_from_rgb_cor = rgb2lab(corrected_rgb);

evaluate_error(lab_from_xyz_ref, lab_from_xyz_cor, test_idx, m, n, 'XYZ', output_folder, '_error_poly2_grey.png');
evaluate_error(lab_ref, corrected_lab, test_idx, m, n, 'Lab', output_folder, '_error_poly2_grey.png');
evaluate_error(lab_from_rgb_ref, lab_from_rgb_cor, test_idx, m, n, 'RGB', output_folder, '_error_poly2_grey.png');




%% Plotting Train and Test patches
% rgb_input = xyz2rgb(xyz_input);
% border_rgb = rgb_input(train_idx, :);
% non_border_rgb = rgb_input(test_idx, :);
% 
% % Plot border patches
% figure;
% imshow(imresize(reshape(border_rgb, [], 1, 3), 50, 'nearest'));
% title('Trin Patches (RGB)');
% 
% % Plot non-border patches
% figure;
% imshow(imresize(reshape(non_border_rgb, [], 1, 3), 50, 'nearest'));
% title('Test Patches (RGB)');