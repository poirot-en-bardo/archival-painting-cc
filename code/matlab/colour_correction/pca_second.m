%% Define parameters
clc; clear; close all;
roof = double(intmax('uint16')); % Histogram normalization factor

% Load spectral data
ill = importdata('../../../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../../../data/CIE2degCMFs_1931.txt');
CMFs = CMFs_1931;
rng(10);

%% Select the reference and the cube to be corrected

cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_120f-velvia-f8.hdr";
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
[m, n, bd] = size(refCUBE);
lincube_ref = reshape(refCUBE, [], bd);

% Interpolate illuminant and CMFs
illIP = interp1(ill(:,1),ill(:,2),bands,'spline');
CMFsIP = [interp1(CMFs(:,1),CMFs(:,2),bands,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,3),bands,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,4),bands,'spline')];
sp_tristREF = CMFsIP .* illIP;

% Compute XYZ values for the input and reference cubes
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2),1);

illIP = interp1(ill(:,1),ill(:,2),bands_ref,'spline');
CMFsIP = [interp1(CMFs(:,1),CMFs(:,2),bands_ref,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,3),bands_ref,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,4),bands_ref,'spline')];
sp_tristREF = CMFsIP .* illIP;

xyz_ref = (lincube_ref * sp_tristREF) ./ sum(sp_tristREF(:,2),1);
lab_input = xyz2lab(xyz_input);
lab_ref = xyz2lab(xyz_ref);
rgb_input = xyz2rgb(xyz_input, 'ColorSpace','prophoto-rgb');
rgb_ref = xyz2rgb(xyz_ref, 'ColorSpace','prophoto-rgb');

%% Train-test split (80% training, 20% testing)
num_samples = size(xyz_input, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx = perm(round(0.8*num_samples)+1:end);

xyz_input_train = xyz_input(train_idx, :);
xyz_ref_train = xyz_ref(train_idx, :);

lab_input_train = lab_input(train_idx, :);
lab_ref_train = lab_ref(train_idx, :);

rgb_input_train = rgb_input(train_idx, :);
rgb_ref_train = rgb_ref(train_idx, :);


%% Apply PCA to reduce dimensions
reshapedData = reshape(inCUBE, [], size(inCUBE, 3));
[coeff, score, ~, ~, explained] = pca(reshapedData);

% Select only the first 10 principal components
numPCs = 9;
reducedData = score(:, 1:numPCs);
reducedData_train = reducedData(train_idx, :);

%% Polynomial Regression using PCA components
X_poly_pc = [reducedData_train, reducedData_train.^2];
coeffs_pc = (X_poly_pc' * X_poly_pc) \ (X_poly_pc' * rgb_ref_train);

% Predictions
X_poly_pc_all = [reducedData, reducedData.^2];
predictions_pc = X_poly_pc_all * coeffs_pc;

% Transform back to original space (optional)
% originalSpacePredictions = predictions_pc * coeff(:, 1:numPCs)';

%% Display results
output_folder = '../../../results/error_maps';
mkdir(output_folder);

lab_from_xyz_ref = xyz2lab(xyz_ref);
lab_from_xyz_cor = rgb2lab(predictions_pc, 'ColorSpace','prophoto-rgb');

[~, img_name, ~] = fileparts(cubeFile);

evaluate_error(lab_from_xyz_ref, lab_from_xyz_cor, test_idx, m, n, 'PCA_', output_folder, img_name + "_pca_spectra_rgb_poly2.png");



