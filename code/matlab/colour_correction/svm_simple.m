%% Define parameters
clc; clear; close all;
roof = double(intmax('uint16')); % Histogram normalization factor

% Load spectral data
ill = importdata('../../../data/CIE_D65.txt'); % Illuminant
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

%% Train-test split (80% training, 20% testing)
num_samples = size(xyz_input, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx = perm(round(0.8*num_samples)+1:end);

xyz_input_train = xyz_input(train_idx, :);
xyz_ref_train = xyz_ref(train_idx, :);
lab_input_train = xyz2lab(xyz_input_train);
lab_ref_train = xyz2lab(xyz_ref_train);
lab_input = xyz2lab(xyz_input);

%% Train an SVM for XYZ-to-XYZ mapping (for each channel)
% X channel (X -> X)
svm_model_X = fitrsvm(lab_input_train, lab_ref_train(:, 1), 'KernelFunction', 'gaussian', 'Standardize', true);

% Y channel (Y -> Y)
svm_model_Y = fitrsvm(lab_input_train, lab_ref_train(:, 2), 'KernelFunction', 'gaussian', 'Standardize', true);

% Z channel (Z -> Z)
svm_model_Z = fitrsvm(lab_input_train, lab_ref_train(:, 3), 'KernelFunction', 'gaussian', 'Standardize', true);

%% Predict XYZ values for the test data using the trained SVMs
xyz_pred_X = predict(svm_model_X, lab_input);
xyz_pred_Y = predict(svm_model_Y, lab_input);
xyz_pred_Z = predict(svm_model_Z, lab_input);

% Combine the predicted channels into one XYZ matrix
lab_pred = [xyz_pred_X, xyz_pred_Y, xyz_pred_Z];


%% Display results for XYZ or Lab (X, Y, Z or L, a, b)
output_folder = '../../../results/error_maps';  % Output folder to save the error maps
mkdir(output_folder);  % Create the output folder if it doesn't exist

% If you want Lab comparison, you can convert:
lab_from_ref = xyz2lab(xyz_ref);
% lab_pred = xyz2lab(lab_pred);

[~, img_name, ~] = fileparts(cubeFile);

% Evaluate and display error maps
evaluate_error(lab_from_ref, lab_pred, test_idx, m, n, 'Lab_SVM_', output_folder, img_name + "_error_svm.png");
