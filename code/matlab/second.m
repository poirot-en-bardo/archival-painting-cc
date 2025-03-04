%% Define parameters
clc; clear; close all;
roof = double(intmax('uint16')); % Histogram normalization factor

% Load spectral data
ill = importdata('../../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../../data/CIE2degCMFs_1931.txt');
CMFs = CMFs_1931;
rng(10);

%% Select the reference and the cube to be corrected

cubeFile = "../../data/colorChecker_SG/cubes/cubeCC_120f-velvia-f8.hdr";
refFile = "../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

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

%% Polynomial Regression Model for XYZ (2nd order)
X_poly = [xyz_input_train, xyz_input_train.^2];
coeffs_xyz = (X_poly' * X_poly) \ (X_poly' * xyz_ref_train);

X_poly_all = [xyz_input, xyz_input.^2];
corrected_xyz = X_poly_all * coeffs_xyz;

%% Convert to Lab and apply 2nd order polynomial correction
lab_input = xyz2lab(xyz_input);
lab_ref = xyz2lab(xyz_ref);

lab_input_train = lab_input(train_idx, :);
lab_ref_train = lab_ref(train_idx, :);

X_poly_lab = [lab_input_train, lab_input_train.^2];
coeffs_lab = (X_poly_lab' * X_poly_lab) \ (X_poly_lab' * lab_ref_train);

X_poly_lab_all = [lab_input, lab_input.^2];
corrected_lab = X_poly_lab_all * coeffs_lab;

%% Polynomial correction in RGB space
rgb_input = xyz2rgb(xyz_input, 'ColorSpace','prophoto-rgb');
rgb_ref = xyz2rgb(xyz_ref, 'ColorSpace','prophoto-rgb');

rgb_input_train = rgb_input(train_idx, :);
rgb_ref_train = rgb_ref(train_idx, :);

X_poly_rgb = [rgb_input_train, rgb_input_train.^2];
coeffs_rgb = (X_poly_rgb' * X_poly_rgb) \ (X_poly_rgb' * rgb_ref_train);

X_poly_rgb_all = [rgb_input, rgb_input.^2];
corrected_rgb = X_poly_rgb_all * coeffs_rgb;

%% Display results for XYZ, Lab, and RGB
output_folder = '../results/error_maps';  % Output folder to save the error maps
mkdir(output_folder);  % Create the output folder if it doesn't exist

lab_from_xyz_ref = xyz2lab(xyz_ref);
lab_from_xyz_cor = xyz2lab(corrected_xyz);

lab_from_rgb_ref = rgb2lab(rgb_ref, 'ColorSpace','prophoto-rgb');
lab_from_rgb_cor = rgb2lab(corrected_rgb, 'ColorSpace','prophoto-rgb');

evaluate_error(lab_from_xyz_ref, lab_from_xyz_cor, test_idx, m, n, 'XYZ', output_folder, '_error_poly2.png');
evaluate_error(lab_ref, corrected_lab, test_idx, m, n, 'Lab', output_folder, '_error_poly2.png');
evaluate_error(lab_from_rgb_ref, lab_from_rgb_cor, test_idx, m, n, 'RGB', output_folder, '_error_poly2.png');

