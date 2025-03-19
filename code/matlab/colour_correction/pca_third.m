clc; clear; close all;
roof = double(intmax('uint16')); % Normalization factor

% Load spectral data
ill = importdata('../../../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../../../data/CIE2degCMFs_1931.txt');
CMFs = CMFs_1931;
rng(10);

%% Load and process the hyperspectral cube
cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_120f-velvia-f8.hdr";
refFile = "../../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd] = size(inCUBE);
lincube = reshape(inCUBE, [], bd); % Flatten spectral data


hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
[m, n, bd_ref] = size(refCUBE);
lincube_ref = reshape(refCUBE, [], bd_ref);

%% PCA on spectral data 
[coeff_input, score_input, ~, ~, explained_input] = pca(lincube);
num_components = find(cumsum(explained_input) >= 99, 1); % 99% variance

[coeff_ref, score_ref, ~, ~, explained_ref] = pca(lincube_ref);
num_components_ref = find(cumsum(explained_ref) >= 99, 1);

% Reduce the spectral data using PCA
input_pca = score_input(:, 1:num_components);
ref_pca = score_ref(:, 1:num_components_ref);

%% Reconstruct the Spectral Data from PCA
reconstructed_input_spectrum = input_pca * coeff_input(:, 1:num_components)';
reconstructed_ref_spectrum = ref_pca * coeff_ref(:, 1:num_components_ref)';

%% Convert Reconstructed Spectral Data to XYZ
% Interpolate illuminant and CMFs
illIP = interp1(ill(:,1),ill(:,2),bands,'spline');
CMFsIP = [interp1(CMFs(:,1),CMFs(:,2),bands,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,3),bands,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,4),bands,'spline')];
sp_tristREF = CMFsIP .* illIP;

% Compute XYZ values for the input and reference cubes
xyz_input = (reconstructed_input_spectrum * sp_tristREF) ./ sum(sp_tristREF(:,2),1);

illIP = interp1(ill(:,1),ill(:,2),bands_ref,'spline');
CMFsIP = [interp1(CMFs(:,1),CMFs(:,2),bands_ref,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,3),bands_ref,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,4),bands_ref,'spline')];
sp_tristREF = CMFsIP .* illIP;

xyz_ref = (reconstructed_ref_spectrum * sp_tristREF) ./ sum(sp_tristREF(:,2),1);


%% Train-Test Split for Regression
num_samples = size(xyz_input, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx = perm(round(0.8*num_samples)+1:end);

xyz_input_train = xyz_input(train_idx, :);
xyz_ref_train = xyz_ref(train_idx, :);


%% Regression XYZ
X_poly_xyz_train = poly3_features(xyz_input_train);  
coeffs_xyz = pinv(X_poly_xyz_train) * xyz_ref_train;  

% Apply the regression to the full dataset
X_poly_xyz = poly3_features(xyz_input);  
corrected_xyz = X_poly_xyz * coeffs_xyz;


%% Convert to Lab and apply 2nd order polynomial correction
lab_input = xyz2lab(xyz_input);
lab_ref = xyz2lab(xyz_ref);

lab_input_train = lab_input(train_idx, :);
lab_ref_train = lab_ref(train_idx, :);

X_poly_lab_train =  poly3_features(lab_input_train);
coeffs_lab = pinv(X_poly_lab_train) * lab_ref_train;

X_poly_lab = poly3_features(lab_input);
corrected_lab = X_poly_lab * coeffs_lab;

%% ----------------- RGB-Based Regression  -----------------
% Convert input and reference XYZ to RGB (using prophoto-rgb)
rgb_input = xyz2rgb(xyz_input, 'ColorSpace', 'prophoto-rgb');
rgb_ref   = xyz2rgb(xyz_ref,   'ColorSpace', 'prophoto-rgb');

rgb_input_train = rgb_input(train_idx, :);
rgb_ref_train   = rgb_ref(train_idx, :);

% Compute 3rd-degree polynomial features for RGB training data
X_poly_rgb_train = poly3_features(rgb_input_train);  

% Compute coefficients using pseudoinverse
coeffs_rgb = pinv(X_poly_rgb_train) * rgb_ref_train;  


% Apply the regression to the full dataset
X_poly_rgb = poly3_features(rgb_input);  
corrected_rgb = X_poly_rgb * coeffs_rgb;



% Convert corrected outputs to Lab for error evaluation
lab_from_xyz_ref = xyz2lab(xyz_ref);
lab_from_xyz_cor = xyz2lab(corrected_xyz);

lab_from_rgb_ref = rgb2lab(rgb_ref, 'ColorSpace', 'prophoto-rgb');
lab_from_rgb_cor = rgb2lab(corrected_rgb, 'ColorSpace', 'prophoto-rgb');

output_folder = '../../../results/error_maps';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

[~, img_name, ~] = fileparts(cubeFile);

evaluate_error(lab_from_xyz_ref, lab_from_xyz_cor, test_idx, m, n, 'XYZ_', output_folder, img_name + "_pca_poly3.png");
evaluate_error(lab_ref, corrected_lab, test_idx, m, n, 'Lab_', output_folder, img_name + "_pca_poly3.png");
evaluate_error(lab_from_rgb_ref, lab_from_rgb_cor, test_idx, m, n, 'RGB_', output_folder, img_name + "_pca_poly3.png");
%% Polynomial Feature Function
function X_poly = poly3_features(input_data)
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);
    
    X_poly = [a, b, c, a.^2, b.^2, c.^2, a.*b, a.*c, b.*c, ...
              a.^3, b.^3, c.^3, a.^2.*b, a.^2.*c, b.^2.*a, b.^2.*c, c.^2.*a, c.^2.*b, a.*b.*c];
end
