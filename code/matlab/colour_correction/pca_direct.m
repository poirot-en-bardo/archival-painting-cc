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

%% Train-test split
num_samples = size(input_pca, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx = perm(round(0.8*num_samples)+1:end);

input_pca_train = input_pca(train_idx, :);
ref_pca_train = ref_pca(train_idx, :);

%% Polynomial Regression between PCA components
X_poly_pca_train = poly3_features(input_pca_train);  
coeffs_pca = pinv(X_poly_pca_train) * ref_pca_train;  

% Apply the regression to the full dataset
X_poly_pca = poly3_features(input_pca);  
pca_corrected = X_poly_pca * coeffs_pca;

%% Reconstruct full spectral data from corrected PCA
reconstructed_spectrum = pca_corrected * coeff_ref(:, 1:num_components_ref)';

%% Convert to XYZ
illIP = interp1(ill(:,1),ill(:,2),bands_ref,'spline');
CMFsIP = [interp1(CMFs(:,1),CMFs(:,2),bands_ref,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,3),bands_ref,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,4),bands_ref,'spline')];
sp_tristREF = CMFsIP .* illIP;

xyz_corrected = (reconstructed_spectrum * sp_tristREF) ./ sum(sp_tristREF(:,2),1);
xyz_ref = (lincube_ref * sp_tristREF) ./ sum(sp_tristREF(:,2),1);

%% Convert to RGB
rgb_corrected = xyz2rgb(xyz_corrected, 'ColorSpace', 'prophoto-rgb');
rgb_ref = xyz2rgb(xyz_ref, 'ColorSpace', 'prophoto-rgb');

%% Convert to Lab
lab_corrected = xyz2lab(xyz_corrected);
lab_ref = xyz2lab(xyz_ref);

%% Evaluate Error
output_folder = '../../../results/error_maps';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

[~, img_name, ~] = fileparts(cubeFile);

evaluate_error(lab_ref, lab_corrected, test_idx, m, n, 'PCA_', output_folder, img_name + "_pca_poly3.png");

%% Polynomial Feature Function
function X_poly = poly3_features(input_data)
    num_cols = size(input_data, 2);
    X_poly = []; 
    for i = 1:num_cols
        X_poly = [X_poly, input_data(:, i)];
    end
    for i = 1:num_cols
        X_poly = [X_poly, input_data(:, i).^2];
        for j = i+1:num_cols
            X_poly = [X_poly, input_data(:, i).*input_data(:, j)];
        end
    end
    for i = 1:num_cols
        X_poly = [X_poly, input_data(:, i).^3];
        for j = i+1:num_cols
            X_poly = [X_poly, input_data(:, i).^2.*input_data(:, j)];
            X_poly = [X_poly, input_data(:, i).*input_data(:, j).^2];
            for k = j+1:num_cols
                X_poly = [X_poly, input_data(:, i).*input_data(:, j).*input_data(:, k)];
            end
        end
    end
end