%% Define some parameters
clc; clear; close all;
roof = double(intmax('uint16'));

histoFACT = 200; % Factor used to find the edges of the histograms
rng(10);

ill = importdata('../data/CIE_D65.txt');
CMFs_1931 = importdata('../data/CIE2degCMFs_1931.txt');
CMFs_2006 = importdata('../data/CIE2degCMFs_2006.txt');
CMFs = CMFs_1931;

%% Select the reference and the cube to be corrected
cube_folder = '../data/colorChecker_SG/colorChecker_SG_Elias';

% cube to be corrected
f = msgbox('Select the cube to be corrected');
movegui(f, 'north');
[cubeFileName, cubePath] = uigetfile(fullfile(cube_folder, '*.hdr'));
close(f);


cubeFile = fullfile(cubePath, cubeFileName);

% reference cube
f = msgbox('Select the reference cube');
movegui(f, 'north');
[refFileName, refPath] = uigetfile(fullfile(cube_folder, '*.hdr'), 'Select Reference Cube');
close(f);

refFile = fullfile(refPath, refFileName);

%% Load and process the cube to be corrected
hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd] = size(inCUBE);
lincube = reshape(inCUBE,[],bd);

% Interpolate illuminant and CMF to captured wavelengths
illIP = interp1(ill(:,1),ill(:,2),bands,'spline');
CMFsIP = [interp1(CMFs_1931(:,1),CMFs_1931(:,2),bands,'spline'), ...
          interp1(CMFs_1931(:,1),CMFs_1931(:,3),bands,'spline'), ...
          interp1(CMFs_1931(:,1),CMFs_1931(:,4),bands,'spline')];
sp_tristREF = CMFsIP .* illIP;

% Compute XYZ values of the cube to be corrected
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2),1);

%% Load and process the reference cube
hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
[m, n, bd] = size(refCUBE);
lincube_ref = reshape(refCUBE,[],bd);

illIP = interp1(ill(:,1),ill(:,2),bands_ref,'spline');
CMFsIP = [interp1(CMFs_1931(:,1),CMFs_1931(:,2),bands_ref,'spline'), ...
          interp1(CMFs_1931(:,1),CMFs_1931(:,3),bands_ref,'spline'), ...
          interp1(CMFs_1931(:,1),CMFs_1931(:,4),bands_ref,'spline')];
sp_tristREF = CMFsIP .* illIP;

% Compute XYZ values for the reference cube
xyz_ref = (lincube_ref * sp_tristREF) ./ sum(sp_tristREF(:,2),1);

%% Polynomial Regression Model for Color Correction
% 2nd degree polynomial model 
degree = 2;
X = [xyz_input, xyz_input.^2];  % Include quadratic terms
Y = xyz_ref;  % target reference values

% Train-test split (80% training, 20% testing)
num_samples = size(X, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx = perm(round(0.8*num_samples)+1:end);

X_train = X(train_idx, :);
Y_train = Y(train_idx, :);
X_test = X(test_idx, :);
Y_test = Y(test_idx, :);

% fit polynomial model (using least squares)
coeffs = (X_train' * X_train) \ (X_train' * Y_train);

% apply correction to the entire input cube
corrected_xyz = X * coeffs;

%%
%% 3rd degree 

% Define 3rd-degree polynomial features
X = [xyz_input, xyz_input.^2, xyz_input.^3]; % Include quadratic and cubic terms

Y = xyz_ref; % Target reference values

% Train-test split (80% training, 20% testing)
num_samples = size(X, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx = perm(round(0.8*num_samples)+1:end);

X_train = X(train_idx, :);
Y_train = Y(train_idx, :);
X_test = X(test_idx, :);
Y_test = Y(test_idx, :);

% Fit the polynomial regression model (using least squares)
coeffs_poly = (X_train' * X_train) \ (X_train' * Y_train);

% Apply correction to the entire input cube
corrected_xyz = [xyz_input, xyz_input.^2, xyz_input.^3] * coeffs_poly;

% Evaluate error
evaluate_error(X_test, Y_test, coeffs_poly, Y, corrected_xyz, m, n);

%% 2nd with interaction
% Define 2nd-degree polynomial features with interaction terms
X = [xyz_input, ... 
     xyz_input.^2, ... % Quadratic terms
     xyz_input(:,1).*xyz_input(:,2), ... % Interaction between feature 1 and 2
     xyz_input(:,1).*xyz_input(:,3), ... % Interaction between feature 1 and 3
     xyz_input(:,2).*xyz_input(:,3)]; % Interaction between feature 2 and 3

Y = xyz_ref; % Target reference values

% Train-test split (80% training, 20% testing)
num_samples = size(X, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx = perm(round(0.8*num_samples)+1:end);

X_train = X(train_idx, :);
Y_train = Y(train_idx, :);
X_test = X(test_idx, :);
Y_test = Y(test_idx, :);

% Fit the polynomial regression model (using least squares)
coeffs_poly_interaction = (X_train' * X_train) \ (X_train' * Y_train);

% Apply correction to the entire input cube using the same expanded X matrix
corrected_xyz = X * coeffs_poly_interaction;

% Evaluate error
evaluate_error(X_test, Y_test, coeffs_poly_interaction, Y, corrected_xyz, m, n);

%% 3rd with interaction
% Define 3rd-degree polynomial features with interaction terms
X = [xyz_input, ...
     xyz_input.^2, ... % Quadratic terms
     xyz_input.^3, ... % Cubic terms
     xyz_input(:,1).*xyz_input(:,2), ... % Interaction between feature 1 and 2
     xyz_input(:,1).*xyz_input(:,3), ... % Interaction between feature 1 and 3
     xyz_input(:,2).*xyz_input(:,3), ... % Interaction between feature 2 and 3
     xyz_input(:,1).^2.*xyz_input(:,2), ... % Interaction between square of feature 1 and feature 2
     xyz_input(:,1).*xyz_input(:,2).^2, ... % Interaction between feature 1 and square of feature 2
     xyz_input(:,2).^2.*xyz_input(:,3), ... % Interaction between square of feature 2 and feature 3
     xyz_input(:,1).*xyz_input(:,2).*xyz_input(:,3)]; % Interaction between all features

Y = xyz_ref; % Target reference values

% Train-test split (80% training, 20% testing)
num_samples = size(X, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx = perm(round(0.8*num_samples)+1:end);

X_train = X(train_idx, :);
Y_train = Y(train_idx, :);
X_test = X(test_idx, :);
Y_test = Y(test_idx, :);

% Fit the polynomial regression model (using least squares)
coeffs_poly_interaction = (X_train' * X_train) \ (X_train' * Y_train);

% Apply correction to the entire input cube using the same expanded X matrix
corrected_xyz = X * coeffs_poly_interaction;

% Evaluate error
evaluate_error(X_test, Y_test, coeffs_poly_interaction, Y, corrected_xyz, m, n);
%% Compute Delta E2000 Error
lab_ref = xyz2lab(Y_test);
lab_corrected = xyz2lab((X_test) * coeffs);

deltaE2000_errors = deltaE2000(lab_corrected, lab_ref);

% error map for visualisation
lab_full_ref = xyz2lab(Y);
lab_full_corrected = xyz2lab(corrected_xyz);
error_map = reshape(deltaE2000(lab_full_corrected, lab_full_ref), m, n);

%% Visualize the Error Map
figure;
imshow(imresize(error_map, 100, 'nearest'), ...
    'DisplayRange', [0 10], 'Colormap', jet(255));
colorbar;
title('ΔE2000 Error Map');

disp(['Mean ΔE2000 Error: ', num2str(mean(deltaE2000_errors))]);
disp(['Max ΔE2000 Error: ', num2str(max(deltaE2000_errors))]);


%% Convert XYZ to RGB
rgb_image = xyz2rgb(xyz_input, 'ColorSpace', 'adobe-rgb-1998');
rgb_image = uint16(reshape(rgb_image, m, n, 3) .* 65535); % Scale to 16-bit

% Create RGB copies for training and test visualization
train_rgb = rgb_image;
test_rgb = rgb_image;

% Create binary masks for training and testing patches
train_mask = zeros(m, n);
test_mask = zeros(m, n);

% Mark training and testing positions in the masks
train_mask(train_idx) = 1; % Mark training samples
test_mask(test_idx) = 1;   % Mark testing samples

% Apply the masks by setting non-selected regions to white (or dim them)
for c = 1:3  % Loop over RGB channels
    train_rgb(:,:,c) = train_rgb(:,:,c) .* uint16(train_mask) + uint16(~train_mask) * 65535;
    test_rgb(:,:,c) = test_rgb(:,:,c) .* uint16(test_mask) + uint16(~test_mask) * 65535;
end

% Resize for better visibility
scale_factor = 100; % Adjust the factor as needed
train_rgb_resized = imresize(train_rgb, scale_factor, 'nearest');
test_rgb_resized = imresize(test_rgb, scale_factor, 'nearest');

% Display training patches
figure, imshow(train_rgb_resized);
title('Training Patches Highlighted');

% Display testing patches
figure, imshow(test_rgb_resized);
title('Testing Patches Highlighted');