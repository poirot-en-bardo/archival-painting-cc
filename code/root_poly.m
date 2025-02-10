%% Define parameters
clc; clear; close all;
roof = double(intmax('uint16'));

histoFACT = 200; % Histogram normalization factor

% Load spectral data
ill = importdata('../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../data/CIE2degCMFs_1931.txt');
CMFs_2006 = importdata('../data/CIE2degCMFs_2006.txt');
CMFs = CMFs_1931;
rng(10);

%% Select the reference and the cube to be corrected
cube_folder = '../data/colorChecker_SG/colorChecker_SG_Elias';

% [cubeFileName, cubePath] = uigetfile(fullfile(cube_folder, '*.hdr'), 'Select Cube to be Corrected');
% cubeFile = fullfile(cubePath, cubeFileName);
cubeFile = "/Volumes/School/Thesis/data/colorChecker_SG/cubes/cubeCC_fuji-frame4.hdr";

% [refFileName, refPath] = uigetfile(fullfile(cube_folder, '*.hdr'), 'Select Reference Cube');
% refFile = fullfile(refPath, refFileName);
refFile = "/Volumes/School/Thesis/data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

%% Load and process the cube to be corrected
hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd] = size(inCUBE);
lincube = reshape(inCUBE,[],bd);

% Interpolate illuminant and CMF to captured wavelengths
illIP = interp1(ill(:,1),ill(:,2),bands,'spline');
CMFsIP = [interp1(CMFs(:,1),CMFs(:,2),bands,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,3),bands,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,4),bands,'spline')];
sp_tristREF = CMFsIP .* illIP;

% Compute XYZ values of the cube to be corrected
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2),1);

%%% Load and process the reference cube
hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
[m, n, bd] = size(refCUBE);
lincube_ref = reshape(refCUBE,[],bd);

illIP = interp1(ill(:,1),ill(:,2),bands_ref,'spline');
CMFsIP = [interp1(CMFs(:,1),CMFs(:,2),bands_ref,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,3),bands_ref,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,4),bands_ref,'spline')];
sp_tristREF = CMFsIP .* illIP;

% Compute XYZ values for the reference cube
xyz_ref = (lincube_ref * sp_tristREF) ./ sum(sp_tristREF(:,2),1);


%% Root-Polynomial Regression Model for Color Correction
degree = 2;

% Define root-polynomial features
X = [xyz_input, sqrt(xyz_input(:,1).*xyz_input(:,2)), ...
     sqrt(xyz_input(:,2).*xyz_input(:,3)), sqrt(xyz_input(:,1).*xyz_input(:,3))]; % Root terms

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

% Fit root-polynomial regression model (using least squares)
% coeffs_root = (X_train' * X_train) \ (X_train' * Y_train);
coeffs_root = (Y_train' * X_train * pinv(X_train' * X_train))'; % Moore-Penrose Inverse

% Apply correction to the entire input cube
corrected_xyz = X * coeffs_root;
evaluate_error(X_test, Y_test, coeffs_root, Y, corrected_xyz, m, n);




%%
function evaluate_error(X_test, Y_test, coeffs, Y, corrected_xyz, m, n)
    % Compute Delta E2000 Error
    lab_ref = xyz2lab(Y_test);
    lab_corrected = xyz2lab((X_test * coeffs));
    
    deltaE2000_errors = deltaE2000(lab_corrected, lab_ref);
    
    % Error map visualization
    lab_full_ref = xyz2lab(Y);
    lab_full_corrected = xyz2lab(corrected_xyz);
    error_map = reshape(deltaE2000(lab_full_corrected, lab_full_ref), m, n);
    
    figure;
    
     % Display the error map with a colorbar
    imagesc(error_map);
    colormap(jet);
    colorbar;
    clim([0 10]); % Adjust the color axis limits as needed
    title('ΔE2000 Error Map');
    grid off;
    
    % Display summary statistics
    disp(['Mean ΔE2000 Error: ', num2str(mean(deltaE2000_errors))]);
    disp(['Max ΔE2000 Error: ', num2str(max(deltaE2000_errors))]);
end
