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

%% Process XYZ values
disp('Processing XYZ values...');
xyz_input = process_xyz(lincube, bands, ill, CMFs);
xyz_ref = process_xyz(lincube_ref, bands_ref, ill, CMFs);

%% Normalize Data
disp('Normalizing data...');
X_min = min(xyz_input, [], 1);
X_max = max(xyz_input, [], 1);
Y_min = min(xyz_ref, [], 1);
Y_max = max(xyz_ref, [], 1);

xyz_input_norm = (xyz_input - X_min) ./ (X_max - X_min);
xyz_ref_norm = (xyz_ref - Y_min) ./ (Y_max - Y_min);

%% Train-Test Split
num_samples = size(xyz_input_norm, 1);
train_ratio = 0.8;
train_idx = randperm(num_samples, round(train_ratio * num_samples));
test_idx = setdiff(1:num_samples, train_idx);

X_train = xyz_input_norm(train_idx, :);
Y_train = xyz_ref_norm(train_idx, :);
X_test = xyz_input_norm(test_idx, :);
Y_test = xyz_ref_norm(test_idx, :);

%% Define and Train Neural Network
disp('Training neural network...');
hiddenLayerSize = 5;
net = feedforwardnet(hiddenLayerSize, 'trainlm');
net.trainParam.epochs = 300;
net.trainParam.lr = 0.01;
net.trainParam.showWindow = false; % Disable training window for speed

% Train network
net = train(net, X_train', Y_train');

%% Evaluate Model on Test Set
disp('Evaluating the model...');
Y_pred_test = net(X_test')';
mse_loss = mean((Y_test - Y_pred_test).^2, 'all');
fprintf('Test MSE Loss: %.4f\n', mse_loss);

%% Apply Correction to Full Image
disp('Applying correction...');
corrected_xyz_norm = net(xyz_input_norm')';
corrected_xyz = corrected_xyz_norm .* (Y_max - Y_min) + Y_min;

%% Compute ΔE2000 Error Map
disp('Computing ΔE2000 error map...');
lab_ref = xyz2lab(xyz_ref);
lab_corr = xyz2lab(corrected_xyz);

error_map = deltaE2000(lab_ref, lab_corr);
mean_error = mean(error_map(:));
max_error = max(error_map(:));
fprintf('Mean ΔE2000 Error: %.4f\n', mean_error);
fprintf('Max ΔE2000 Error: %.4f\n', max_error);


%% Display Error Map
figure;
imagesc(reshape(error_map, m, n));
colormap jet;
colorbar;
clim([0 10]); % Set color limits for better visualization
title('\DeltaE2000 Error Map');
xlabel('Width');
ylabel('Height');

%% Function to Process XYZ Data
function xyz = process_xyz(lincube, bands, ill, CMFs)
    illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
    CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
    sp_tristREF = CMFsIP .* illIP;
    xyz = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);
end

