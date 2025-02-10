%% Define parameters
clc; clear; close all;

% Load spectral data
ill = importdata('../data/CIE_D65.txt'); % Illuminant
CMFs = importdata('../data/CIE2degCMFs_1931.txt'); % Color matching functions
rng(10);

%% Load the reference and test cubes
cubeFile = "/Volumes/School/Thesis/data/colorChecker_SG/cubes/cubeCC_fuji-frame4.hdr";
refFile = "/Volumes/School/Thesis/data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd_in] = size(inCUBE);
lincube = reshape(inCUBE,[],bd_in);

hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
lincube_ref = reshape(refCUBE,[],size(refCUBE,3));

% Process XYZ values
xyz_input = process_xyz(lincube, bands, ill, CMFs);
xyz_ref = process_xyz(lincube_ref, bands_ref, ill, CMFs);

%% Neural Network Model for Color Correction
hiddenLayerSize = 5; % Simple single hidden layer with 5 neurons
net = feedforwardnet(hiddenLayerSize, 'trainlm'); % Fast training method

% Set training parameters
net.trainParam.epochs = 300;  % Train for fewer epochs
net.trainParam.goal = 1e-4;   % Stop early if MSE is small

% Train-test split
num_samples = size(xyz_input, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx = perm(round(0.8*num_samples)+1:end);

% Normalize data (important for NN)
X_min = min(xyz_input, [], 1);
X_max = max(xyz_input, [], 1);
Y_min = min(xyz_ref, [], 1);
Y_max = max(xyz_ref, [], 1);

X_train = (xyz_input(train_idx, :) - X_min) ./ (X_max - X_min);
Y_train = (xyz_ref(train_idx, :) - Y_min) ./ (Y_max - Y_min);
X_test = (xyz_input(test_idx, :) - X_min) ./ (X_max - X_min);
Y_test = (xyz_ref(test_idx, :) - Y_min) ./ (Y_max - Y_min);

% Train the neural network
net = train(net, X_train', Y_train');

% Apply correction to the entire input cube
corrected_xyz = net(((xyz_input - X_min) ./ (X_max - X_min))')';

% Rescale output back to original range
corrected_xyz = corrected_xyz .* (Y_max - Y_min) + Y_min;

% Evaluate performance
evaluate_error(X_test, Y_test, net, xyz_ref, corrected_xyz, m, n);


%% Function to Process XYZ Data
function xyz = process_xyz(lincube, bands, ill, CMFs)
    illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
    CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
    sp_tristREF = CMFsIP .* illIP;
    xyz = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);
end

%% Function to Evaluate ΔE2000 Error
function evaluate_error(X_test, Y_test, net, Y, corrected_xyz, m, n)
    % Ensure correct Nx3 format for xyz2lab conversion
    lab_ref = xyz2lab(reshape(Y_test', [], 3));
    lab_corrected = xyz2lab(reshape(net(X_test')', [], 3));
    
    deltaE2000_errors = deltaE2000(lab_corrected, lab_ref);
    
    lab_full_ref = xyz2lab(reshape(Y, [], 3));
    lab_full_corrected = xyz2lab(reshape(corrected_xyz, [], 3));
    error_map = reshape(deltaE2000(lab_full_corrected, lab_full_ref), m, n);
    
    figure;
    imagesc(error_map);
    colormap(jet);
    colorbar;
    clim([0 10]);
    title('\DeltaE2000 Error Map');
    grid off;
    
    disp(['Mean \DeltaE2000 Error: ', num2str(mean(deltaE2000_errors))]);
    disp(['Max \DeltaE2000 Error: ', num2str(max(deltaE2000_errors))]);
end

