%% Define parameters
clc; clear; close all;

% Load spectral data
ill = importdata('../data/CIE_D65.txt'); % Illuminant
CMFs = importdata('../data/CIE2degCMFs_1931.txt'); % Color matching functions
rng(10);

%% Load the reference and test cubes
cubeFile = "../data/colorChecker_SG/cubes/cubeCC_fuji-frame4.hdr";
refFile = "../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

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
hiddenLayerSize = 5; % Single hidden layer with 5 neurons
net = feedforwardnet(hiddenLayerSize, 'trainlm'); % Fast training method


% Training parameters
num_epochs = 300;
learning_rate = 0.01;  % Small learning rate for stability

% Train-test split
num_samples = size(xyz_input, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx = perm(round(0.8*num_samples)+1:end);

% Normalize data
X_min = min(xyz_input, [], 1);
X_max = max(xyz_input, [], 1);
Y_min = min(xyz_ref, [], 1);
Y_max = max(xyz_ref, [], 1);

X_train = (xyz_input(train_idx, :) - X_min) ./ (X_max - X_min);
Y_train = (xyz_ref(train_idx, :) - Y_min) ./ (Y_max - Y_min);
X_test = (xyz_input(test_idx, :) - X_min) ./ (X_max - X_min);
Y_test = (xyz_ref(test_idx, :) - Y_min) ./ (Y_max - Y_min);

% Transpose X_train to match expected input size
X_train = X_train'; % Size should be (num_features, num_samples)
Y_train = Y_train'; % Size should be (num_targets, num_samples)
% Set the network input size explicitly
net.inputs{1}.size = size(X_train, 1);
% Confirm the size of the training data
disp(['Input size (X_train): ', num2str(size(X_train))]);  % Expected: (num_features, num_samples)
disp(['Output size (Y_train): ', num2str(size(Y_train))]);  % Expected: (num_targets, num_samples)
disp(['Network input size: ', num2str(net.inputs{1}.size)]);

%% Custom Training Loop Using ΔE2000 Loss
for epoch = 1:num_epochs
    % Forward Pass: Get predicted XYZ from the NN
    Y_pred = net(X_train);

    % Convert XYZ to Lab
    lab_ref = xyz2lab(Y_train');
    lab_pred = xyz2lab(Y_pred);

    % Compute ΔE2000 Loss
    loss = mean(deltaE2000(lab_pred, lab_ref));

    % Compute Gradients Numerically
    gradients = computeGradient(net, X_train, lab_ref);

    % Update Weights using Gradient Descent
    net.IW{1} = net.IW{1} - learning_rate * gradients;

    % Display loss
    disp(['Epoch: ', num2str(epoch), '  ΔE2000 Loss: ', num2str(loss)]);
end

%% Apply Correction and Evaluate Performance
corrected_xyz = net(((xyz_input - X_min) ./ (X_max - X_min))')';

% Rescale output back to original range
corrected_xyz = corrected_xyz .* (Y_max - Y_min) + Y_min;

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

%% Function to Compute Numerical Gradients for ΔE2000
function gradients = computeGradient(net, X_train, lab_ref)
    epsilon = 1e-5;
    gradients = zeros(size(net.IW{1}));

    for i = 1:numel(net.IW{1})
        net_temp = net;
        
        % Compute f(x+h)
        net_temp.IW{1}(i) = net_temp.IW{1}(i) + epsilon;
        lab_pred_plus = xyz2lab(net_temp(X_train')');
        loss_plus = mean(deltaE2000(lab_pred_plus, lab_ref));

        % Compute f(x-h)
        net_temp.IW{1}(i) = net_temp.IW{1}(i) - 2 * epsilon;
        lab_pred_minus = xyz2lab(net_temp(X_train')');
        loss_minus = mean(deltaE2000(lab_pred_minus, lab_ref));

        % Compute numerical gradient
        gradients(i) = (loss_plus - loss_minus) / (2 * epsilon);
    end
end

%% Function to Evaluate ΔE2000 Error
function evaluate_error(X_test, Y_test, net, Y, corrected_xyz, m, n)
    % Convert to Lab for ΔE2000
    lab_ref = xyz2lab(reshape(Y_test', [], 3));
    lab_corrected = xyz2lab(reshape(net(X_test')', [], 3));
    
    % Compute ΔE2000 Errors
    deltaE2000_errors = deltaE2000(lab_corrected, lab_ref);
    
    % Compute full ΔE2000 error map
    lab_full_ref = xyz2lab(reshape(Y, [], 3));
    lab_full_corrected = xyz2lab(reshape(corrected_xyz, [], 3));
    error_map = reshape(deltaE2000(lab_full_corrected, lab_full_ref), m, n);
    
    % Display Error Map
    figure;
    imagesc(error_map);
    colormap(jet);
    colorbar;
    clim([0 10]);
    title('\DeltaE2000 Error Map');
    
    % Print error metrics
    disp(['Mean \DeltaE2000 Error: ', num2str(mean(deltaE2000_errors))]);
    disp(['Max \DeltaE2000 Error: ', num2str(max(deltaE2000_errors))]);
end
