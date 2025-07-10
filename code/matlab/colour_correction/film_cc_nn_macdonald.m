clear; clc; close all;

%% Load input and reference patch data
inputFile     = '/home/oem/eliza/data/xyz_lab_rgb/colorchecker/yoda_halogen_fuji_underexp_balanced_colorchecker.mat';
referenceFile = '/home/oem/eliza/data/xyz_lab_rgb/reference/xrite_cc_reference_official.mat';

rng(42);

input     = load(inputFile);
reference = load(referenceFile);

RGB_input = input.patchRGB_lin;     % gamma-encoded RGB
XYZ_ref   = reference.XYZ_ref;      % target XYZ

%% Normalize inputs
RGB_input = RGB_input ./ max(RGB_input(:));
XYZ_ref   = XYZ_ref ./ max(XYZ_ref(:));

%% Convert to dlarray
X = dlarray(single(RGB_input'), 'CB');  % (features x batch)
Y = dlarray(single(XYZ_ref'), 'CB');

%% Define network
layers = [
    featureInputLayer(3, 'Normalization', 'none', 'Name', 'input')
    fullyConnectedLayer(79, 'Name', 'fc1')
    eluLayer('Name', 'elu1')
    fullyConnectedLayer(36, 'Name', 'fc2')
    eluLayer('Name', 'elu2')
    fullyConnectedLayer(3, 'Name', 'fc3')
];

lgraph = layerGraph(layers);
net = dlnetwork(lgraph);

%% Training options
numEpochs = 1000;
miniBatchSize = 8;
learnRate = 0.001;
numSamples = size(X, 2);
numIterationsPerEpoch = floor(numSamples / miniBatchSize);

trailingAvg = [];
trailingAvgSq = [];

%% Training loop
for epoch = 1:numEpochs
    idx = randperm(numSamples);
    for i = 1:numIterationsPerEpoch
        batchIdx = idx((i-1)*miniBatchSize+1:i*miniBatchSize);
        XBatch = X(:, batchIdx);
        YBatch = Y(:, batchIdx);

        [loss, gradients] = dlfeval(@modelLoss, net, XBatch, YBatch);
        [net, trailingAvg, trailingAvgSq] = adamupdate(net, gradients, ...
            trailingAvg, trailingAvgSq, epoch, learnRate);
    end

    if mod(epoch, 50) == 0
        fprintf('Epoch %d, Loss = %.4f\n', epoch, extractdata(loss));
    end
end

%% Predict and evaluate
YPred = predict(net, X);

Lab_pred = xyz2lab_dl(YPred);
Lab_true = xyz2lab_dl(Y);


Lab_pred = extractdata(Lab_pred)';
Lab_true = extractdata(Lab_true)';
deltaE = deltaE2000(Lab_pred, Lab_true);
fprintf('\nΔE2000 Mean = %.2f | Max = %.2f\n', mean(deltaE), max(deltaE));

%% Loss function using Lab difference (ΔE approximation)
function [loss, gradients] = modelLoss(net, X, Y)
    % Forward pass
    YPred = forward(net, X);
    YPred = max(YPred, 0);  % Clamp to avoid negative XYZ

    % Convert XYZ to Lab (approximate, differentiable)
    Lab_pred = xyz2lab_dl(YPred);
    Lab_true = xyz2lab_dl(Y);

    % Compute Euclidean distance in Lab space (ΔE approximation)
    % deltaE = sqrt(sum((Lab_pred - Lab_true).^2, 1));  % 1 = batch dim
    % loss = mean(deltaE);

     loss = mean(sum((Lab_pred - Lab_true).^2, 1));


    % Compute gradients
    gradients = dlgradient(loss, net.Learnables);
end

function Lab = xyz2lab_dl(XYZ)
    % XYZ: dlarray of size (3, batch)
    % Normalize by D65 white point
    whitePoint = dlarray([0.95047; 1.00000; 1.08883]);  % D65
    XYZ = XYZ ./ whitePoint;

    % Constants
    epsilon = 216/24389;
    kappa = 24389/27;

    % f(t) function
    f = @(t) (t > epsilon) .* t.^(1/3) + (t <= epsilon) .* ((kappa * t + 16) / 116);

    fx = f(XYZ(1,:));
    fy = f(XYZ(2,:));
    fz = f(XYZ(3,:));

    L = 116 * fy - 16;
    a = 500 * (fx - fy);
    b = 200 * (fy - fz);

    Lab = [L; a; b];
end
