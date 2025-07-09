clear; clc; close all;

%% Load input and reference patch data
inputFile     = '/home/oem/eliza/data/xyz_lab_rgb/colorchecker/yoda_halogen_fuji_underexp_balanced_colorchecker.mat';
referenceFile = '/home/oem/eliza/data/xyz_lab_rgb/reference/xrite_cc_reference_official.mat';

input     = load(inputFile);
reference = load(referenceFile);

RGB_input = input.patchRGB_lin;     % gamma-encoded (required)
XYZ_ref   = reference.XYZ_ref;  % target

%% Split data
N = size(RGB_input, 1);
K = 6;
rng(42);
cv = cvpartition(N, 'KFold', K);

mean_de = zeros(K,1);
max_de  = zeros(K,1);

for k = 1:K
    train_idx = training(cv, k);
    test_idx  = test(cv, k);

    X_train = RGB_input(train_idx, :);
    Y_train = XYZ_ref(train_idx, :);

    X_test  = RGB_input(test_idx, :);
    Y_test  = XYZ_ref(test_idx, :);

    %% Define NN: 3-79-36-3
    layers = [
        featureInputLayer(3)
        fullyConnectedLayer(79)
        reluLayer
        fullyConnectedLayer(36)
        reluLayer
        fullyConnectedLayer(3)
        regressionLayer
    ];

    options = trainingOptions('adam', ...
        'MaxEpochs', 500, ...
        'MiniBatchSize', 8, ...
        'InitialLearnRate', 0.001, ...
        'Shuffle', 'every-epoch', ...
        'ValidationData', {X_train, Y_train}, ...
        'ValidationFrequency', 10, ...
        'Verbose', false, ...
        'Plots', 'none', ...
        'ExecutionEnvironment','cpu', ...
        'ValidationPatience', 100);

    net = trainNetwork(X_train, Y_train, layers, options);

    %% Predict
    Y_pred = predict(net, RGB_input);  % full set
    Lab_pred = xyz2lab_custom(Y_pred);
    Lab_ref  = xyz2lab_custom(XYZ_ref);

    deltaE = deltaE2000(Lab_pred, Lab_ref);
    mean_de(k) = mean(deltaE(test_idx));
    max_de(k)  = max(deltaE(test_idx));
end

%% Report
fprintf('\n---- Macdonald Neural Network (RGB → XYZ → Lab) ----\n');
fprintf('ΔE2000 Mean = %.2f ± %.2f | Max = %.2f ± %.2f\n', ...
    mean(mean_de), std(mean_de), ...
    mean(max_de),  std(max_de));
