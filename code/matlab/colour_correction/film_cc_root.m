clear; close all;

%% Load data
inputFile = '/home/oem/eliza/data/xyz_lab_rgb/colorchecker/yoda_halogen_fuji_underexp_balanced_colorchecker.mat';  
refFile   = '/home/oem/eliza/masters-thesis/data/CCpassport_spectra.txt';
ill       = importdata('../../../data/CIE_D50.txt');
CMFs      = importdata('../../../data/CIE2degCMFs_1931.txt');
S         = load(inputFile);

RGB_lin_input = S.patchRGB;

%% Load spectral reference and compute Lab
refData     = readmatrix(refFile, 'NumHeaderLines', 1);  
ref_wl      = refData(:,1);  
ref_spectra = refData(:,2:end);  

valid_idx     = ref_wl >= 380 & ref_wl <= 780;
ref_wl        = ref_wl(valid_idx);
ref_spectra   = ref_spectra(valid_idx, :);

illIP   = interp1(ill(:,1), ill(:,2), ref_wl, 'spline');
CMFsIP  = interp1(CMFs(:,1), CMFs(:,2:4), ref_wl, 'spline');  % x̄ ȳ z̄
sp_tristREF = CMFsIP .* illIP;
k            = 100 / sum(sp_tristREF(:,2));
XYZ_ref      = k * (ref_spectra' * sp_tristREF);
Lab_ref      = xyz2lab_custom(XYZ_ref);  % <--- regression target is Lab!

%% Root-polynomial regression (RGB_lin → Lab_ref)
N = size(RGB_lin_input, 1);
K = 6;
rng(42);
cv = cvpartition(N, 'KFold', K);
mean_de = zeros(K, 1);
max_de  = zeros(K, 1);

lambda = 0.001;

for k = 1:K
    train_idx = training(cv, k);
    test_idx  = test(cv, k);

    X_train = root_poly_features(RGB_lin_input(train_idx, :));
    Y_train = Lab_ref(train_idx, :);

    % Ridge regression
    coeffs = (X_train' * X_train + lambda * eye(size(X_train,2))) \ (X_train' * Y_train);

    % Apply model
    X_all  = root_poly_features(RGB_lin_input);
    Lab_pred = X_all * coeffs;

    deltaE = deltaE2000(Lab_pred, Lab_ref);
    mean_de(k) = mean(deltaE(test_idx));
    max_de(k)  = max(deltaE(test_idx));
end

%% Results
fprintf('\n---- Finlayson Root-Polynomial (RGB_lin → Lab) ----\n');
fprintf('λ = %.4f\n', lambda);
fprintf('ΔE2000 Mean = %.2f ± %.2f | Max = %.2f ± %.2f\n', ...
    mean(mean_de), std(mean_de), ...
    mean(max_de), std(max_de));

%% --- Root-polynomial feature expansion
function X_rootpoly = root_poly_features(RGB)
    R = max(RGB(:,1), 0);
    G = max(RGB(:,2), 0);
    B = max(RGB(:,3), 0);

    X_rootpoly = [ ...
        ones(size(R)), ...
        R, G, B, ...
        sqrt(R), sqrt(G), sqrt(B), ...
        R.*G, R.*B, G.*B, ...
        sqrt(R.*G), sqrt(R.*B), sqrt(G.*B), ...
        R.^2, G.^2, B.^2, ...
        sqrt(R.^2 + G.^2 + B.^2) ];
end
