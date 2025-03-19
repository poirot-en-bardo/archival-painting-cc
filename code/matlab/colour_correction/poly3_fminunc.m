clc; clear; close all;
roof = double(intmax('uint16'));
histoFACT = 200; % Histogram normalization factor

%% Load spectral data
ill = importdata('../../../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../../../data/CIE2degCMFs_1931.txt');
CMFs_2006 = importdata('../../../data/CIE2degCMFs_2006.txt');
CMFs = CMFs_1931;
rng(10);

%% Select the reference and the cube to be corrected
cube_folder = '../../../data/colorChecker_SG/colorChecker_SG_Elias';
cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_fuji-frame4.hdr";
refFile  = "../../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

%% Load and process the cube to be corrected (input cube)
hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd] = size(inCUBE);
lincube = reshape(inCUBE, [], bd);

% Interpolate illuminant and CMFs to captured wavelengths
illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
sp_tristREF = CMFsIP .* illIP;
% Compute XYZ values for the input cube (normalized by Y)
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

%%% Load and process the reference cube
hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
[m, n, bd] = size(refCUBE);
lincube_ref = reshape(refCUBE, [], bd);

illIP = interp1(ill(:,1), ill(:,2), bands_ref, 'spline');
CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands_ref, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,3), bands_ref, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,4), bands_ref, 'spline')];
sp_tristREF = CMFsIP .* illIP;
% Compute XYZ values for the reference cube
xyz_ref = (lincube_ref * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

%% Train-test split (80% training, 20% testing)
num_samples = size(xyz_input, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx  = perm(round(0.8*num_samples)+1:end);

xyz_input_train = xyz_input(train_idx, :);
xyz_ref_train   = xyz_ref(train_idx, :);

%% ----------------- XYZ-Based Regression via fminunc -----------------
% Compute 3rd-degree polynomial features from the training set (XYZ mode)
X_poly_xyz_train = poly3_features(xyz_input_train);  % N_train x 19

% Define objective function that computes total ΔE2000 error (XYZ mode)
objFun_xyz = @(p) deltaE2000_loss(p, X_poly_xyz_train, xyz_ref_train, 'xyz');

% Initialize parameters (vectorized regression coefficients for a 19x3 matrix)
initialParams_xyz = rand(size(X_poly_xyz_train,2)*3, 1);
options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton', 'MaxFunctionEvaluations', 10000);
optParams_xyz = fminunc(objFun_xyz, initialParams_xyz, options);

% Reshape optimized parameters to obtain regression matrix (19 x 3)
coeffs_xyz = reshape(optParams_xyz, size(X_poly_xyz_train,2), 3);

% Apply the optimized XYZ regression to the full dataset
X_poly_xyz = poly3_features(xyz_input);  % N x 19 matrix
corrected_xyz = X_poly_xyz * coeffs_xyz;

%% ----------------- Lab-Based Regression via fminunc -----------------

% Convert XYZ to Lab
lab_input = xyz2lab(xyz_input);
lab_ref = xyz2lab(xyz_ref);

lab_input_train = lab_input(train_idx, :);
lab_ref_train   = lab_ref(train_idx, :);

% Compute 3rd-degree polynomial features for Lab training data
X_poly_lab_train = poly3_features(lab_input_train);  % N_train x 19

% Define objective function for Lab mode
objFun_lab = @(p) deltaE2000_loss(p, X_poly_lab_train, lab_ref_train, 'lab');

initialParams_lab = rand(size(X_poly_lab_train,2)*3, 1);
optParams_lab = fminunc(objFun_lab, initialParams_lab, options);
coeffs_lab = reshape(optParams_lab, size(X_poly_lab_train,2), 3);

% Apply the optimized Lab regression to the full dataset
X_poly_lab = poly3_features(lab_input);  % N x 19 matrix
corrected_lab = X_poly_lab * coeffs_lab;

%% ----------------- RGB-Based Regression via fminunc -----------------
% Convert input and reference XYZ to RGB (using prophoto-rgb)
rgb_input = xyz2rgb(xyz_input, 'ColorSpace', 'prophoto-rgb');
rgb_ref   = xyz2rgb(xyz_ref,   'ColorSpace', 'prophoto-rgb');

rgb_input_train = rgb_input(train_idx, :);
rgb_ref_train   = rgb_ref(train_idx, :);

% Compute 3rd-degree polynomial features for RGB training data
X_poly_rgb_train = poly3_features(rgb_input_train);  % N_train x 19

% Define objective function for RGB mode
objFun_rgb = @(p) deltaE2000_loss(p, X_poly_rgb_train, rgb_ref_train, 'rgb');

initialParams_rgb = rand(size(X_poly_rgb_train,2)*3, 1);
optParams_rgb = fminunc(objFun_rgb, initialParams_rgb, options);
coeffs_rgb = reshape(optParams_rgb, size(X_poly_rgb_train,2), 3);

% Apply the optimized RGB regression to the full dataset
X_poly_rgb = poly3_features(rgb_input);  % N x 19 matrix
corrected_rgb = X_poly_rgb * coeffs_rgb;

%% ----------------- Evaluate & Display Results -----------------
% Convert corrected outputs to Lab for error evaluation
lab_from_xyz_ref = xyz2lab(xyz_ref);
lab_from_xyz_cor = xyz2lab(corrected_xyz);

lab_from_rgb_ref = rgb2lab(rgb_ref, 'ColorSpace', 'prophoto-rgb');
lab_from_rgb_cor = rgb2lab(corrected_rgb, 'ColorSpace', 'prophoto-rgb');

output_folder = '../../../results/error_maps';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

[~, img_name, ~] = fileparts(cubeFile);
evaluate_error(lab_from_xyz_ref, lab_from_xyz_cor, test_idx, m, n, 'XYZ_', output_folder, strcat(img_name, "_poly3_fminunc.png"));
evaluate_error(lab_ref, corrected_lab, test_idx, m, n, 'Lab_', output_folder, strcat(img_name, "_poly3_fminunc.png"));
evaluate_error(lab_from_rgb_ref, lab_from_rgb_cor, test_idx, m, n, 'RGB_', output_folder, strcat(img_name, "_poly3_fminunc.png"));

%% ===================== Local Helper Functions =====================

function X_poly = poly3_features(input_data)
    % Constructs a 3rd-degree polynomial expansion (without constant term)
    % for a 3-column input. For input columns a, b, c, this returns:
    % [ a, b, c, a.^2, b.^2, c.^2, a.*b, a.*c, b.*c, ...
    %   a.^3, b.^3, c.^3, a.^2.*b, a.^2.*c, b.^2.*a, b.^2.*c, c.^2.*a, c.^2.*b, a.*b.*c ]
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);
    
    % Degree 1
    feat1 = a;
    feat2 = b;
    feat3 = c;
    
    % Degree 2
    feat4 = a.^2;
    feat5 = b.^2;
    feat6 = c.^2;
    feat7 = a.*b;
    feat8 = a.*c;
    feat9 = b.*c;
    
    % Degree 3
    feat10 = a.^3;
    feat11 = b.^3;
    feat12 = c.^3;
    feat13 = a.^2 .* b;
    feat14 = a.^2 .* c;
    feat15 = b.^2 .* a;
    feat16 = b.^2 .* c;
    feat17 = c.^2 .* a;
    feat18 = c.^2 .* b;
    feat19 = a .* b .* c;
    
    X_poly = [feat1, feat2, feat3, feat4, feat5, feat6, feat7, feat8, feat9, ...
              feat10, feat11, feat12, feat13, feat14, feat15, feat16, feat17, feat18, feat19];
end

function loss = deltaE2000_loss(p, X_poly, target, mode)
    % p: vectorized regression coefficients (numFeatures*3 x 1)
    % X_poly: N x numFeatures feature matrix
    % target: N x 3 target values (XYZ, RGB, or Lab)
    % mode: 'xyz', 'rgb', or 'lab'
    
    coeffs = reshape(p, size(X_poly,2), 3);
    pred = X_poly * coeffs; % Predicted colors (N x 3)
    
    if strcmpi(mode, 'rgb')
        pred = min(max(pred,0),1); % clip RGB predictions to [0,1]
    end
    
    % Convert predicted and target values to Lab
    switch lower(mode)
        case 'xyz'
            lab_pred = xyz2lab(pred);
            lab_target = xyz2lab(target);
        case 'rgb'
            lab_pred = rgb2lab(pred, 'ColorSpace', 'prophoto-rgb');
            lab_target = rgb2lab(target, 'ColorSpace', 'prophoto-rgb');
        case 'lab'
            lab_pred = pred; % Already in Lab space
            lab_target = target;
        otherwise
            error('Unknown mode: %s', mode);
    end
    
    % Use the provided deltaE2000 function
    DE = deltaE2000(lab_target, lab_pred, [1, 1, 1]);
    loss = sum(DE);
end
