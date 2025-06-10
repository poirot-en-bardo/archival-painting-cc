%% Full Script: Pseudoinverse for 4th Degree Polynomial Regression Minimizing ΔE2000
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

%% ----------------- XYZ-Based Regression -----------------
% Compute 4th-degree polynomial features from the training set (XYZ mode)
X_poly_xyz_train = poly4_features(xyz_input_train);  

% Compute coefficients using pseudoinverse
coeffs_xyz = pinv(X_poly_xyz_train) * xyz_ref_train;  

% Apply the regression to the full dataset
X_poly_xyz = poly4_features(xyz_input);  
corrected_xyz = X_poly_xyz * coeffs_xyz;

%% Convert to Lab and apply 4th order polynomial correction
lab_input = xyz2lab(xyz_input);
lab_ref = xyz2lab(xyz_ref);

lab_input_train = lab_input(train_idx, :);
lab_ref_train = lab_ref(train_idx, :);

X_poly_lab_train =  poly4_features(lab_input_train);
coeffs_lab = pinv(X_poly_lab_train) * lab_ref_train;

X_poly_lab = poly4_features(lab_input);
corrected_lab = X_poly_lab * coeffs_lab;

%% ----------------- RGB-Based Regression  -----------------
% Convert input and reference XYZ to RGB (using prophoto-rgb)
rgb_input = xyz2rgb(xyz_input ./100, 'ColorSpace', 'prophoto-rgb');
rgb_ref   = xyz2rgb(xyz_ref ./100,   'ColorSpace', 'prophoto-rgb');

rgb_input_train = rgb_input(train_idx, :);
rgb_ref_train   = rgb_ref(train_idx, :);

% Compute 4th-degree polynomial features for RGB training data
X_poly_rgb_train = poly4_features(rgb_input_train);  

% Compute coefficients using pseudoinverse
coeffs_rgb = pinv(X_poly_rgb_train) * rgb_ref_train;  

% Apply the regression to the full dataset
X_poly_rgb = poly4_features(rgb_input);  
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

evaluate_error(lab_from_xyz_ref, lab_from_xyz_cor, test_idx, m, n, 'XYZ_', output_folder, img_name + "_error_poly4.png");
evaluate_error(lab_ref, corrected_lab, test_idx, m, n, 'Lab_', output_folder, img_name + "_error_poly4.png");
evaluate_error(lab_from_rgb_ref, lab_from_rgb_cor, test_idx, m, n, 'RGB_', output_folder, img_name + "_error_poly4.png");

%% function for creating the 4th-degree polynomial expression
function X_poly = poly4_features(input_data)
    % Constructs a 4th-degree polynomial expansion (without constant term)
    % for a 3-column input. For input columns a, b, c, this returns 34 terms:
    %
    % Degree 1 (3 terms):
    %   a, b, c
    %
    % Degree 2 (6 terms):
    %   a.^2, b.^2, c.^2, a.*b, a.*c, b.*c
    %
    % Degree 3 (10 terms):
    %   a.^3, b.^3, c.^3, a.^2.*b, a.^2.*c, a.*b.^2, a.*c.^2, b.^2.*c, b.*c.^2, a.*b.*c
    %
    % Degree 4 (15 terms):
    %   a.^4, b.^4, c.^4, a.^3.*b, a.^3.*c, a.*b.^3, a.*c.^3, b.^3.*c, b.*c.^3, ...
    %   a.^2.*b.^2, a.^2.*c.^2, b.^2.*c.^2, a.^2.*b.*c, a.*b.^2.*c, a.*b.*c.^2
    %
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);

    constant = ones(size(a));
    
    % Degree 1
    feat1  = a;
    feat2  = b;
    feat3  = c;
    
    % Degree 2
    feat4  = a.^2;
    feat5  = b.^2;
    feat6  = c.^2;
    feat7  = a.*b;
    feat8  = a.*c;
    feat9  = b.*c;
    
    % Degree 3
    feat10  = a.^3;
    feat11  = b.^3;
    feat12  = c.^3;
    feat13  = a.^2 .* b;
    feat14  = a.^2 .* c;
    feat15  = a .* b.^2;
    feat16  = a .* c.^2;
    feat17  = b.^2 .* c;
    feat18  = b .* c.^2;
    feat19  = a .* b .* c;
    
    % Degree 4
    feat20 = a.^4;
    feat21 = b.^4;
    feat22 = c.^4;
    feat23 = a.^3 .* b;
    feat24 = a.^3 .* c;
    feat25 = a .* b.^3;
    feat26 = a .* c.^3;
    feat27 = b.^3 .* c;
    feat28 = b .* c.^3;
    feat29 = a.^2 .* b.^2;
    feat30 = a.^2 .* c.^2;
    feat31 = b.^2 .* c.^2;
    feat32 = a.^2 .* b .* c;
    feat33 = a .* b.^2 .* c;
    feat34 = a .* b .* c.^2;
    
    X_poly = [constant, feat1, feat2, feat3, feat4, feat5, feat6, feat7, feat8, feat9, ...
              feat10, feat11, feat12, feat13, feat14, feat15, feat16, feat17, feat18, feat19, ...
              feat20, feat21, feat22, feat23, feat24, feat25, feat26, feat27, feat28, ...
              feat29, feat30, feat31, feat32, feat33, feat34];
end
