%% Define parameters
clc; clear; close all;
roof = double(intmax('uint16'));

histoFACT = 200; % Histogram normalization factor

% Load spectral data
ill = importdata('../../../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../../../data/CIE2degCMFs_1931.txt');
CMFs_2006 = importdata('../../../data/CIE2degCMFs_2006.txt');
CMFs = CMFs_1931;
rng(10);

%% Select the reference and the cube to be corrected
cube_folder = '../../../data/colorChecker_SG/colorChecker_SG_Elias';

% [cubeFileName, cubePath] = uigetfile(fullfile(cube_folder, '*.hdr'), 'Select Cube to be Corrected');
% cubeFile = fullfile(cubePath, cubeFileName);
cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_120f-velvia-f8.hdr";

% [refFileName, refPath] = uigetfile(fullfile(cube_folder, '*.hdr'), 'Select Reference Cube');
% refFile = fullfile(refPath, refFileName);
refFile = "../../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

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


%% Train-test split (80% training, 20% testing)
num_samples = size(xyz_input, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx = perm(round(0.8*num_samples)+1:end);

xyz_input_train = xyz_input(train_idx, :);
xyz_ref_train = xyz_ref(train_idx, :);

%% XYZ Root-Polynomial Regression Model for Color Correction
% Define 2nd degree root-polynomial features
X_poly_xyz = [xyz_input, sqrt(xyz_input(:,1).*xyz_input(:,2)), ...
     sqrt(xyz_input(:,2).*xyz_input(:,3)), ... 
     sqrt(xyz_input(:,1).*xyz_input(:,3))]; 

% Assume xyz_input is an N x 3 matrix: [X, Y, Z]
X = xyz_input(:,1);
Y = xyz_input(:,2);
Z = xyz_input(:,3);

% Linear terms (degree 1)
feat1 = X;
feat2 = Y;
feat3 = Z;

% 2nd degree terms: square roots of pairwise products
feat4 = sqrt(X .* Y);
feat5 = sqrt(X .* Z);
feat6 = sqrt(Y .* Z);

% 3rd degree terms: cube roots of selected monomials
feat7  = (X .* Y.^2).^(1/3);
feat8  = (X .* Z.^2).^(1/3);
feat9  = (Y .* X.^2).^(1/3);
feat10 = (Y .* Z.^2).^(1/3);
feat11 = (Z .* X.^2).^(1/3);
feat12 = (Z .* Y.^2).^(1/3);
feat13 = (X .* Y .* Z).^(1/3);

% 4th degree terms: 4th roots of selected monomials
feat14 = (X.^3 .* Y).^(1/4);
feat15 = (X.^3 .* Z).^(1/4);
feat16 = (Y.^3 .* X).^(1/4);
feat17 = (Y.^3 .* Z).^(1/4);
feat18 = (Z.^3 .* X).^(1/4);
feat19 = (Z.^3 .* Y).^(1/4);
feat20 = (X.^2 .* Y.^2).^(1/4);
feat21 = (X.^2 .* Z.^2).^(1/4);
feat22 = (Y.^2 .* Z.^2).^(1/4);

% Combine all features into one matrix (N x 22)
X_poly_xyz = [feat1, feat2, feat3, ...
                   feat4, feat5, feat6, ...
                   feat7, feat8, feat9, feat10, feat11, feat12, feat13, ...
                   feat14, feat15, feat16, feat17, feat18, feat19, ...
                   feat20, feat21, feat22];


X_poly_xyz_train = X_poly_xyz(train_idx, :);

% Fit root-polynomial regression model (using least squares)
% coeffs_xyz = (X_poly_xyz_train' * X_poly_xyz_train) \ (X_poly_xyz_train' * xyz_ref_train);
% coeffs_xyz = (xyz_ref_train' * X_poly_xyz_train * pinv(X_poly_xyz_train' * X_poly_xyz_train))'; % Moore-Penrose Inverse
coeffs_xyz = pinv(X_poly_xyz_train) * xyz_ref_train;


% Apply correction to the entire input cube
corrected_xyz = X_poly_xyz * coeffs_xyz;


%% RGB conversion and regression
rgb_input = xyz2rgb(xyz_input, 'ColorSpace','prophoto-rgb');
rgb_ref = xyz2rgb(xyz_ref, 'ColorSpace','prophoto-rgb');

rgb_input_train = rgb_input(train_idx, :);
rgb_ref_train = rgb_ref(train_idx, :);
% 
% X_poly_rgb = [rgb_input, sqrt(rgb_input(:,1).*rgb_input(:,2)), ...
%      sqrt(rgb_input(:,2).*rgb_input(:,3)), ... 
%      sqrt(rgb_input(:,1).*rgb_input(:,3))];

% Assume xyz_input is an N x 3 matrix: [r, g, b]
r = rgb_input(:,1);
g = rgb_input(:,2);
b = rgb_input(:,3);

% Linear terms (degree 1)
feat1 = r;
feat2 = g;
feat3 = b;

% 2nd degree terms: square roots of pairwise products
feat4 = sqrt(r .* g);
feat5 = sqrt(r .* b);
feat6 = sqrt(g .* b);

% 3rd degree terms: cube roots of selected monomials
feat7  = (r .* g.^2).^(1/3);
feat8  = (r .* b.^2).^(1/3);
feat9  = (g .* r.^2).^(1/3);
feat10 = (g .* b.^2).^(1/3);
feat11 = (b .* r.^2).^(1/3);
feat12 = (b .* g.^2).^(1/3);
feat13 = (r .* g .* b).^(1/3);

% 4th degree terms: 4th roots of selected monomials
feat14 = (r.^3 .* g).^(1/4);
feat15 = (r.^3 .* b).^(1/4);
feat16 = (g.^3 .* r).^(1/4);
feat17 = (g.^3 .* b).^(1/4);
feat18 = (b.^3 .* r).^(1/4);
feat19 = (b.^3 .* g).^(1/4);
feat20 = (r.^2 .* g.^2).^(1/4);
feat21 = (r.^2 .* b.^2).^(1/4);
feat22 = (g.^2 .* b.^2).^(1/4);

% Combine all features into one matrix (N x 22)
X_poly_rgb = [feat1, feat2, feat3, ...
               feat4, feat5, feat6, ...
               feat7, feat8, feat9, feat10, feat11, feat12, feat13, ...
               feat14, feat15, feat16, feat17, feat18, feat19, ...
               feat20, feat21, feat22];


X_poly_rgb_train = X_poly_rgb(train_idx, :);
% coeffs_rgb = (X_poly_rgb_train' * X_poly_rgb_train) \ (X_poly_rgb_train' * rgb_ref_train);
coeffs_rgb = pinv(X_poly_rgb_train) * rgb_ref_train;


corrected_rgb = X_poly_rgb * coeffs_rgb;



%% Display results for XYZ, Lab, and RGB
output_folder = '../../../results/error_maps';  % Output folder to save the error maps
mkdir(output_folder);  % Create the output folder if it doesn't exist

lab_from_xyz_ref = xyz2lab(xyz_ref);
lab_from_xyz_cor = xyz2lab(corrected_xyz);

lab_from_rgb_ref = rgb2lab(rgb_ref, 'ColorSpace','prophoto-rgb');
lab_from_rgb_cor = rgb2lab(corrected_rgb, 'ColorSpace','prophoto-rgb');

[~, img_name, ~] = fileparts(cubeFile);

evaluate_error(lab_from_xyz_ref, lab_from_xyz_cor, test_idx, m, n, 'XYZ_', output_folder, img_name + "_error_root.png");
evaluate_error(lab_from_rgb_ref, lab_from_rgb_cor, test_idx, m, n, 'RGB_', output_folder, img_name + "_error_root.png");

