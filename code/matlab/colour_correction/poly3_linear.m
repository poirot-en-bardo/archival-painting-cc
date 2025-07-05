%% Full Script: Pseudoinverse for 3rd Degree Polynomial Regression Minimizing ΔE2000
clc; clear; close all;
roof = double(intmax('uint16'));

%% Load spectral data
ill = importdata('../../../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../../../data/CIE2degCMFs_1931.txt');
CMFs_2006 = importdata('../../../data/CIE2degCMFs_2006.txt');
CMFs = CMFs_1931;
rng(10);

%% Select the reference and the cube to be corrected
% cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_fuji-frame4.hdr";


% Set folders and get list of cube files
cube_folder = '../../../data/colorChecker_SG/cubes';
refFile  = "../../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";
output_folder = '../../../results/error_maps/poly3';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

% List all .hdr files in the folder (skip reference cube)
cube_files_struct = dir(fullfile(cube_folder, '*.hdr'));

for k = 1:length(cube_files_struct)
    cubeFile = fullfile(cube_folder, cube_files_struct(k).name);

    % Skip reference file if in same folder
    if strcmpi(cubeFile, refFile)
        continue
    end

    % Load and process the cube to be corrected (input cube)
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
    xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

    % Load and process the reference cube
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
    xyz_ref = (lincube_ref * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

    % Train-test split (80% training, 20% testing)
    num_samples = size(xyz_input, 1);
    perm = randperm(num_samples);
    train_idx = perm(1:round(0.8*num_samples));
    test_idx  = perm(round(0.8*num_samples)+1:end);

    xyz_input_train = xyz_input(train_idx, :);
    xyz_ref_train   = xyz_ref(train_idx, :);

    % XYZ Regression
    X_poly_xyz_train = poly3_features(xyz_input_train);  
    coeffs_xyz = pinv(X_poly_xyz_train) * xyz_ref_train;  
    X_poly_xyz = poly3_features(xyz_input);  
    corrected_xyz = X_poly_xyz * coeffs_xyz;

    % Lab Regression
    lab_input = xyz2lab(xyz_input);
    lab_ref = xyz2lab(xyz_ref);

    lab_input_train = lab_input(train_idx, :);
    lab_ref_train = lab_ref(train_idx, :);

    X_poly_lab_train = poly3_features(lab_input_train);
    coeffs_lab = pinv(X_poly_lab_train) * lab_ref_train;
    X_poly_lab = poly3_features(lab_input);
    corrected_lab = X_poly_lab * coeffs_lab;

    % RGB Regression (linear RGB)
    rgb_input = xyz2rgb(xyz_input ./100, 'ColorSpace', 'linear-rgb', 'WhitePoint', 'd50');
    rgb_ref   = xyz2rgb(xyz_ref ./100,   'ColorSpace', 'linear-rgb', 'WhitePoint', 'd50');


    rgb_input_train = rgb_input(train_idx, :);
    rgb_ref_train   = rgb_ref(train_idx, :);

    X_poly_rgb_train = poly3_features(rgb_input_train);
    coeffs_rgb = pinv(X_poly_rgb_train) * rgb_ref_train;
    X_poly_rgb = poly3_features(rgb_input);
    corrected_rgb = X_poly_rgb * coeffs_rgb;

    % Evaluation
    corrected_rgb_srgb = linear2srgb(corrected_rgb);
    rgb_ref_srgb_eval  = linear2srgb(rgb_ref);

    lab_from_rgb_ref = rgb2lab(rgb_ref_srgb_eval, 'ColorSpace', 'srgb');
    lab_from_rgb_cor = rgb2lab(corrected_rgb_srgb, 'ColorSpace', 'srgb');
    lab_from_xyz_ref = xyz2lab(xyz_ref);
    lab_from_xyz_cor = xyz2lab(corrected_xyz);

    [~, img_name, ~] = fileparts(cubeFile);

    evaluate_error(lab_from_xyz_ref, lab_from_xyz_cor, test_idx, m, n, 'XYZ_', output_folder, img_name + "_error_poly3.png");
    evaluate_error(lab_ref, corrected_lab, test_idx, m, n, 'Lab_', output_folder, img_name + "_error_poly3.png");
    evaluate_error(lab_from_rgb_ref, lab_from_rgb_cor, test_idx, m, n, 'RGB_', output_folder, img_name + "_error_poly3.png");
end

%% --------- Helper Functions ----------

function rgb_linear = srgb2linear(rgb)
    % sRGB to linear RGB conversion
    % Assumes input in [0,1]
    threshold = 0.04045;
    rgb_linear = zeros(size(rgb));
    mask = rgb <= threshold;
    rgb_linear(mask) = rgb(mask) / 12.92;
    rgb_linear(~mask) = ((rgb(~mask) + 0.055) / 1.055) .^ 2.4;
end

function rgb_srgb = linear2srgb(rgb)
    % linear RGB to sRGB conversion
    % Assumes input in [0,1]
    threshold = 0.0031308;
    rgb_srgb = zeros(size(rgb));
    mask = rgb <= threshold;
    rgb_srgb(mask) = 12.92 * rgb(mask);
    rgb_srgb(~mask) = 1.055 * rgb(~mask).^(1/2.4) - 0.055;
end

function X_poly = poly3_features(input_data)
    % 3rd-degree polynomial expansion (WITH constant term)
    % For input [a b c], returns [1, a, b, c, ..., a*b*c]
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);

    X_poly = [ ...
        ones(size(a)), ...         
        a, b, c, ...
        a.^2, b.^2, c.^2, ...
        a.*b, a.*c, b.*c, ...
        a.^3, b.^3, c.^3, ...
        a.^2.*b, a.^2.*c, ...
        b.^2.*a, b.^2.*c, ...
        c.^2.*a, c.^2.*b, ...
        a.*b.*c];
end
