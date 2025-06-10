clc; clear; close all;
roof = double(intmax('uint16'));

%% Load spectral data
ill = importdata('../../../data/CIE_D50.txt');
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

illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
sp_tristREF = CMFsIP .* illIP;
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

%% Load and process the reference cube
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

%% ----------------- Cross-Validation -----------------
K = 6;
num_samples = size(xyz_input, 1);
cv = cvpartition(num_samples, 'KFold', K);

lab_input = xyz2lab(xyz_input);
lab_ref = xyz2lab(xyz_ref);
rgb_input = xyz2prophoto(xyz_input ./100, true);
rgb_ref   = xyz2prophoto(xyz_ref ./100, true);
lab_from_rgb_input = rgb2lab(rgb_input, 'ColorSpace', 'prophoto-rgb', 'WhitePoint','d50');
lab_from_rgb_ref   = rgb2lab(rgb_ref,   'ColorSpace', 'prophoto-rgb', 'WhitePoint','d50');

output_folder = '../../../results/error_maps';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
[~, img_name, ~] = fileparts(cubeFile);

mean_errors_xyz = zeros(K, 1);
max_errors_xyz = zeros(K, 1);
mean_errors_lab = zeros(K, 1);
max_errors_lab = zeros(K, 1);
mean_errors_rgb = zeros(K, 1);
max_errors_rgb = zeros(K, 1);
mean_errors_labFromRGB = zeros(K, 1);
max_errors_labFromRGB = zeros(K, 1);

for k = 1:K
    fprintf('Fold %d of %d\n', k, K);
    num_samples = size(xyz_input, 1);
    perm = randperm(num_samples);
    train_idx = perm(1:round(0.8 * num_samples));
    test_idx = perm(round(0.8 * num_samples) + 1:end);

     %% XYZ Regression
    X_poly_xyz_train = poly3_features(xyz_input(train_idx, :));
    coeffs_xyz = pinv(X_poly_xyz_train) * xyz_ref(train_idx, :);
    corrected_xyz = poly3_features(xyz_input) * coeffs_xyz;
    lab_from_xyz_cor = xyz2lab(corrected_xyz);
    [mean_errors_xyz(k), max_errors_xyz(k)] = evaluate_error(xyz2lab(xyz_ref), lab_from_xyz_cor, test_idx, m, n, 'XYZ_', output_folder, img_name + "_fold" + k + "_error_poly3.png");

    %% Lab Regression
    X_poly_lab_train = poly3_features(lab_input(train_idx, :));
    coeffs_lab = pinv(X_poly_lab_train) * lab_ref(train_idx, :);
    corrected_lab = poly3_features(lab_input) * coeffs_lab;
    [mean_errors_lab(k), max_errors_lab(k)] = evaluate_error(lab_ref, corrected_lab, test_idx, m, n, 'Lab_', output_folder, img_name + "_fold" + k + "_error_poly3.png");

    %% RGB Regression
    X_poly_rgb_train = poly3_features(rgb_input(train_idx, :));
    coeffs_rgb = pinv(X_poly_rgb_train) * rgb_ref(train_idx, :);
    corrected_rgb = poly3_features(rgb_input) * coeffs_rgb;
    lab_from_rgb_cor = rgb2lab(corrected_rgb, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');
    [mean_errors_rgb(k), max_errors_rgb(k)] = evaluate_error(rgb2lab(rgb_ref, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50'), lab_from_rgb_cor, test_idx, m, n, 'RGB_', output_folder, img_name + "_fold" + k + "_error_poly3.png");

    %% RGB to Lab Regression
    X_poly_lab_rgb = poly3_features(lab_from_rgb_input(train_idx, :));
    coeffs_lab_rgb = pinv(X_poly_lab_rgb) * lab_from_rgb_ref(train_idx, :);
    corrected_lab_from_rgb = poly3_features(lab_from_rgb_input) * coeffs_lab_rgb;
    [mean_errors_labFromRGB(k), max_errors_labFromRGB(k)] = evaluate_error(lab_from_rgb_ref, corrected_lab_from_rgb, test_idx, m, n, 'LabFromRGB_', output_folder, img_name + "_fold" + k + "_error_poly3.png");
end

fprintf('\nAverage ΔE2000 over %d folds:\n', K);
fprintf('XYZ Regression     : Mean = %.2f, Std = %.2f\n', mean(mean_errors_xyz), std(mean_errors_xyz));
fprintf('Lab Regression     : Mean = %.2f, Std = %.2f\n', mean(mean_errors_lab), std(mean_errors_lab));
fprintf('RGB Regression     : Mean = %.2f, Std = %.2f\n', mean(mean_errors_rgb), std(mean_errors_rgb));
fprintf('Lab-from-RGB Only  : Mean = %.2f, Std = %.2f\n', mean(mean_errors_labFromRGB), std(mean_errors_labFromRGB));


%% 3rd-degree polynomial feature expansion
function X_poly = poly3_features(input_data)
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



%% Delta E Evaluation Function
function [mean_deltaE, max_deltaE] = evaluate_error(ref_lab, corrected_lab, test_idx, m, n, space, output_folder, file_name)
    % Compute Delta E2000 Error
    deltaE2000_errors = deltaE2000(corrected_lab, ref_lab);
    
    % Calculate mean and max Delta E2000 for the test set
    mean_deltaE = mean(deltaE2000_errors(test_idx));
    max_deltaE = max(deltaE2000_errors(test_idx));
    
    % Compute the Delta E2000 error map for the full set (error map for entire dataset)
    error_map = reshape(deltaE2000(corrected_lab, ref_lab), m, n);
    
    % Display the results
    disp([space ' - Mean ΔE2000: ', num2str(mean_deltaE), '; Max ΔE2000: ', num2str(max_deltaE)]);
    % disp([space ' - Max ΔE2000: ', num2str(max_deltaE)]);
    
    % Plot the error map
    figure;
    imagesc(error_map);
    colormap(jet);
    colorbar;
    clim([0 5]);  % Set the color axis limits for clarity
    title([space ' ΔE2000', ' (Mean: ', num2str(mean_deltaE, '%.2f'), ', Max: ', ...
        num2str(max_deltaE, '%.2f'), ')'], Interpreter="none");
    grid off;

    [rows, cols] = ind2sub([m, n], test_idx); % Convert linear indices to (row, col)
    hold on;  % Keep the error map displayed
    plot(cols, rows, 'w.', 'MarkerSize', 10);  % Plot white dots at training locations
    hold off;
    
    % Save the error map figure
    % saveas(gcf, fullfile(output_folder, space + file_name));
end