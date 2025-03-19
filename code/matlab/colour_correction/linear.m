%% Define parameters
clc; clear; close all;
roof = double(intmax('uint16')); % Histogram normalization factor

% Load spectral data
ill = importdata('../../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../../data/CIE2degCMFs_1931.txt');
CMFs = CMFs_1931;
rng(10);

%% Select the reference and the cube to be corrected

cubeFile = "../../data/colorChecker_SG/cubes/cubeCC_fuji-frame4.hdr";
refFile = "../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

%% Load and process the cubes
hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd] = size(inCUBE);
lincube = reshape(inCUBE, [], bd);

hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
[m, n, bd] = size(refCUBE);
lincube_ref = reshape(refCUBE, [], bd);

% Interpolate illuminant and CMFs
illIP = interp1(ill(:,1),ill(:,2),bands,'spline');
CMFsIP = [interp1(CMFs(:,1),CMFs(:,2),bands,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,3),bands,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,4),bands,'spline')];
sp_tristREF = CMFsIP .* illIP;

% Compute XYZ values for the input and reference cubes
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2),1);

illIP = interp1(ill(:,1),ill(:,2),bands_ref,'spline');
CMFsIP = [interp1(CMFs(:,1),CMFs(:,2),bands_ref,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,3),bands_ref,'spline'), ...
          interp1(CMFs(:,1),CMFs(:,4),bands_ref,'spline')];
sp_tristREF = CMFsIP .* illIP;

xyz_ref = (lincube_ref * sp_tristREF) ./ sum(sp_tristREF(:,2),1);

%% Train-test split (80% training, 20% testing)
num_samples = size(xyz_input, 1);
perm = randperm(num_samples);
train_idx = perm(1:round(0.8*num_samples));
test_idx = perm(round(0.8*num_samples)+1:end);

xyz_input_train = xyz_input(train_idx, :);
xyz_ref_train = xyz_ref(train_idx, :);

%% Polynomial Regression Model for XYZ
% Fit polynomial regression for XYZ
coeffs_xyz = (xyz_input_train' * xyz_input_train) \ (xyz_input_train' * xyz_ref_train);
corrected_xyz = xyz_input * coeffs_xyz;

%% Convert to Lab and apply linear correction
lab_input = xyz2lab(xyz_input);
lab_ref = xyz2lab(xyz_ref);

lab_input_train = lab_input(train_idx, :);
lab_ref_train = lab_ref(train_idx, :);

coeffs_lab = (lab_input_train' * lab_input_train) \ (lab_input_train' * lab_ref_train);
corrected_lab = lab_input * coeffs_lab;

%% Linear correction in RGB space
rgb_input = xyz2rgb(xyz_input, 'ColorSpace', 'prophoto-rgb');
rgb_ref = xyz2rgb(xyz_ref, 'ColorSpace', 'prophoto-rgb');

rgb_input_train = rgb_input(train_idx, :);
rgb_ref_train = rgb_ref(train_idx, :);


coeffs_rgb = (rgb_input_train' * rgb_input_train) \ (rgb_input_train' * rgb_ref_train);
corrected_rgb = rgb_input * coeffs_rgb;

%% Display results for XYZ, Lab, and RGB
output_folder = '../results/error_maps';  % Output folder to save the error maps
mkdir(output_folder);  % Create the output folder if it doesn't exist

lab_from_xyz_ref = xyz2lab(xyz_ref);
lab_from_xyz_cor = xyz2lab(corrected_xyz);

lab_from_rgb_ref = rgb2lab(rgb_ref, 'ColorSpace','prophoto-rgb');
lab_from_rgb_cor = rgb2lab(corrected_rgb, 'ColorSpace','prophoto-rgb');

evaluate_error(lab_from_xyz_ref, lab_from_xyz_cor, test_idx, m, n, 'XYZ', output_folder, '_error_linear.png');
evaluate_error(lab_ref, corrected_lab, test_idx, m, n, 'Lab', output_folder, '_error_linear.png');
evaluate_error(lab_from_rgb_ref, lab_from_rgb_cor, test_idx, m, n, 'RGB', output_folder, '_error_linear.png');

%% Delta E Evaluation Function
function evaluate_error(ref_lab, corrected_lab, test_idx, m, n, space, output_folder, file_name)
    % Compute Delta E2000 Error
    deltaE2000_errors = deltaE2000(corrected_lab, ref_lab);
    
    % Calculate mean and max Delta E2000 for the test set
    mean_deltaE = mean(deltaE2000_errors(test_idx));
    max_deltaE = max(deltaE2000_errors(test_idx));
    
    % Compute the Delta E2000 error map for the full set (error map for entire dataset)
    error_map = reshape(deltaE2000(corrected_lab, ref_lab), m, n);
    
    % Display the results
    disp([space ' - Mean ΔE2000 Error: ', num2str(mean_deltaE)]);
    disp([space ' - Max ΔE2000 Error: ', num2str(max_deltaE)]);
    
    % Plot the error map
    figure;
    imagesc(error_map);
    colormap(jet);
    colorbar;
    clim([0 10]);  % Set the color axis limits for clarity
    title([space ' ΔE2000', ' (Mean: ', num2str(mean_deltaE, '%.2f'), ', Max: ', ...
        num2str(max_deltaE, '%.2f'), ')']);
    grid off;

    [rows, cols] = ind2sub([m, n], test_idx); % Convert linear indices to (row, col)
    hold on;  % Keep the error map displayed
    plot(cols, rows, 'w.', 'MarkerSize', 10);  % Plot white dots at training locations
    hold off;
    
    % Save the error map figure
    saveas(gcf, fullfile(output_folder, [space, file_name]));
end