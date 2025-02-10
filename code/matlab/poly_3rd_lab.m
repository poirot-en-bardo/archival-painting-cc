clc;
clear;
close all;

%% Define parameters
cube_folder = '../data/colorChecker_SG/cubes';  % Folder with HDR files
reference_cube_path = '../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr';  % Path to reference cube (same for all)
output_folder = '../results/error_maps/third_lab';  % Output folder for Lab-based error maps
rng(10);

% Create output folder if it doesn't exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Load spectral data
ill = importdata('../data/CIE_D65.txt'); % Illuminant
CMFs = importdata('../data/CIE2degCMFs_1931.txt');  % CIE 1931 CMFs

% Load the reference cube (same for all)
hcube_ref = hypercube(reference_cube_path);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;

[m, n, bd] = size(refCUBE);
lincube_ref = reshape(refCUBE, [], bd);

illIP = interp1(ill(:, 1), ill(:, 2), bands_ref, 'spline');
CMFsIP = [interp1(CMFs(:, 1), CMFs(:, 2), bands_ref, 'spline'), ...
          interp1(CMFs(:, 1), CMFs(:, 3), bands_ref, 'spline'), ...
          interp1(CMFs(:, 1), CMFs(:, 4), bands_ref, 'spline')];
sp_tristREF = CMFsIP .* illIP;

% Compute XYZ values for the reference cube
xyz_ref = (lincube_ref * sp_tristREF) ./ sum(sp_tristREF(:, 2), 1);

% Convert reference XYZ to Lab
lab_ref = xyz2lab(xyz_ref);

%% Loop through all HDR files in the folder (excluding the reference)
hdr_files = dir(fullfile(cube_folder, '*.hdr'));  % Get all .hdr files
for i = 1:length(hdr_files)
    if strcmp(hdr_files(i).name, 'reference_cube.hdr')  % Skip reference cube
        continue;
    end

    % Load the current cube to be corrected
    cubeFile = fullfile(hdr_files(i).folder, hdr_files(i).name);
    hcube = hypercube(cubeFile);
    inCUBE = hcube.DataCube;
    bands = hcube.Wavelength;

    [m, n, bd] = size(inCUBE);
    lincube = reshape(inCUBE, [], bd);

    % Interpolate illuminant and CMF to captured wavelengths
    illIP = interp1(ill(:, 1), ill(:, 2), bands, 'spline');
    CMFsIP = [interp1(CMFs(:, 1), CMFs(:, 2), bands, 'spline'), ...
              interp1(CMFs(:, 1), CMFs(:, 3), bands, 'spline'), ...
              interp1(CMFs(:, 1), CMFs(:, 4), bands, 'spline')];
    sp_tristREF = CMFsIP .* illIP;

    % Compute XYZ values of the current cube
    xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:, 2), 1);

    % Convert XYZ to Lab for correction
    lab_input = xyz2lab(xyz_input);

    %% **Third-Degree Polynomial Regression Model (Lab Space)**
    degree = 3;
    X = [lab_input, lab_input.^2, lab_input.^3, ...
         lab_input(:, 1).*lab_input(:, 2), lab_input(:, 1).*lab_input(:, 3), lab_input(:, 2).*lab_input(:, 3), ...
         lab_input(:, 1).^2 .* lab_input(:, 2), lab_input(:, 1).*lab_input(:, 2).^2, ...
         lab_input(:, 2).^2 .* lab_input(:, 3), lab_input(:, 1).*lab_input(:, 2).*lab_input(:, 3)];

    Y = lab_ref;  % Target Lab values

    % Train-test split (80% training, 20% testing)
    num_samples = size(X, 1);
    perm = randperm(num_samples);
    train_idx = perm(1:round(0.8 * num_samples));
    test_idx = perm(round(0.8 * num_samples) + 1:end);

    X_train = X(train_idx, :);
    Y_train = Y(train_idx, :);
    X_test = X(test_idx, :);
    Y_test = Y(test_idx, :);

    % Fit polynomial regression model (using least squares)
    % coeffs = (X_train' * X_train) \ (X_train' * Y_train);
    coeffs = (Y_train' * X_train * pinv(X_train' * X_train))';

    % Apply correction in Lab space
    corrected_lab = X * coeffs;

    % Clamp Lab values to valid ranges
    corrected_lab(:,1) = min(max(corrected_lab(:,1), 0), 100);  % L in [0,100]
    corrected_lab(:,2) = min(max(corrected_lab(:,2), -128), 127);  % a in [-128,127]
    corrected_lab(:,3) = min(max(corrected_lab(:,3), -128), 127);  % b in [-128,127]

    % Convert corrected Lab back to XYZ
    corrected_xyz = lab2xyz(corrected_lab);

    % Evaluate error for the current cube
    evaluate_error(X_test, Y_test, coeffs, Y, corrected_xyz, m, n, hdr_files(i).name, output_folder);
end

function evaluate_error(X_test, Y_test, coeffs, Y, corrected_xyz, m, n, file_name, output_folder)
    % Compute Delta E2000 Error
    lab_ref = xyz2lab(Y_test);
    lab_corrected = xyz2lab(X_test * coeffs);
    
    deltaE2000_errors = deltaE2000(lab_corrected, lab_ref);
    
    % Error map visualization
    lab_full_ref = xyz2lab(Y);
    lab_full_corrected = xyz2lab(corrected_xyz);
    error_map = reshape(deltaE2000(lab_full_corrected, lab_full_ref), m, n);
    
    % Calculate the min and max Delta E values
    min_deltaE = min(deltaE2000_errors);
    max_deltaE = max(deltaE2000_errors);

    % Extract the part of the filename after the first underscore
    [~, name, ~] = fileparts(file_name);
    underscore_pos = strfind(name, '_');
    if ~isempty(underscore_pos)
        name_part = name(underscore_pos(1)+1:end);  % Get substring after the first "_"
    else
        name_part = name;  % If no underscore, use the full filename
    end
    
    % Escape underscores in the name_part for the title
    name_part = strrep(name_part, '_', '\_');
    
    % Display the error map
    figure;
    imagesc(error_map);
    colormap(jet);
    colorbar;
    clim([0 10]); % Adjust the color axis limits
    grid off;
    
    % Add the file name part and Delta E values to the figure title (with 2 decimals)
    title(['ΔE2000 Error Map - ', name_part, ...
        ' (Min: ', num2str(min_deltaE, '%.2f'), ', Max: ', ...
        num2str(max_deltaE, '%.2f'), ')']);
    
    % Save the figure
    saveas(gcf, fullfile(output_folder, [name_part, '_error_map.png']));
    
    % Display summary statistics
    disp(['Mean ΔE2000 Error: ', num2str(mean(deltaE2000_errors))]);
    disp(['Max ΔE2000 Error: ', num2str(max(deltaE2000_errors))]);
end
