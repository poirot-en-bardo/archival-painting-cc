clc; clear; close all;
rng(10);

% Define the folder containing the hyperspectral data cubes
dataFolder = "../../../data/colorChecker_SG/cubes/";

% Get a list of all .hdr files in the folder
cubeFiles = dir(fullfile(dataFolder, '*.hdr'));

% Reference file (assumed to be the same for all cubes)
refFile = "../../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

% Load reference hypercube
hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
[mr, nr, bd_ref] = size(refCUBE);

% Flatten the reference spectral data
lincube_ref = reshape(refCUBE, [], bd_ref);

% Load illuminant and color matching functions (CMFs)
ill = importdata("../../../data/CIE_D65.txt");
CMFs_1931 = importdata("../../../data/CIE2degCMFs_1931.txt");
CMFs = CMFs_1931;

% Interpolate illuminant and CMFs for reference bands
illIP_ref = interp1(ill(:,1), ill(:,2), bands_ref, 'spline');
CMFsIP_ref = [interp1(CMFs(:,1), CMFs(:,2), bands_ref, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,3), bands_ref, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,4), bands_ref, 'spline')];
sp_tristREF_ref = CMFsIP_ref .* illIP_ref;

% Compute XYZ for reference cube
xyz_ref = (lincube_ref * sp_tristREF_ref) ./ sum(sp_tristREF_ref(:,2), 1);

% Reshape to 10x14 for patch indexing
xyz_ref_2d = reshape(xyz_ref, 10, 14, 3);

% Convert reference XYZ to RGB
rgb_ref = xyz2rgb(xyz_ref, 'ColorSpace', 'prophoto-rgb');
rgb_ref_2d = reshape(rgb_ref, 10, 14, 3);

% Loop through each hyperspectral data cube in the folder
for k = 1:length(cubeFiles)
    % Construct the full path to the current cube file
    cubeFile = fullfile(cubeFiles(k).folder, cubeFiles(k).name);

    % Load input hypercube
    hcube = hypercube(cubeFile);
    inCUBE = hcube.DataCube;
    bands = hcube.Wavelength;
    [m, n, bd] = size(inCUBE);

    % Flatten the input spectral data
    lincube = reshape(inCUBE, [], bd);

    % Interpolate illuminant and CMFs for input bands
    illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
    CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
    sp_tristREF = CMFsIP .* illIP;

    % Compute XYZ for input cube
    xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

    % Reshape to 10x14 for patch indexing
    xyz_input_2d = reshape(xyz_input, 10, 14, 3);

    % Convert input XYZ to RGB
    rgb_input = xyz2rgb(xyz_input, 'ColorSpace', 'prophoto-rgb');
    rgb_input_2d = reshape(rgb_input, 10, 14, 3);

    % Flatten the RGB data
    rgb_input_flat = reshape(rgb_input_2d, [], 3);
    rgb_ref_flat = reshape(rgb_ref_2d, [], 3);

    % Train-test split in RGB
    num_samples = size(rgb_input_flat, 1);
    perm = randperm(num_samples);
    trainCount = round(0.8 * num_samples);
    train_idx = perm(1:trainCount);
    initial_test_idx = perm(trainCount+1:end);

    rgb_input_train = rgb_input_flat(train_idx, :);
    rgb_ref_train = rgb_ref_flat(train_idx, :);

    % Filter test set: Keep only points inside training volume
    shp = alphaShape(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), Inf);
    insideLogical = inShape(shp, rgb_input_flat(initial_test_idx,1), ...
                                  rgb_input_flat(initial_test_idx,2), ...
                                  rgb_input_flat(initial_test_idx,3));
    filtered_test_idx = initial_test_idx(insideLogical);
    test_idx = filtered_test_idx;

    % Build 3D LUT in RGB space using training data
    F_R = scatteredInterpolant(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), rgb_ref_train(:,1), 'linear', 'nearest');
    F_G = scatteredInterpolant(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), rgb_ref_train(:,2), 'linear', 'nearest');
    F_B = scatteredInterpolant(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), rgb_ref_train(:,3), 'linear', 'nearest');

    % Apply LUT to entire data
    rgb_corrected_flat = zeros(size(rgb_input_flat));
    for i = 1:size(rgb_input_flat,1)
        rgb_corrected_flat(i,1) = F_R(rgb_input_flat(i,1), rgb_input_flat(i,2), rgb_input_flat(i,3));
        rgb_corrected_flat(i,2) = F_G(rgb_input_flat(i,1), rgb_input_flat(i,2), rgb_input_flat(i,3));
        rgb_corrected_flat(i,3) = F_B(rgb_input_flat(i,1), rgb_input_flat(i,2), rgb_input_flat(i,3));
    end

    % Evaluate correction
    output_folder = '../../../results/error_maps';
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end
    [~, img_name, ~] = fileparts(cubeFile);
    evaluate_error(rgb_ref_flat, rgb_corrected_flat, test_idx, 10, 14, ['RGB_3dLUT_', img_name], output_folder, [img_name, '_RGB3dLUT_correction.png']);

    % % Display correction results
    % figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    % tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');

    % nexttile;
    % imagesc(rgb_input_2d);
    % axis equal off;
    % title('Original Input (RGB)');
    % 
    % nexttile;
    % imagesc(rgb_ref_2d);
    % axis equal off;
    % title('Reference (RGB)');
    % 
    % % Reshape the corrected flat data back to a 10x14x3 image.
    % rgb_corrected_2d = reshape(rgb_corrected_flat, 10, 14, 3);
    % nexttile;
    % imagesc(rgb_corrected_2d);
    % axis equal off;
    % title('Corrected (3D LUT in RGB)');

    %% ---------------- Visualize Training Volume and Test Points in RGB Space ----------------
    % figure('Name', ['Training vs Test Points: ', img_name]);
    % % Plot training points in blue
    % scatter3(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), 80, 'b', 'filled');
    % hold on;
    % % Plot test points in red
    % scatter3(rgb_input_flat(test_idx,1), rgb_input_flat(test_idx,2), rgb_input_flat(test_idx,3), 80, 'r', 'filled');
    % % Compute the convex hull of the training set and plot it
    % K = convhull(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3));
    % trisurf(K, rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), ...
    %     'FaceAlpha', 0.3, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
    % xlabel('R'); ylabel('G'); zlabel('B');
    % legend('Training Set', 'Test Set', 'Training Volume', 'Location', 'best');
    % title(['Training vs Test Points in RGB Space: ', img_name]);
    % grid on;
    % view(3);
    % hold off;

end % End of loop through all hyperspectral cubes

   
 
