clc; clear; close all;
rng(10);

%% ---------------- Load Hyperspectral Data ----------------
cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_120f-ekta100-f14.hdr";
refFile  = "../../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

hcube = hypercube(cubeFile);
inCUBE = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd] = size(inCUBE);

hcube_ref = hypercube(refFile);
refCUBE = hcube_ref.DataCube;
bands_ref = hcube_ref.Wavelength;
[mr, nr, bd_ref] = size(refCUBE);

% Flatten the spectral data
lincube = reshape(inCUBE, [], bd);
lincube_ref = reshape(refCUBE, [], bd_ref);

%% ---------------- Compute XYZ from Spectral Data ----------------
% Load illuminant and CMFs
ill = importdata("../../../data/CIE_D65.txt");
CMFs_1931 = importdata("../../../data/CIE2degCMFs_1931.txt");
CMFs = CMFs_1931;

% For input cube:
illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');
CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...
          interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];
sp_tristREF = CMFsIP .* illIP;
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

% For reference cube:
illIP_ref = interp1(ill(:,1), ill(:,2), bands_ref, 'spline');
CMFsIP_ref = [interp1(CMFs(:,1), CMFs(:,2), bands_ref, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,3), bands_ref, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,4), bands_ref, 'spline')];
sp_tristREF_ref = CMFsIP_ref .* illIP_ref;
xyz_ref = (lincube_ref * sp_tristREF_ref) ./ sum(sp_tristREF_ref(:,2), 1);

%% ---------------- Reshape to 10x14 for Patch Indexing ----------------
% Here we assume the chart is exactly 10 rows by 14 columns.
xyz_input_2d = reshape(xyz_input, 10, 14, 3);
xyz_ref_2d   = reshape(xyz_ref,   10, 14, 3);

%% ---------------- Convert XYZ to RGB ----------------
% Convert using the same color space as your display.
rgb_input = xyz2rgb(xyz_input, 'ColorSpace', 'prophoto-rgb');
rgb_ref   = xyz2rgb(xyz_ref,   'ColorSpace', 'prophoto-rgb');

% Reshape into 10x14 for patch indexing
rgb_input_2d = reshape(rgb_input, 10, 14, 3);
rgb_ref_2d   = reshape(rgb_ref,   10, 14, 3);

%% ---------------- Train-Test Split in RGB ----------------
% Flatten the RGB data (each row is one patch with [R, G, B]).
rgb_input_flat = reshape(rgb_input_2d, [], 3);
rgb_ref_flat   = reshape(rgb_ref_2d, [], 3);

num_samples = size(rgb_input_flat, 1);
perm = randperm(num_samples);
trainCount = round(0.8 * num_samples);
train_idx = perm(1:trainCount);
initial_test_idx  = perm(trainCount+1:end);

rgb_input_train = rgb_input_flat(train_idx, :);
rgb_ref_train   = rgb_ref_flat(train_idx, :);

%% ----- Filter Test Set: Keep Only Points Inside Training Volume -----
% Create an alphaShape object of the training set. Using Inf ensures it equals the convex hull.
shp = alphaShape(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), Inf);

% Check which initial test points fall inside the training volume.
insideLogical = inShape(shp, rgb_input_flat(initial_test_idx,1), ...
                              rgb_input_flat(initial_test_idx,2), ...
                              rgb_input_flat(initial_test_idx,3));
filtered_test_idx = initial_test_idx(insideLogical);

% Optionally, if you want to retain approximately the same number as before, you can sample.
test_idx = filtered_test_idx;

%% ---------------- Build 3D LUT in RGB Space using Training Data ----------------
% Create scattered interpolants that map each training input [R,G,B] to the corresponding reference [R,G,B].
F_R = scatteredInterpolant(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), rgb_ref_train(:,1), 'linear', 'nearest');
F_G = scatteredInterpolant(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), rgb_ref_train(:,2), 'linear', 'nearest');
F_B = scatteredInterpolant(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), rgb_ref_train(:,3), 'linear', 'nearest');

%% ---------------- Apply LUT to Entire Data ----------------
rgb_corrected_flat = zeros(size(rgb_input_flat));
for i = 1:size(rgb_input_flat,1)
    rgb_corrected_flat(i,1) = F_R(rgb_input_flat(i,1), rgb_input_flat(i,2), rgb_input_flat(i,3));
    rgb_corrected_flat(i,2) = F_G(rgb_input_flat(i,1), rgb_input_flat(i,2), rgb_input_flat(i,3));
    rgb_corrected_flat(i,3) = F_B(rgb_input_flat(i,1), rgb_input_flat(i,2), rgb_input_flat(i,3));
end

%% ---------------- Evaluate Correction ----------------
output_folder = '../../../results/error_maps';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
[~, img_name, ~] = fileparts(cubeFile);
evaluate_error(rgb_ref_flat, rgb_corrected_flat, test_idx, 10, 14, 'RGB_3dLUT_', output_folder, [img_name, '_RGB3dLUT_correction.png']);

%% ---------------- Display Correction Results ----------------
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
imagesc(rgb_input_2d);
axis equal off;
title('Original Input (RGB)');

nexttile;
imagesc(rgb_ref_2d);
axis equal off;
title('Reference (RGB)');

% Reshape the corrected flat data back to a 10x14x3 image.
rgb_corrected_2d = reshape(rgb_corrected_flat, 10, 14, 3);
nexttile;
imagesc(rgb_corrected_2d);
axis equal off;
title('Corrected (3D LUT in RGB)');

sgtitle('Single-Step 3D LUT Correction in RGB Space');

%% ---------------- Visualize Training Volume and Test Points in RGB Space ----------------
figure('Name','Training vs Test Points with Volume in RGB Space');
% Plot training points in blue
scatter3(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), 80, 'b', 'filled');
hold on;
% Plot test points in red
scatter3(rgb_input_flat(test_idx,1), rgb_input_flat(test_idx,2), rgb_input_flat(test_idx,3), 80, 'r', 'filled');
% Compute the convex hull of the training set and plot it
K = convhull(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3));
trisurf(K, rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), ...
    'FaceAlpha', 0.3, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
xlabel('R'); ylabel('G'); zlabel('B');
legend('Training Set', 'Test Set', 'Training Volume','Location','best');
title('Training vs Test Points with Training Volume in RGB Space');
grid on;
view(3);
hold off;
