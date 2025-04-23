clc; clear; close all;
rng(10);

%% ---------------- Load Hyperspectral Data ----------------
cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_120f-velvia-f8.hdr";
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

%% ---------------- Convert XYZ to Lab ----------------
lab_input = xyz2lab(xyz_input);
lab_ref = xyz2lab(xyz_ref);

% Reshape into 10x14 for patch indexing
lab_input_2d = reshape(lab_input, 10, 14, 3);
lab_ref_2d   = reshape(lab_ref,   10, 14, 3);

%% ---------------- Train-Test Split ----------------
% Flatten the Lab data (each row is one patch)
lab_input_flat = reshape(lab_input_2d, [], 3);
lab_ref_flat   = reshape(lab_ref_2d, [], 3);

% Use 80% of the data for training
num_samples = size(lab_input_flat, 1);
perm = randperm(num_samples);
trainCount = round(0.8 * num_samples);
train_idx = perm(1:trainCount);
test_idx  = perm(trainCount+1:end);

lab_input_train = lab_input_flat(train_idx, :);
lab_ref_train   = lab_ref_flat(train_idx, :);

%% ---------------- Build 3D LUT in Lab Space using Training Data ----------------
% Create a scattered interpolant for each Lab channel mapping input -> reference.
F_L = scatteredInterpolant(lab_input_train(:,1), lab_input_train(:,2), lab_input_train(:,3), lab_ref_train(:,1), 'linear', 'nearest');
F_a = scatteredInterpolant(lab_input_train(:,1), lab_input_train(:,2), lab_input_train(:,3), lab_ref_train(:,2), 'linear', 'nearest');
F_b = scatteredInterpolant(lab_input_train(:,1), lab_input_train(:,2), lab_input_train(:,3), lab_ref_train(:,3), 'linear', 'nearest');

%% ---------------- Apply LUT to Entire Data ----------------
lab_corrected_flat = zeros(size(lab_input_flat));
for i = 1:size(lab_input_flat,1)
    lab_corrected_flat(i,1) = F_L(lab_input_flat(i,1), lab_input_flat(i,2), lab_input_flat(i,3));
    lab_corrected_flat(i,2) = F_a(lab_input_flat(i,1), lab_input_flat(i,2), lab_input_flat(i,3));
    lab_corrected_flat(i,3) = F_b(lab_input_flat(i,1), lab_input_flat(i,2), lab_input_flat(i,3));
end

%% ---------------- Evaluate Correction ----------------
% Evaluate error by comparing the corrected Lab to the reference Lab
% (Assuming evaluate_error is a function that computes ΔE and displays an error map)
output_folder = '../../../results/error_maps';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
[~, img_name, ~] = fileparts(cubeFile);
evaluate_error(lab_ref_flat, lab_corrected_flat, test_idx, 10, 14, 'Lab_3dLUT_', output_folder, [img_name, '_Lab3dLUT_correction.png']);

%% ---------------- Display Results ----------------
% For visualization, convert the corrected Lab back to RGB.
% If lab2rgb is available, use it; otherwise, convert Lab->XYZ->RGB.

rgb_corrected = lab2rgb(reshape(lab_corrected_flat, 10, 14, 3), 'ColorSpace', 'prophoto-rgb');

rgb_input_disp = xyz2rgb(xyz_input_2d, 'ColorSpace', 'prophoto-rgb');
rgb_ref_disp   = xyz2rgb(xyz_ref_2d, 'ColorSpace', 'prophoto-rgb');

figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
imagesc(rgb_input_disp);
axis equal off;
title('Original Input (RGB)');

nexttile;
imagesc(rgb_ref_disp);
axis equal off;
title('Reference (RGB)');

nexttile;
imagesc(rgb_corrected);
axis equal off;
title('Corrected (3D LUT in Lab)');

sgtitle('Single-Step 3D LUT Correction in Lab Space');



%% ---------------- Visualize Training vs Test Points ----------------
figure('Name','Training vs Test Points in Lab Space');
% Plot training points in blue
scatter3(lab_input_train(:,1), lab_input_train(:,2), lab_input_train(:,3), 80, 'b', 'filled');
hold on;
% Plot test points in red
scatter3(lab_input_flat(test_idx,1), lab_input_flat(test_idx,2), lab_input_flat(test_idx,3), 80, 'r', 'filled');
xlabel('L*');
ylabel('a*');
zlabel('b*');
legend('Training Set', 'Test Set','Location','best');
title('Visualization of Training vs Test Points in Lab Space');
grid on;
view(3);

%% ---------------- Visualize Training vs Test Points with Training Volume ----------------
figure('Name','Training vs Test Points in Lab Space with Volume');
% Plot training points in blue
scatter3(lab_input_train(:,1), lab_input_train(:,2), lab_input_train(:,3), 80, 'b', 'filled');
hold on;
% Plot test points in red
scatter3(lab_input_flat(test_idx,1), lab_input_flat(test_idx,2), lab_input_flat(test_idx,3), 80, 'r', 'filled');

% Compute the convex hull of the training set in Lab space
K = convhull(lab_input_train(:,1), lab_input_train(:,2), lab_input_train(:,3));
% Plot the convex hull as a semi-transparent patch
trisurf(K, lab_input_train(:,1), lab_input_train(:,2), lab_input_train(:,3), ...
    'FaceAlpha', 0.3, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');

xlabel('L*');
ylabel('a*');
zlabel('b*');
legend('Training Set', 'Test Set','Training Volume','Location','best');
title('Training vs Test Points - Lab Space');
grid on;
view(3);
hold off;
