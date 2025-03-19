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

% Flatten spectral data
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
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2),1);

% For reference cube:
illIP_ref = interp1(ill(:,1), ill(:,2), bands_ref, 'spline');
CMFsIP_ref = [interp1(CMFs(:,1), CMFs(:,2), bands_ref, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,3), bands_ref, 'spline'), ...
              interp1(CMFs(:,1), CMFs(:,4), bands_ref, 'spline')];
sp_tristREF_ref = CMFsIP_ref .* illIP_ref;
xyz_ref = (lincube_ref * sp_tristREF_ref) ./ sum(sp_tristREF_ref(:,2),1);

%% ---------------- Reshape to 10x14 ----------------
% The chart is 10x14 patches.
xyz_input_2d = reshape(xyz_input, 10, 14, 3);
xyz_ref_2d   = reshape(xyz_ref, 10, 14, 3);

%% ---------------- Convert XYZ to RGB ----------------
% Use built-in conversion (here using 'prophoto-rgb')
rgb_input = xyz2rgb(xyz_input, 'ColorSpace', 'prophoto-rgb');
rgb_ref   = xyz2rgb(xyz_ref,   'ColorSpace', 'prophoto-rgb');

% Reshape to 10x14 for display if needed
rgb_input_2d = reshape(rgb_input, 10, 14, 3);
rgb_ref_2d   = reshape(rgb_ref, 10, 14, 3);

%% ---------------- Train-Test Split ----------------
% Flatten the Lab data (each row is one patch)
rgb_input_flat = reshape(rgb_input_2d, [], 3);
rgb_ref_flat   = reshape(rgb_ref_2d, [], 3);

% Use 80% of the data for training
num_samples = size(rgb_input_flat, 1);
perm = randperm(num_samples);
trainCount = round(0.8 * num_samples);
train_idx = perm(1:trainCount);
test_idx  = perm(trainCount+1:end);

rgb_input_train = rgb_input_flat(train_idx, :);
rgb_ref_train   = rgb_ref_flat(train_idx, :);

%% ---------------- Build 3D LUT in Lab Space using Training Data ----------------
% Create a scattered interpolant for each Lab channel mapping input -> reference.
F_R = scatteredInterpolant(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), rgb_ref_train(:,1), 'linear', 'nearest');
F_G = scatteredInterpolant(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), rgb_ref_train(:,2), 'linear', 'nearest');
F_B = scatteredInterpolant(rgb_input_train(:,1), rgb_input_train(:,2), rgb_input_train(:,3), rgb_ref_train(:,3), 'linear', 'nearest');

%% ---------------- Apply LUT to Entire Data ----------------
rgb_corrected = zeros(size(rgb_input_flat));
for i = 1:size(rgb_input_flat,1)
    rgb_corrected(i,1) = F_R(rgb_input_flat(i,1), rgb_input_flat(i,2), rgb_input_flat(i,3));
    rgb_corrected(i,2) = F_G(rgb_input_flat(i,1), rgb_input_flat(i,2), rgb_input_flat(i,3));
    rgb_corrected(i,3) = F_B(rgb_input_flat(i,1), rgb_input_flat(i,2), rgb_input_flat(i,3));
end

%% ---------------- Evaluate Correction ----------------
% For evaluation, convert both corrected and reference XYZ to Lab.
lab_corrected = rgb2lab(rgb_corrected, 'ColorSpace', 'prophoto-rgb');
lab_ref = xyz2lab(xyz_ref);

% Evaluate error (assumes evaluate_error exists)
% Here we use full dataset evaluation.
output_folder = '../../../results/error_maps';
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
[~, img_name, ~] = fileparts(cubeFile);
evaluate_error(lab_ref, lab_corrected, test_idx, 10, 14, '3dLUT_', output_folder, [img_name, '_3dLUT_correction.png']);

%% ---------------- Display Results ----------------
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

nexttile;
imagesc(reshape(rgb_corrected, 10, 14, 3));
axis equal off;
title('Corrected (3D LUT)');

sgtitle('3D LUT Correction from Input to Reference');



%% Visualize 3D LUT as a Color Cloud in RGB Space
% Determine the range of input RGB values (assuming they lie between 0 and 1)
min_val = min(rgb_input_train);
max_val = max(rgb_input_train);
numSamples = 17;  % Adjust resolution of grid sampling

% Create a 3D grid of sample points in the input RGB space
[r_grid, g_grid, b_grid] = ndgrid( linspace(min_val(1), max_val(1), numSamples), ...
                                   linspace(min_val(2), max_val(2), numSamples), ...
                                   linspace(min_val(3), max_val(3), numSamples) );
samplePoints = [r_grid(:), g_grid(:), b_grid(:)];

% Evaluate the LUT mapping at these sample points using the scattered interpolants
mapped_R = F_R(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3));
mapped_G = F_G(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3));
mapped_B = F_B(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3));
mappedPoints = [mapped_R, mapped_G, mapped_B];

% Visualize the input sample points and the LUT-mapped points side-by-side
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);

subplot(1,2,1);
scatter3(samplePoints(:,1), samplePoints(:,2), samplePoints(:,3), 36, samplePoints, 'filled');
xlabel('R'); ylabel('G'); zlabel('B');
title('Input Color Space Samples');
axis([0 1 0 1 0 1]); grid on; colorbar;

subplot(1,2,2);
scatter3(mappedPoints(:,1), mappedPoints(:,2), mappedPoints(:,3), 36, mappedPoints, 'filled');
xlabel('R'); ylabel('G'); zlabel('B');
title('LUT-Mapped Color Space');
axis([0 1 0 1 0 1]); grid on; colorbar;

sgtitle('3D LUT Visualization: Input vs. Mapped Colors');

%%

% Flatten the 10x14 images to N-by-3 matrices for visualization:
rgb_input_flat = reshape(rgb_input_2d, [], 3);
rgb_ref_flat   = reshape(rgb_ref_2d, [], 3);
% If rgb_corrected is not already flattened, do so:
rgb_corrected_flat = reshape(rgb_corrected, [], 3);

% Visualize the three point clouds in one figure:
figure('Units','normalized','OuterPosition',[0 0 1 1]);

subplot(1,2,1);
scatter3(rgb_input_flat(:,1), rgb_input_flat(:,2), rgb_input_flat(:,3), 36, rgb_input_flat, 'filled');
xlabel('R'); ylabel('G'); zlabel('B');
title('Input RGB Cloud');
axis([0 1 0 1 0 1]); grid on; colorbar;

% subplot(1,3,2);
% scatter3(rgb_ref_flat(:,1), rgb_ref_flat(:,2), rgb_ref_flat(:,3), 36, rgb_ref_flat, 'filled');
% xlabel('R'); ylabel('G'); zlabel('B');
% title('Reference RGB Cloud');
% axis([0 1 0 1 0 1]); grid on; colorbar;

subplot(1,2,2);
scatter3(rgb_corrected_flat(:,1), rgb_corrected_flat(:,2), rgb_corrected_flat(:,3), 36, rgb_corrected_flat, 'filled');
xlabel('R'); ylabel('G'); zlabel('B');
title('Corrected RGB Cloud');
axis([0 1 0 1 0 1]); grid on; colorbar;

sgtitle('RGB Color Clouds: Input, Reference, Corrected');
