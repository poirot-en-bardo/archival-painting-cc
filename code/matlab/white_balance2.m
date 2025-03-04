%% Define parameters
clc; clear; close all;
roof = double(intmax('uint16')); % Histogram normalization factor

% Load spectral data
ill = importdata('../../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../../data/CIE2degCMFs_1931.txt');
CMFs = CMFs_1931;
rng(10);

%% Select the reference and the cube to be corrected

cubeFile = "../../data/colorChecker_SG/cubes/cubeCC_ekta100-studio13.hdr";
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

rgb_input =  xyz2rgb(xyz_input, 'ColorSpace','prophoto-rgb');

%% Assume xyz_input is an N×3 matrix (with N = m*n) obtained from your data
% and that m and n are the spatial dimensions of the color checker.

%% Assume xyz_input (N×? matrix) is computed from your data, and m and n are the spatial dimensions.
% Convert XYZ to RGB (using ProPhoto-RGB in this example)
rgb_image = reshape(xyz2rgb(xyz_input, 'ColorSpace','prophoto-rgb'), [m, n, 3]);

%% Define the white patch positions (based on your known pattern)
% White patches are at specific rows and columns:
white_patch_rows = [1, 4, 7, 10];
white_patch_cols = [1, 4, 7, 10, 14];

[gridRows, gridCols] = ndgrid(1:m, 1:n);
white_mask = ismember(gridRows, white_patch_rows) & ismember(gridCols, white_patch_cols);

% Get coordinates of the white patches:
[row_white, col_white] = find(white_mask);

% Extract the RGB values at those white patches:
numWhite = numel(row_white);
white_rgb_values = zeros(numWhite, 3);
for k = 1:numWhite
    white_rgb_values(k, :) = squeeze(rgb_image(row_white(k), col_white(k), :))';
end

%% Compute the mean white value for each channel
mean_white = mean(white_rgb_values, 1);  % [mean_R, mean_G, mean_B]

%% Compute relative gain factors at the white patch positions
gain_r_values = mean_white(1) ./ white_rgb_values(:,1);
gain_g_values = mean_white(2) ./ white_rgb_values(:,2);
gain_b_values = mean_white(3) ./ white_rgb_values(:,3);

%% Interpolate gain maps over the entire image
% Use scatteredInterpolant to interpolate the gains from the white patches.
%% Use scatteredInterpolant to interpolate the gains from the white patches.
F_r = scatteredInterpolant(col_white, row_white, gain_r_values, 'linear', 'nearest');
F_g = scatteredInterpolant(col_white, row_white, gain_g_values, 'linear', 'nearest');
F_b = scatteredInterpolant(col_white, row_white, gain_b_values, 'linear', 'nearest');

% Interpolation grid: we need to transpose col_grid and row_grid so they match the gain map's orientation
[col_grid, row_grid] = meshgrid(1:n, 1:m);  % col_grid = 1:n, row_grid = 1:m
gain_r_map = F_r(col_grid, row_grid);  % Interpolated gain map for red channel
gain_g_map = F_g(col_grid, row_grid);  % Interpolated gain map for green channel
gain_b_map = F_b(col_grid, row_grid);  % Interpolated gain map for blue channel

%% Apply spatial compensation mask (smooth transition for row/column balancing)
% Apply a smooth mask that compensates for rows and columns by using interpolation

% Compensation factors for rows and columns (scale factors that gradually change)
row_compensation = linspace(1, 1.05, m); % Compensation factors for rows (1 to 1.05 over m rows)
col_compensation = linspace(1, 1.05, n); % Compensation factors for columns (1 to 1.05 over n columns)

% Create a meshgrid that expands these factors across the entire image
[row_comp_mask, col_comp_mask] = meshgrid(row_compensation, col_compensation);

% Ensure the mask sizes match the image size (m x n)
% No need to transpose because now the gain maps and the compensation masks are aligned.

%% Apply compensation to gain maps
gain_r_map = gain_r_map .* row_comp_mask .* col_comp_mask;
gain_g_map = gain_g_map .* row_comp_mask .* col_comp_mask;
gain_b_map = gain_b_map .* row_comp_mask .* col_comp_mask;

%% Clamp the gain maps to avoid extreme corrections (tweak the limits as needed)
gain_min = 0.8; gain_max = 1.2;
gain_r_map = max(min(gain_r_map, gain_max), gain_min);
gain_g_map = max(min(gain_g_map, gain_max), gain_min);
gain_b_map = max(min(gain_b_map, gain_max), gain_min);




%% Apply the gain maps to the image (per channel)
balanced_rgb = zeros(size(rgb_image));
balanced_rgb(:,:,1) = rgb_image(:,:,1) .* gain_r_map;
balanced_rgb(:,:,2) = rgb_image(:,:,2) .* gain_g_map;
balanced_rgb(:,:,3) = rgb_image(:,:,3) .* gain_b_map;

% Clip values to [0, 1]
balanced_rgb = min(max(balanced_rgb, 0), 1);

%% Enlarge the image for visualization (e.g., using nearest neighbor interpolation)
scaleFactor = 100;  % Adjust scaleFactor as needed
enlarged_image = imresize(balanced_rgb, scaleFactor, 'nearest');

figure;
imshow(enlarged_image);
title('Spatially Variant White Balanced Image');
