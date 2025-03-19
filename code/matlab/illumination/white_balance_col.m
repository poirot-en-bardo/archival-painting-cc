%% Define parameters
clc; clear; close all;
roof = double(intmax('uint16')); % Histogram normalization factor

% Load spectral data
ill = importdata('../../../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../../../data/CIE2degCMFs_1931.txt');
CMFs = CMFs_1931;
rng(10);

%% Select the reference and the cube to be corrected
cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_120f-velvia-f8.hdr";
refFile = "../../../data/colorChecker_SG/cubeCC_DigitalSG_REF.hdr";

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

% Convert XYZ to RGB (using ProPhoto-RGB in this example)
rgb_input =  xyz2rgb(xyz_input, 'ColorSpace','prophoto-rgb');
rgb_input = reshape(rgb_input, [m, n, 3]); % Ensure this is a 3D RGB matrix

%% Define the white, grey, and black patch positions (based on your known pattern)
% Define the color checker pattern for white, grey, black patches
patch_rows = [1, 4, 7, 10, 1, 4, 7, 10, 1, 4, 7, 10, 1, 4]; % Example rows for each patch
patch_cols = [1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12, 13, 14]; % Example columns for each patch

% Define the expected RGB values for white, grey, and black patches
% White = [1, 1, 1], Grey = [0.5, 0.5, 0.5], Black = [0, 0, 0]
patch_rgb_values = zeros(length(patch_rows), 3);

for i = 1:length(patch_rows)
    patch_rgb_values(i, :) = squeeze(rgb_input(patch_rows(i), patch_cols(i), :))';
end

% Set up the expected RGB values for white, grey, and black patches
expected_rgb_white = [1, 1, 1];  % white patch RGB
expected_rgb_grey = [0.5, 0.5, 0.5];  % grey patch RGB
expected_rgb_black = [0, 0, 0];  % black patch RGB

% Now compute the gain factors for each patch type (white, grey, black)
gain_factors = zeros(size(patch_rgb_values));  % To store gain factors for each patch

for i = 1:length(patch_rows)
    if all(patch_rgb_values(i, :) == expected_rgb_white)
        gain_factors(i, :) = expected_rgb_white ./ patch_rgb_values(i, :);  % White patches
    elseif all(patch_rgb_values(i, :) == expected_rgb_grey)
        gain_factors(i, :) = expected_rgb_grey ./ patch_rgb_values(i, :);  % Grey patches
    elseif all(patch_rgb_values(i, :) == expected_rgb_black)
        gain_factors(i, :) = [1, 1, 1];  % Black patches (no correction)
    end
end

% Visualize gain factors for each channel
figure;
subplot(1, 3, 1);
imagesc(gain_factors(:, 1));
colorbar;
title('Red Gain Factors');

subplot(1, 3, 2);
imagesc(gain_factors(:, 2));
colorbar;
title('Green Gain Factors');

subplot(1, 3, 3);
imagesc(gain_factors(:, 3));
colorbar;
title('Blue Gain Factors');

%% Interpolate gain maps (without multiplying across rows and columns)
% Use scatteredInterpolant to interpolate gain values for each color channel
F_r = scatteredInterpolant(patch_cols(:), patch_rows(:), gain_factors(:, 1), 'linear', 'nearest');
F_g = scatteredInterpolant(patch_cols(:), patch_rows(:), gain_factors(:, 2), 'linear', 'nearest');
F_b = scatteredInterpolant(patch_cols(:), patch_rows(:), gain_factors(:, 3), 'linear', 'nearest');

[col_grid, row_grid] = meshgrid(1:n, 1:m);

gain_r_map = F_r(col_grid, row_grid);
gain_g_map = F_g(col_grid, row_grid);
gain_b_map = F_b(col_grid, row_grid);

% Normalize gain maps to ensure no extreme values
gain_min = 0.1; 
gain_max = 2;  % Increase the maximum gain to make the image brighter
gain_r_map = max(min(gain_r_map, gain_max), gain_min);
gain_g_map = max(min(gain_g_map, gain_max), gain_min);
gain_b_map = max(min(gain_b_map, gain_max), gain_min);

%% Apply the gain maps independently to each RGB channel
balanced_rgb = zeros(size(rgb_input));
balanced_rgb(:,:,1) = rgb_input(:,:,1) .* gain_r_map;
balanced_rgb(:,:,2) = rgb_input(:,:,2) .* gain_g_map;
balanced_rgb(:,:,3) = rgb_input(:,:,3) .* gain_b_map;

% Clip values to [0, 1] to avoid overflow/underflow
balanced_rgb = min(max(balanced_rgb, 0), 1);

%% Visualization
scaleFactor = 100;
enlarged_image = imresize(balanced_rgb, scaleFactor, 'nearest');

figure;
imshow(enlarged_image);
title('Balanced Image with White, Grey, and Black Patches');
