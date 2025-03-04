clc; clear; close all;
roof = double(intmax('uint16'));

%% Load two images
img_path1 = '/Volumes/School/Thesis/thesis-repo/data/scream/mod/modified_blend_0.5_0.8_exp.png';
img_path2 = '/Volumes/School/Thesis/thesis-repo/data/scream/scream_rgb.png';

% Read the images (assuming they are in uint16 format)
img1 = imread(img_path1);
img2 = imread(img_path2);

% Convert images to grayscale (if needed) using the standard luminance formula
img1_grey = rgb2gray(img1);
img2_grey = rgb2gray(img2);

% Normalize the grayscale images to the [0, 1] range
img1_grey = double(img1_grey) / double(roof);
img2_grey = double(img2_grey) / double(roof);

%% Fourier Transform for both images
F1 = fft2(img1_grey);
F2 = fft2(img2_grey);

% Shift zero frequency to center
F1_shifted = fftshift(F1);
F2_shifted = fftshift(F2);

%% Extract low and high frequencies
[rows, cols] = size(F1_shifted);

% Radius for low-frequency cutoff
low_freq_radius = min(rows, cols) / 500;

% Create grid of distance values from the center
[X, Y] = meshgrid(1:cols, 1:rows);
centerX = floor(cols / 2);
centerY = floor(rows / 2);
distance = sqrt((X - centerX).^2 + (Y - centerY).^2);

% Create low-frequency mask (circle)
low_freq_mask = distance <= low_freq_radius;

% Create high-frequency mask (everything else)
high_freq_mask = distance > low_freq_radius;

% Apply the masks to extract the low and high frequencies
low_freq1 = F1_shifted .* low_freq_mask;
high_freq1 = F1_shifted .* high_freq_mask;

low_freq2 = F2_shifted .* low_freq_mask;
high_freq2 = F2_shifted .* high_freq_mask;

% Inverse FFT to recover the images from the frequency components
low_freq1_spatial = ifft2(ifftshift(low_freq1));  % Inverse FFT
high_freq1_spatial = ifft2(ifftshift(high_freq1)); % Inverse FFT

low_freq2_spatial = ifft2(ifftshift(low_freq2));  % Inverse FFT
high_freq2_spatial = ifft2(ifftshift(high_freq2)); % Inverse FFT

% Take the magnitude of the recovered spatial domain images
low_freq1_spatial = abs(low_freq1_spatial);
high_freq1_spatial = abs(high_freq1_spatial);

low_freq2_spatial = abs(low_freq2_spatial);
high_freq2_spatial = abs(high_freq2_spatial);

%% Compute the differences in the spatial domain
low_freq_diff_spatial = abs(low_freq1_spatial - low_freq2_spatial);
high_freq_diff_spatial = abs(high_freq1_spatial - high_freq2_spatial);

%% Normalize the differences to [0, 1] for display
low_freq_diff_spatial_norm = mat2gray(low_freq_diff_spatial);
high_freq_diff_spatial_norm = mat2gray(high_freq_diff_spatial);

%% Display the differences as heatmaps

% Create figure for displaying low frequency difference
figure;
imshow(low_freq_diff_spatial_norm, []);  % Display low-frequency difference
colorbar;
title('Difference in Low Frequencies');
axis equal;

% Save the low frequency difference heatmap with title and colorbar
saveas(gcf, 'stripe_low_freq_diff.png');  % Save figure as PNG

% Create figure for displaying high frequency difference
figure;
imshow(high_freq_diff_spatial_norm, []);  % Display high-frequency difference
colorbar;
title('Difference in High Frequencies');
axis equal;

% Save the high frequency difference heatmap with title and colorbar
saveas(gcf, 'stripe_high_freq_diff.png');  % Save figure as PNG
