clc; clear; close all;
roof = double(intmax('uint16'));

%% Loading data
ill = importdata('../../../data/CIE_D65.txt'); % Illuminant
CMFs = importdata('../../../data/CIE2degCMFs_1931.txt'); % Color matching functions

path_cube = '/Volumes/School/Thesis/scream_hsi/inpainted_scream.hdr';

hcube = hypercube(path_cube);
cube = hcube.DataCube;
bands = hcube.Wavelength;
[m, n, bd] = size(cube);
cube_lin = reshape(cube, [], bd);

visible_idx = (bands >= 380) & (bands <= 780);
bands_vis = bands(visible_idx);

% Subset the hyperspectral cube to the visible wavelengths
cube_vis = cube(:,:,visible_idx);
cube_vis_lin = reshape(cube_vis, [], sum(visible_idx));

%% Computing the RGB image
% Interpolate illuminant and CMFs to captured wavelengths
ill_interp = interp1(ill(:,1), ill(:,2), bands_vis, 'spline');
CMFs_interp = interp1(CMFs(:,1), CMFs(:,2:4), bands_vis, 'spline');
sp_tristREF = CMFs_interp .* ill_interp;

% Calculate XYZ tristimulus values
xyz = (cube_vis_lin * sp_tristREF) ./ sum(sp_tristREF(:,2),1);
rgb = xyz2rgb(xyz, 'ColorSpace', 'adobe-rgb-1998');
scream_rgb = uint16(reshape(rgb, m, n, 3) .* roof);

%% Saving the RGB image
% imwrite(scream_rgb, 'scream_rgb.png', 'BitDepth', 16);

%% Fourier Analysis

img_path = '/Volumes/School/Thesis/thesis-repo/code/matlab/img_diff/deltae.png';
scream_rgb = imread(img_path);

%%
% Convert RGB to greyscale (using standard luminance formula)
scream_grey = rgb2gray(scream_rgb);

% Normalise to double for FFT ([0, 1])
scream_grey = double(scream_grey) / double(roof);
%%
img_path = '/Volumes/School/Thesis/thesis-repo/code/matlab/img_diff/deltae.png';
scream_grey = imread(img_path);
F = fft2(scream_grey);
F_shifted = fftshift(F); % Shift zero frequency to center


magnitude_spectrum = mat2gray(log(abs(F_shifted)));

% % Displaying the spectrum
% figure;
% imshow(magnitude_spectrum);
% title('Fourier Spectrum of The Scream');


% imwrite(uint8(255 * magnitude_spectrum), 'fourier_spectrum.png');



% Extracting low and high frequencies

[rows, cols] = size(F_shifted);

% radius for low-frequency cutoff (in pixels)
low_freq_radius = min(rows, cols) / 100;  

% grid of distance values from the center
[X, Y] = meshgrid(1:cols, 1:rows);
centerX = floor(cols / 2);
centerY = floor(rows / 2);
distance = sqrt((X - centerX).^2 + (Y - centerY).^2);

% Create low-frequency mask (circle)
low_freq_mask = distance <= low_freq_radius;

% Create high-frequency mask (everything else)
high_freq_mask = distance > low_freq_radius;

% Apply masks to get low and high-frequency components directly from F_shifted
low_freq = F_shifted .* low_freq_mask;
high_freq = F_shifted .* high_freq_mask;

% Get the magnitude spectra of low and high frequencies
low_freq_magnitude = log(1 + abs(low_freq));
high_freq_magnitude = log(1 + abs(high_freq));

% Normalize to [0, 1] for display
low_freq_magnitude_norm = mat2gray(low_freq_magnitude);
high_freq_magnitude_norm = mat2gray(high_freq_magnitude);

% Display and save the low-frequency and high-frequency images
% figure;
% imshow(low_freq_magnitude_norm);
% title('Low Frequencies');
% 
% figure;
% imshow(high_freq_magnitude_norm);
% title('High Frequencies');

%%
% Apply inverse FFT to get the spatial domain images 
low_freq_recovered = ifft2(ifftshift(low_freq));  % Inverse FFT 
low_freq_recovered = abs(low_freq_recovered);    % Take the magnitude 
low_freq_recovered = mat2gray(low_freq_recovered);

high_freq_recovered = ifft2(ifftshift(high_freq)); 
high_freq_recovered = abs(high_freq_recovered);    
high_freq_recovered = mat2gray(high_freq_recovered); % Normalize to [0, 1] range

% Display the recovered images
figure;
imshow(low_freq_recovered, []);
title('Recovered Image from Low Frequencies');

figure;
imshow(high_freq_recovered, []);
title('Recovered Image from High Frequencies');

% imwrite(high_freq_magnitude_norm, 'high_freq_spectrum.png');  
% imwrite(low_freq_magnitude_norm, 'low_freq_spectrum.png');  
% imwrite(low_freq_recovered, 'blend_0.5_low.png'); 
% imwrite(high_freq_recovered, 'blend_0.5_high.png');  