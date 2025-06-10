clc; clear; close all;

%% Load spectral data (same as before)
ill = importdata('../../../data/CIE_D65.txt'); % Illuminant
CMFs_1931 = importdata('../../../data/CIE2degCMFs_1931.txt'); % CMFs
CMFs = CMFs_1931; % Color Matching Functions

% Set random seed for reproducibility
rng(10);

%% Select the reference and the cube to be corrected
cubeFile = "../../../data/colorChecker_SG/cubes/cubeCC_no6-ekta100-expo3.hdr";
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

%% Interpolate illuminant and CMFs for the input cube
illIP = interp1(ill(:,1), ill(:,2), bands, 'spline');  % Interpolate the illuminant
CMFsIP = [interp1(CMFs(:,1), CMFs(:,2), bands, 'spline'), ...  % X channel
          interp1(CMFs(:,1), CMFs(:,3), bands, 'spline'), ...  % Y channel
          interp1(CMFs(:,1), CMFs(:,4), bands, 'spline')];    % Z channel

% Compute tristimulus values for the input cube
sp_tristREF = CMFsIP .* illIP;

% Compute XYZ values for the input cube
xyz_input = (lincube * sp_tristREF) ./ sum(sp_tristREF(:,2), 1);

% Convert XYZ to RGB (ProPhoto-RGB)
rgb_input = xyz2rgb(xyz_input, 'ColorSpace', 'prophoto-rgb');
rgb_image = reshape(rgb_input, [m, n, 3]);

disp(['Input image size: ', num2str(size(rgb_image))]);

%% Fourier-Based Enhancement Inspired by FourLLIE

path1 = '/Volumes/School/Thesis/data/captures/registered/prophoto/cactus3_reg_fuji_led_after.png';
% Read and convert image to HSV
img = imread(path1); % uint16 image
rgb_image = im2double(img); % Convert to [0, 1] range for processing
%%
% Parameters
D0 = 5;        % Frequency threshold for low frequencies
gamma = 1.5;    % Amplification factor for low frequencies

correctedImgRGB = zeros(size(rgb_image));

for channel = 1:3
    % Extract current channel
    channelImg = rgb_image(:, :, channel);
    
    % Compute the Fourier transform and shift DC to center
    F = fft2(channelImg);
    Fshift = fftshift(F);
    
    % Decompose into amplitude and phase
    amplitude = abs(Fshift);
    phase = angle(Fshift);
    
    % Create distance map from center
    [Mch, Nch] = size(channelImg);
    [x, y] = meshgrid(1:Nch, 1:Mch);
    D = sqrt((x - Nch/2).^2 + (y - Mch/2).^2);
    
    % Create a low-frequency mask: where frequency is below threshold D0
    lowFreqMask = (D < D0);
    
    % Enhance the amplitude in low frequencies
    amplitude_enhanced = amplitude;
    amplitude_enhanced(lowFreqMask) = amplitude(lowFreqMask) * gamma;
    
    % Reconstruct the enhanced Fourier spectrum
    F_enhanced = amplitude_enhanced .* exp(1i * phase);
    F_enhanced = ifftshift(F_enhanced);
    
    % Inverse Fourier transform to get the enhanced channel
    channel_enhanced = real(ifft2(F_enhanced));
    
    % Normalize the enhanced channel to [0,1]
    channel_enhanced = mat2gray(channel_enhanced);
    
    % Save enhanced channel
    correctedImgRGB(:, :, channel) = channel_enhanced;
    
    % Debugging: Display intermediate steps for this channel
    % figure;
    % subplot(2,3,1); imshow(mat2gray(channelImg)); title(sprintf('Original Channel %d', channel));
    % subplot(2,3,2); imshow(mat2gray(amplitude)); title('Original Amplitude');
    % subplot(2,3,3); imshow(mat2gray(amplitude_enhanced)); title('Enhanced Amplitude');
    % subplot(2,3,4); imshow(mat2gray(phase)); title('Phase');
    % subplot(2,3,5); imshow(channel_enhanced); title('Enhanced Channel');
    % subplot(2,3,6); imshow(mat2gray(D)); title('Distance Map');
end

% Final Display
figure;
subplot(1,2,1); imshow(rgb_image); title('Original RGB Image');
subplot(1,2,2); imshow(correctedImgRGB); title('Fourier Enhanced RGB Image');
