img_path = "/Volumes/School/Thesis/thesis-repo/data/scream/mod/scream_verm_nonunif2.png";
rgb_image = imread(img_path);

%% Fourier-Based Enhancement Inspired by FourLLIE
% Parameters
D0 = 10;        % Frequency threshold for low frequencies
gamma = 1.5;    % Amplification factor for low frequencies

% Create a zero matrix for the corrected RGB image
correctedImgRGB = zeros(size(rgb_image), 'uint16');  % Initialize as uint16

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
    
    % Normalize the enhanced channel to the [0, 1] range
    channel_enhanced = mat2gray(channel_enhanced);
    
    % Convert the normalized channel back to uint16 (range [0, 65535])
    correctedImgRGB(:, :, channel) = uint16(channel_enhanced * 65535);
    
    % Debugging: Display intermediate steps for this channel
    figure;
    subplot(2,3,1); imshow(mat2gray(channelImg)); title(sprintf('Original Channel %d', channel));
    subplot(2,3,2); imshow(mat2gray(amplitude)); title('Original Amplitude');
    subplot(2,3,3); imshow(mat2gray(amplitude_enhanced)); title('Enhanced Amplitude');
    subplot(2,3,4); imshow(mat2gray(phase)); title('Phase');
    subplot(2,3,5); imshow(channel_enhanced); title('Enhanced Channel');
    subplot(2,3,6); imshow(mat2gray(D)); title('Distance Map');
end

%% Final Display
figure;
subplot(1,2,1); imshow(rgb_image); title('Original RGB Image');
subplot(1,2,2); imshow(correctedImgRGB); title('Fourier Enhanced RGB Image');

%% Save the final enhanced image as uint16
% imwrite(correctedImgRGB, 'scream_verm_high_corr.png');  % Saving as uint16 image

%%

img_path = "/Volumes/School/Thesis/thesis-repo/code/matlab/img_diff/results/deltaE/scream_low_exp_deltaE_map.png";
diff_img = imread(img_path);

F = fft2(diff_img);
F_shifted = fftshift(F);

magnitude_spectrum = mat2gray(log(abs(F_shifted)+1));

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
 %Apply inverse FFT to get the spatial domain images 
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


% imwrite(uint8(255 * magnitude_spectrum), 'fourier_spectrum.png');


    