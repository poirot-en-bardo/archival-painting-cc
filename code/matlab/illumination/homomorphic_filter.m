% Read and convert image to HSV
path1 = '/Volumes/School/Thesis/thesis-repo/data/scream/mod/scream_verm_nonunif.png';
path2 = '/Volumes/School/Thesis/thesis-repo/data/scream/mod/scream_verm_nonunif2.png';

path1 = '/Volumes/School/Thesis/data/captures/registered/prophoto/yoda2_reg_kodak_halogen_after.png';

% Read and convert image to HSV
img = imread(path1); % uint16 image
img = im2double(img); % Convert to [0, 1] range for processing
hsv_img = rgb2hsv(img);
V = hsv_img(:,:,3); % Extract brightness (Value) channel

% Log transformation to separate illumination and reflectance
log_V = log(double(V) + eps); % Avoid log(0) issues with small constant

% Fourier Transform
F = fft2(log_V);
F_shifted = fftshift(F);

% Gaussian low-pass filter for illumination estimation
[rows, cols] = size(V);
[X, Y] = meshgrid(1:cols, 1:rows);
centerX = ceil(cols / 2);
centerY = ceil(rows / 2);
D = sqrt((X - centerX).^2 + (Y - centerY).^2);

D0 = 0.2; % Cutoff frequency 
H = exp(-(D.^2) / (2 * D0^2)); % Gaussian low-pass filter

% Apply the filter to estimate illumination component
F_filtered = F_shifted .* H;

% Inverse Fourier Transform to get illumination in log space
F_unshifted = ifftshift(F_filtered);
log_illumination = real(ifft2(F_unshifted));

% Reflectance = log(V) - log(illumination)
log_reflectance = log_V - log_illumination;

% Convert back from log space to linear space
V_corrected = exp(log_reflectance);

% Normalize to maintain mean brightness
mean_V = mean(V(:));
V_corrected = V_corrected * (mean_V / mean(V_corrected(:)));

% Ensure V_corrected is within [0, 1] range
V_corrected = max(0, min(1, V_corrected));

% Update V channel in HSV and convert back to RGB
hsv_img(:,:,3) = V_corrected;
corrected_img = hsv2rgb(hsv_img);

% Convert back to uint16 for output
corrected_img_uint16 = im2uint16(corrected_img);

% Display result
figure; imshow(corrected_img_uint16); title('Illumination Corrected Image');

output_path = '../../../results/illum_corrected1.png';

% Save the uint16 corrected image
% imwrite(corrected_img_uint16, output_path);
