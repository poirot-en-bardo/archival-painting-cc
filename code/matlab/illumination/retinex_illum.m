% Read the image
path1 = '/Volumes/School/Thesis/thesis-repo/data/scream/mod/scream_verm_nonunif.png';
img = imread(path1);

img = im2double(img); % Convert to double for processing

% Define Retinex parameters
sigma_list = [15, 80, 250]; % Multi-scale Gaussian sigmas
weight = 1 / numel(sigma_list); % Equal weight for each scale

% Initialize the output image
[m, n, c] = size(img);
MSR = zeros(m, n, c);

% Compute Retinex for each color channel
for k = 1:c
    retinex = zeros(m, n);
    for sigma = sigma_list
        % Gaussian filter
        gauss_filter = fspecial('gaussian', [m n], sigma);
        img_blur = imfilter(img(:,:,k), gauss_filter, 'replicate');
        
        % Compute Retinex: log(I) - log(I * blur)
        retinex = retinex + weight * (log(img(:,:,k) + eps) - log(img_blur + eps));
    end
    % Normalize to [0, 1] range
    MSR(:,:,k) = mat2gray(retinex);
end

% Display the result
figure;
imshow(MSR);
title('Multi-Scale Retinex Result');
