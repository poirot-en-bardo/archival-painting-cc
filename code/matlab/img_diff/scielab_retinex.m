clear; close all;

addpath('SCIELAB-1996');
savepath; 

% Load image pair

% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus2_ref_hsi_fuji_led_after.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus2_reg_fuji_led_after.png'));
% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda2_ref_hsi_kodak_halogen_after.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda2_reg_kodak_halogen_after.png'));
% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/icc_manual_reference_ProPhoto.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/icc_manual_registered_ProPhoto.png'));
% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus1_ref_hsi_after.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus1_reg_hsi_before.png'));
img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda1_ref_hsi_after.png'));
img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda1_reg_hsi_before.png'));

% Parameters
gamma = 2.0;
d50 = [0.9642, 1.0000, 0.8251];
samples_per_degree = 100;

% Illumination-Invariant Normalization (STAR-inspired)
refl1 = chroma_structure_normalization(img1, gamma);
refl2 = chroma_structure_normalization(img2, gamma);

% Convert to XYZ for s-CIELAB
disp('Converting to XYZ...');
XYZ1 = rgb2xyz(refl1, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');
XYZ2 = rgb2xyz(refl2, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');

% Compute s-CIELAB difference
dE_map = scielab(samples_per_degree, XYZ1, XYZ2, d50, 'XYZ');

% Visualization
figure;
subplot(1,3,1); imshow(img1); title('Original Image 1'); axis off;
subplot(1,3,2); imshow(img2); title('Original Image 2'); axis off;
subplot(1,3,3); imagesc(dE_map); title('s-CIELAB \DeltaE'); axis image off; colormap jet; colorbar;

%%
% Thresholding
change_threshold = prctile(dE_map(:), 70);
change_mask = dE_map > change_threshold;
figure; imshow(change_mask); title('Change Mask');

% --- Function ---
function refl = chroma_structure_normalization(img_rgb, gamma)
    img_rgb = im2double(img_rgb);
    gray = rgb2gray(img_rgb);

    dx = imfilter(gray, [-1 1], 'replicate');
    dy = imfilter(gray, [-1; 1], 'replicate');
    grad = sqrt(dx.^2 + dy.^2);

    structure = grad .^ gamma;
    texture = grad .^ (1/gamma);

    % Approximate illumination using a Gaussian blur
    sigma = 15;
    illum = imgaussfilt(gray, sigma);
    illum(illum <= 0) = eps;

    % Reflectance estimate
    refl_gray = gray ./ illum;
    refl_gray = refl_gray - min(refl_gray(:));
    refl_gray = refl_gray / max(refl_gray(:));

    % Apply chromaticity scaling
    chroma = img_rgb ./ repmat(gray + eps, 1, 1, 3);
    refl = refl_gray .* chroma;

    % Normalize reflectance
    refl = refl - min(refl(:));
    refl = refl / max(refl(:));
end
