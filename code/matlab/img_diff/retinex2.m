clear; close all;

addpath('SCIELAB-1996');
savepath;

img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus2_ref_hsi_fuji_led_after.png'));
img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus2_reg_fuji_led_after.png'));
% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda2_ref_hsi_kodak_halogen_after.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda2_reg_kodak_halogen_after.png'));
% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus1_ref_hsi_after.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus1_reg_hsi_before.png'));
% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda1_ref_hsi_after.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda1_reg_hsi_before.png'));
% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda4_reg_fuji_led_after.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda4_ref_hsi_fuji_led_after.png'));

% Multi-scale Retinex parameters
sigma_values = [1, 2, 5];

% --- Path 1: Direct deltaE2000 ---
img1_smooth = imgaussfilt(img1, 0.5);
img2_smooth = imgaussfilt(img2, 0.5);
XYZ1_direct = rgb2xyz(img1_smooth, 'ColorSpace', 'prophoto-rgb', 'WhitePoint','d50');
XYZ2_direct = rgb2xyz(img2_smooth, 'ColorSpace', 'prophoto-rgb', 'WhitePoint','d50');
Lab1_direct = xyz2lab(XYZ1_direct, [0.9642, 1.0000, 0.8251]);
Lab2_direct = xyz2lab(XYZ2_direct, [0.9642, 1.0000, 0.8251]);
Lab1_col = reshape(Lab1_direct, [], 3);
Lab2_col = reshape(Lab2_direct, [], 3);
dE_map_direct = reshape(deltaE2000(Lab1_col, Lab2_col), size(Lab1_direct,1), size(Lab1_direct,2));

% --- Path 2: Retinex + deltaE2000 ---
img1_norm = robust_normalize(img1);
img2_norm = robust_normalize(img2);
refl1 = apply_multiscale_retinex(img1_norm, sigma_values);
refl2 = apply_multiscale_retinex(img2_norm, sigma_values);
dE_map_retinex = compute_reflectance_deltaE_direct(refl1, refl2);

% --- Selection based on standard deviation (contrast measure) ---
% Use high-percentile contrast measure
retinex_contrast = prctile(dE_map_retinex(:), 90) - median(dE_map_retinex(:));
direct_contrast  = prctile(dE_map_direct(:), 90) - median(dE_map_direct(:));

if retinex_contrast > direct_contrast
    dE_map = dE_map_retinex;
    disp('Using Retinex-based ΔE map');
else
    dE_map = dE_map_direct;
    disp('Using Direct ΔE map');
end

%%
dE_map = dE_map_direct;
dE_map = dE_map_retinex;
% Visualization
figure;
subplot(1,3,1); imshow(img1); title('Image 1'); axis off;
subplot(1,3,2); imshow(img2); title('Image 2'); axis off;
subplot(1,3,3);
imagesc(dE_map);
axis image off;
title('ΔE Difference');
colormap jet;
colorbar;

% Compute change mask
change_percentile = 70;
change_threshold = prctile(dE_map(:), change_percentile);
change_mask = dE_map > change_threshold;

figure; imshow(change_mask); title('Change Mask');

%% Helper Functions
function img_out = robust_normalize(img_in)
    img_out = zeros(size(img_in));
    for c = 1:3
        channel = img_in(:,:,c);
        channel_mean = mean(channel(:));
        channel_std = std(channel(:));
        img_out(:,:,c) = (channel - channel_mean) / (channel_std + eps);
    end
    img_out = img_out - min(img_out(:));
    img_out = img_out / max(img_out(:));
end

function reflectance = apply_multiscale_retinex(img_rgb, sigma_values)
    img_rgb(img_rgb <= 0) = eps;
    reflectance_sum = zeros(size(img_rgb));
    for sigma = sigma_values
        illum = imgaussfilt(rgb2gray(img_rgb), sigma);
        illum(illum <= 0) = eps;
        for c = 1:3
            reflectance_sum(:,:,c) = reflectance_sum(:,:,c) + log(img_rgb(:,:,c)) - log(illum);
        end
    end
    reflectance = real(reflectance_sum) / numel(sigma_values);
    reflectance = reflectance - min(reflectance(:));
    reflectance = reflectance / max(reflectance(:));
end

function dE_map = compute_reflectance_deltaE_direct(refl1, refl2)
    d50 = [0.9642, 1.0000, 0.8251];
    XYZ1 = rgb2xyz(refl1, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');
    XYZ2 = rgb2xyz(refl2, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');
    Lab1 = xyz2lab(XYZ1, d50);
    Lab2 = xyz2lab(XYZ2, d50);
    Lab1_col = reshape(Lab1, [], 3);
    Lab2_col = reshape(Lab2, [], 3);
    dE_map = reshape(deltaE2000(Lab1_col, Lab2_col), size(Lab1,1), size(Lab1,2));
end
