clear; close all;

addpath('SCIELAB-1996');
savepath;

% Load XYZ images directly
img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus2_ref_hsi_fuji_led_after.png'));
img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus2_reg_fuji_led_after.png'));
img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda2_ref_hsi_kodak_halogen_after.png'));
img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda2_reg_kodak_halogen_after.png'));
% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/icc_manual_reference_ProPhoto.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/icc_manual_registered_ProPhoto.png'));
% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus1_ref_hsi_after.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus1_reg_hsi_before.png'));
% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda1_ref_hsi_after.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda1_reg_hsi_before.png'));

XYZ1 = rgb2xyz(img1, 'ColorSpace', 'prophoto-rgb', 'WhitePoint','d50');
XYZ2 = rgb2xyz(img2, 'ColorSpace', 'prophoto-rgb', 'WhitePoint','d50');

% % Ensure images are 3-channel XYZ
% if size(XYZ1, 3) ~= 3 || size(XYZ2, 3) ~= 3
%     error('Input images must be 3-channel XYZ data.');
% end

% Multi-scale Retinex parameters
sigma_values = [1, 2, 20];

% --- Path 1: Direct deltaE2000 ---
XYZ1_smooth = imgaussfilt(XYZ1, 0.5);
XYZ2_smooth = imgaussfilt(XYZ2, 0.5);
Lab1_direct = xyz2lab(XYZ1_smooth, [0.9642, 1.0000, 0.8251]);
Lab2_direct = xyz2lab(XYZ2_smooth, [0.9642, 1.0000, 0.8251]);
Lab1_col = reshape(Lab1_direct, [], 3);
Lab2_col = reshape(Lab2_direct, [], 3);
dE_map_direct = reshape(deltaE2000(Lab1_col, Lab2_col), size(Lab1_direct,1), size(Lab1_direct,2));

% --- Path 2: Retinex + deltaE2000 ---
XYZ1_norm = robust_normalize(XYZ1);
XYZ2_norm = robust_normalize(XYZ2);
refl1 = apply_multiscale_retinex(XYZ1_norm, sigma_values);
refl2 = apply_multiscale_retinex(XYZ2_norm, sigma_values);
dE_map_retinex = compute_reflectance_deltaE_from_XYZ(refl1, refl2);

% --- Selection based on contrast metric ---
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
dE_map = dE_map_retinex;
dE_map = dE_map_direct;

% Visualization
figure;
subplot(1,3,1); imshow(XYZ1); title('Image 1'); axis off;
subplot(1,3,2); imshow(XYZ2); title('Image 2'); axis off;
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

function reflectance = apply_multiscale_retinex(XYZ_img, sigma_values)
    XYZ_img(XYZ_img <= 0) = eps;
    reflectance_sum = zeros(size(XYZ_img));
    for sigma = sigma_values
        illum = imgaussfilt(rgb2gray(XYZ_img), sigma);
        illum(illum <= 0) = eps;
        for c = 1:3
            reflectance_sum(:,:,c) = reflectance_sum(:,:,c) + log(XYZ_img(:,:,c)) - log(illum);
        end
    end
    reflectance = real(reflectance_sum) / numel(sigma_values);
    reflectance = reflectance - min(reflectance(:));
    reflectance = reflectance / max(reflectance(:));
end

function dE_map = compute_reflectance_deltaE_from_XYZ(XYZ1, XYZ2)
    d50 = [0.9642, 1.0000, 0.8251];
    Lab1 = xyz2lab(XYZ1, d50);
    Lab2 = xyz2lab(XYZ2, d50);
    Lab1_col = reshape(Lab1, [], 3);
    Lab2_col = reshape(Lab2, [], 3);
    dE_map = reshape(deltaE2000(Lab1_col, Lab2_col), size(Lab1,1), size(Lab1,2));
end
