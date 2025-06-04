clear; close all;

addpath('SCIELAB-1996');
savepath; 



img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus2_ref_hsi_fuji_led_after.png'));
img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus2_reg_fuji_led_after.png'));
% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda2_ref_hsi_kodak_halogen_after.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda2_reg_kodak_halogen_after.png'));
% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/icc_manual_reference_ProPhoto.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/icc_manual_registered_ProPhoto.png'));
% img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus1_ref_hsi_after.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/cactus1_reg_hsi_before.png'));
% img1 = aim2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda1_ref_hsi_after.png'));
% img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/prophoto/yoda1_reg_hsi_before.png'));


sigma = 1; 

d50 = [0.9642, 1.0000, 0.8251];

ssim_val = ssim(rgb2gray(img1), rgb2gray(img2));

if ssim_val > 0.7
    img2 = imhistmatch(img2, img1);
end
%%
dE_map = compute_reflectance_deltaE(img1, img2, sigma);

%%
figure;

subplot(1,3,1); imshow(img1); title('Image 1'); axis off;

subplot(1,3,2); imshow(img2); title('Image 2'); axis off;
subplot(1,3,3);
imagesc(dE_map); 
axis image off;
title('ΔE Difference');
colormap jet;  % or 'parula'
colorbar;

%%
% figure;
% imshowpair(img1, img2);
%%
% dE_map = im2double(imread('/Volumes/School/Thesis/thesis-repo/code/matlab/img_diff/results/deltaE/yoda1_ref_hsi_after_deltaE.png'));
%%

dE_vector = dE_map(:);


change_percentile = 70;

change_threshold = prctile(dE_vector, change_percentile);
change_mask = dE_map > change_threshold;


% otsu's threshold
% level = graythresh(mat2gray(dE_map));  % Returns a threshold in [0,1]
% change_mask = mat2gray(dE_map) > level;


figure; imshow(change_mask); 



%%

function reflectance = apply_retinex(img_rgb, sigma)
    
    img_rgb = im2double(img_rgb);
    img_rgb(img_rgb <= 0) = eps;

    % Illumination estimation (blurred grayscale image)
    illum = imgaussfilt(rgb2gray(img_rgb), sigma);
    illum(illum <= 0) = eps; % Avoid log(0) 

    reflectance = zeros(size(img_rgb));
    for c = 1:3
        reflectance(:,:,c) = log(img_rgb(:,:,c)) - log(illum);
    end

    reflectance = real(reflectance);

    % Normalise reflectance to [0,1]
    reflectance = reflectance - min(reflectance(:));
    reflectance = reflectance / max(reflectance(:));

    % reflectance_sum = zeros(size(img_rgb));

    % for s = sigma
    %     illum = imgaussfilt(rgb2gray(img_rgb), s);
    %     illum(illum <= 0) = eps;
    % 
    %     tmp_refl = zeros(size(img_rgb));
    %     for c = 1:3
    %         tmp_refl(:,:,c) = log(img_rgb(:,:,c)) - log(illum);
    %     end
    %     reflectance_sum = reflectance_sum + real(tmp_refl);
    % end
    % 
    % % Average and normalize to [0,1]
    % reflectance = reflectance_sum / numel(sigma);
    % reflectance = reflectance - min(reflectance(:));
    % reflectance = reflectance / max(reflectance(:));
end


function dE_map = compute_reflectance_deltaE(img1, img2, sigma)
    % Apply Retinex to both images
    % sigma = [1, 50, 400];

    % grayDiff = abs(rgb2gray(img1) - rgb2gray(img2));
    % [~, S] = imgradient(grayDiff); % gradient magnitude
    % gradVal = median(S(:));  % typical edge strength (more robust than mean to noise)
    % scaleFactor = 1;
    % sigma1 = max(1, round(0.2 * gradVal * scaleFactor));
    % sigma2 = max(1, round(0.5 * gradVal * scaleFactor));
    % sigma3 = max(1, round(1.2 * gradVal * scaleFactor));
    % sigma = unique([sigma1, sigma2, sigma3]);

  
    refl1 = apply_retinex(img1, sigma);
    refl2 = apply_retinex(img2, sigma);

    
    d50 = [0.9642, 1.0000, 0.8251];

    % figure; imshow(refl1,[]); title('Reflectance 1');
    % figure; imshow(refl2.^0.8, []); title('Reflectance 2');
    % 
    % Convert reflectance RGB to XYZ
    XYZ1 = rgb2xyz(refl1, 'ColorSpace', 'prophoto-rgb', 'WhitePoint','d50');
    XYZ2 = rgb2xyz(refl2, 'ColorSpace', 'prophoto-rgb', 'WhitePoint','d50');
    
    % % Convert XYZ to Lab
    Lab1 = xyz2lab(XYZ1, d50);
    Lab2 = xyz2lab(XYZ2, d50);
    % 
    % % Compute deltaE2000
    % dE_map = deltaE(Lab1, Lab2, isInputLab=true);
    % dE_map = deltaE(refl1, refl2);

    % dE_map = imcolordiff(Lab1,Lab2,"Standard","CIEDE2000", 'isInputLab', true);
    % 
    sz = size(Lab1);
    Lab1_col = reshape(Lab1, [], 3);
    Lab2_col = reshape(Lab2, [], 3);

    % Compute deltaE2000 (vectorized)
    DE_col = deltaE2000(Lab1_col, Lab2_col);

    % Reshape ΔE back to original image shape
    dE_map = reshape(DE_col, sz(1), sz(2));
    
end


function dE_map = compute_reflectance_deltaE_sCIELAB(img1, img2, sigma)
    % Apply Retinex
    refl1 = apply_retinex(img1, sigma);
    refl2 = apply_retinex(img2, sigma);

    % Convert to XYZ using ProPhoto RGB with D50
    XYZ1 = rgb2xyz(refl1, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');
    XYZ2 = rgb2xyz(refl2, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');

    % Define viewing conditions
    samples_per_degree = 60;  % adjust based on viewing setup
    whitepoint = [0.9642, 1.0000, 0.8251];  % D50

    % Compute spatial color difference
    dE_map = scielab(samples_per_degree, XYZ1, XYZ2, whitepoint, 'XYZ');
end

