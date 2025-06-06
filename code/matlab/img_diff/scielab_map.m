clear; 

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



% Convert to XYZ
XYZ1 = rgb2xyz(img1, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');
XYZ2 = rgb2xyz(img2, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');

% Set s-CIELAB parameters
samples_per_degree = 60;
whitepoint = [0.9642, 1.0000, 0.8251];  % D50

% Compute s-CIELAB map
dE_map = scielab(samples_per_degree, XYZ1, XYZ2, whitepoint, 'XYZ');
%%
% Display results
figure;
subplot(1,3,1); imshow(img1); title('Image 1'); axis off;
subplot(1,3,2); imshow(img2); title('Image 2'); axis off;
subplot(1,3,3); imagesc(dE_map); title('s-CIELAB ΔE'); axis image off;
colormap jet; colorbar;

% Binary mask with Otsu threshold
change_mask = mat2gray(dE_map) > graythresh(mat2gray(dE_map));
figure; imshow(change_mask); title('Change Mask');

%%
dE_vector = dE_map(:);


change_percentile = 70;

change_threshold = prctile(dE_vector, change_percentile);
change_mask = dE_map > change_threshold;

figure; imshow(change_mask); 

