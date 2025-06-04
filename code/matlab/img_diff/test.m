img1 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/yoda_reg_hsi_before.png'));
img2 = im2double(imread('/Volumes/School/Thesis/data/captures/registered/yoda_ref_hsi_after.png'));


gray1 = rgb2gray(img1);
gray2 = rgb2gray(img2);

[ssim_index, ssim_map] = ssim(gray1, gray2);
unchanged_mask = ssim_map > 0.9;  

% Visualize
figure; imshow(unchanged_mask); title('Unchanged Regions');
