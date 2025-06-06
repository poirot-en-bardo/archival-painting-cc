% Read the images
% img_path1 = "/Volumes/School/Thesis/thesis-repo/data/scream/mod/scream_verm_exp.png";
% img_path2 = "/Volumes/School/Thesis/thesis-repo/data/scream/mod/scream_low_exp.png";

img_path1 = '/Volumes/School/Thesis/data/captures/registered/prophoto/yoda1_ref_hsi_after.png';
img_path2 = '/Volumes/School/Thesis/data/captures/registered/prophoto/yoda1_reg_hsi_before.png';
% img_path1 = '/Volumes/School/Thesis/data/captures/registered/prophoto/cactus1_ref_hsi_after.png';
% img_path2 = '/Volumes/School/Thesis/data/captures/registered/prophoto/cactus1_reg_hsi_before.png';
% img_path1 = '/Volumes/School/Thesis/data/captures/registered/prophoto/yoda3_reg_kodak_halogen_before.png';
% img_path2 = '/Volumes/School/Thesis/data/captures/registered/prophoto/yoda3_ref_hsi_kodak_halogen_before.png';
% %
save_folder = './results/deltaE';

img1 = imread(img_path1);
img2 = imread(img_path2);
imshowpair(img1, img2);


% Ensure images are in uint16 format
% img1 = uint16(img1);
% img2 = uint16(img2);

% Scale images to [0, 1] range for rgb2lab conversion
img1_norm = double(img1) / double(intmax('uint16'));
img2_norm = double(img2) / double(intmax('uint16'));

% Exposure normalization
% mean1 = mean(img1_norm(:));
% mean2 = mean(img2_norm(:));
% img2_norm = img2_norm * (mean1 / mean2);

% Convert RGB images to L*a*b* color space
lab1 = rgb2lab(img1_norm, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');
lab2 = rgb2lab(img2_norm, 'ColorSpace', 'prophoto-rgb', 'WhitePoint', 'd50');

% Compute ΔE 2000 differences
deltaE_vector = deltaE2000(reshape(lab1, [], 3), reshape(lab2, [], 3));

% Reshape the ΔE 2000 differences to match the original image dimensions
[rows, cols, ~] = size(img1);
deltaE_map = reshape(deltaE_vector, rows, cols);


% Normalize the ΔE map to [0, 1] range
deltaE_map_normalized = mat2gray(deltaE_map);  % Converts to [0, 1]
% Optionally, scale the ΔE map to uint16 (range [0, 65535])
deltaE_map_uint16 = uint16(deltaE_map_normalized * 65535);

% Save as a grayscale image (uint16)
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
[~, img_name, ~] = fileparts(img_path2);
% imwrite(deltaE_map_uint16, fullfile(save_folder, img_name + '_deltaE_map.png'));  


% or
% deltaE_map = imcolordiff(lab1, lab2, 'Standard', 'CIEDE2000');
%%
% Display the ΔE 2000 difference as a heat map
figure;
imagesc(deltaE_map);
clim([min(deltaE_map(:)) 10]);
colorbar;
title('ΔE00 Difference Heat Map');
axis image off;
colormap jet;

% Save figure
[~, img_name, ~] = fileparts(img_path1);
% saveas(gcf, fullfile(save_folder, [img_name '_deltaE_norm.png']));
%%

deltaE_map_norm = mat2gray(deltaE_map);

% Scale to 16-bit grayscale
deltaE_uint16 = uint16(deltaE_map_norm * 65535);

% Save
% imwrite(deltaE_uint16, fullfile(save_folder, [img_name '_deltaE.png']));

%%
% Choose a threshold for ΔE00 difference
thresh = 5;  % Try values like 2, 3, 4, 5

% Option 1: Apply threshold to raw ΔE map (recommended)
mask = deltaE_map > thresh;

% Option 2: Apply threshold to normalized ΔE map (not typical, but possible)
% mask = deltaE_map_normalized > (thresh / max(deltaE_map(:))); % Not needed if using original units

% Optional: Automatic thresholding (Otsu's method, if you want)
level = graythresh(mat2gray(deltaE_map));
mask_otsu = mat2gray(deltaE_map) > level;

% Display the mask
figure;
subplot(1,2,1); imshow(mask); title(['ΔE00 > ' num2str(thresh) ' (Manual Threshold)']);
subplot(1,2,2); imshow(mask_otsu); title('ΔE00 Mask (Otsu)');
%%

dE_vector = deltaE_map(:);


change_percentile = 70;

change_threshold = prctile(dE_vector, change_percentile);
change_mask = deltaE_map > change_threshold;


% otsu's threshold
% level = graythresh(mat2gray(dE_map));  % Returns a threshold in [0,1]
% change_mask = mat2gray(dE_map) > level;


figure; imshow(change_mask); 

