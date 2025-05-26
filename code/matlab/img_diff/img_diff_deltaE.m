% Read the images
% img_path1 = "/Volumes/School/Thesis/thesis-repo/data/scream/mod/scream_verm_exp.png";
% img_path2 = "/Volumes/School/Thesis/thesis-repo/data/scream/mod/scream_low_exp.png";
img_path1 = '/Volumes/School/Thesis/thesis-repo/code/matlab/illumination/scream_low_corr.png';
img_path2 = '/Volumes/School/Thesis/thesis-repo/code/matlab/illumination/scream_verm_high_corr.png';
img_path1 = '/Volumes/School/Thesis/thesis-repo/code/matlab/img_diff/results/prophoto/icc_manual_reference_ProPhoto.png';
img_path2  = '/Volumes/School/Thesis/thesis-repo/code/matlab/img_diff/results/prophoto/icc_manual_registered_ProPhoto.png';
save_folder = './results/deltaE';

img1 = imread(img_path1);
img2 = imread(img_path2);


% Ensure images are in uint16 format
% img1 = uint16(img1);
% img2 = uint16(img2);

% Scale images to [0, 1] range for rgb2lab conversion
img1_norm = double(img1) / double(intmax('uint16'));
img2_norm = double(img2) / double(intmax('uint16'));

% Exposure normalization
mean1 = mean(img1_norm(:));
mean2 = mean(img2_norm(:));
% img2_norm = img2_norm * (mean1 / mean2);

% Convert RGB images to L*a*b* color space
lab1 = rgb2lab(img1_norm, 'ColorSpace', 'prophoto-rgb');
lab2 = rgb2lab(img2_norm, 'ColorSpace', 'prophoto-rgb');

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
clim([0 20]);
colorbar;
title('ΔE 2000 Difference Heat Map');
axis image off;

% Save figure
[~, img_name, ~] = fileparts(img_path2);
% saveas(gcf, fullfile(save_folder, [img_name '_deltaE_norm.png']));

