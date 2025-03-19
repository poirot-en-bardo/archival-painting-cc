% Read the images
img_path1 = '../../data/scream/scream_rgb.png';
% img_path2 = ['../../data/scream/mod/scream_high_exp.png'];
img_path2 = '../../results/illum_corrected1.png';
save_folder = '../../results/deltaE';

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
img2_norm = img2_norm * (mean1 / mean2);

% Convert RGB images to L*a*b* color space
lab1 = rgb2lab(img1_norm, 'ColorSpace', 'prophoto-rgb');
lab2 = rgb2lab(img2_norm, 'ColorSpace', 'prophoto-rgb');

% Compute ΔE 2000 differences
deltaE_vector = deltaE2000(reshape(lab1, [], 3), reshape(lab2, [], 3));

% Reshape the ΔE 2000 differences to match the original image dimensions
[rows, cols, ~] = size(img1);
deltaE_map = reshape(deltaE_vector, rows, cols);

% or
% deltaE_map = imcolordiff(lab1, lab2, 'Standard', 'CIEDE2000');

% Display the ΔE 2000 difference as a heat map
figure;
imagesc(deltaE_map);
% clim([0 10]);
colorbar;
title('ΔE 2000 Difference Heat Map');
axis image off;

% Save figure
[~, img_name, ~] = fileparts(img_path2);
saveas(gcf, fullfile(save_folder, [img_name '_deltaE_norm.png']));

