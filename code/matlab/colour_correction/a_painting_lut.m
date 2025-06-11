clear; close all;

painting_path = '/home/oem/eliza/mac-shared/registered_output/yoda_halogen_fuji_exp0_before/yoda_halogen_fuji_exp0_and_yoda_before_registered.mat';
cc_ref_path = '/home/oem/eliza/mac-shared/colorchecker/reference_cc.mat';
cc_film_path = '/home/oem/eliza/mac-shared/colorchecker/yoda_halogen_fuji_exp0_colorchecker.mat';

painting_data = load(painting_path);
cc_ref = load(cc_ref_path);
cc_film = load(cc_film_path);

xyz_ref = cc_ref.patchXYZ;
xyz_film = cc_film.patchXYZ;

painting_film_xyz = painting_data.reg_xyz;   % h x w x 3
painting_ref_xyz = painting_data.ref_cropped;
[h, w, ~] = size(painting_film_xyz);

%% --- Convert colorchecker patches to ProPhoto RGB (range [0,1]) ---
% Use double precision and ProPhoto RGB conversion
rgb_film = xyz2rgb(xyz_film ./ 100, 'ColorSpace', 'prophoto-rgb');
rgb_ref  = xyz2rgb(xyz_ref  ./ 100, 'ColorSpace', 'prophoto-rgb');

%% --- Build scattered interpolant LUT: rgb_film → rgb_ref ---
F_R = scatteredInterpolant(rgb_film(:,1), rgb_film(:,2), rgb_film(:,3), rgb_ref(:,1), 'linear', 'nearest');
F_G = scatteredInterpolant(rgb_film(:,1), rgb_film(:,2), rgb_film(:,3), rgb_ref(:,2), 'linear', 'nearest');
F_B = scatteredInterpolant(rgb_film(:,1), rgb_film(:,2), rgb_film(:,3), rgb_ref(:,3), 'linear', 'nearest');

%% --- Convert painting measured XYZ to ProPhoto RGB ---
painting_film_rgb = xyz2rgb(reshape(painting_film_xyz, [], 3) ./ 100, 'ColorSpace', 'prophoto-rgb'); % (h*w) x 3

%% --- Apply LUT to painting RGB pixelwise ---
rgb_corr = zeros(size(painting_film_rgb));
for i = 1:size(painting_film_rgb,1)
    rgb_corr(i,1) = F_R(painting_film_rgb(i,1), painting_film_rgb(i,2), painting_film_rgb(i,3));
    rgb_corr(i,2) = F_G(painting_film_rgb(i,1), painting_film_rgb(i,2), painting_film_rgb(i,3));
    rgb_corr(i,3) = F_B(painting_film_rgb(i,1), painting_film_rgb(i,2), painting_film_rgb(i,3));
end

%% --- Optionally convert corrected RGB back to XYZ for DeltaE analysis ---
% Use xyz2rgb's inverse: rgb2xyz (if you have it), else evaluate DeltaE in RGB
if exist('rgb2xyz', 'file')
    corrected_painting_xyz = rgb2xyz(rgb_corr, 'ColorSpace', 'prophoto-rgb') * 100; % convert back to XYZ, in same scale as reference
else
    warning('rgb2xyz function not found. DeltaE analysis in Lab will not be computed.');
end

%% --- DeltaE2000 evaluation ---
if exist('rgb2xyz', 'file')
    painting_ref_xyz_reshaped = reshape(painting_ref_xyz, [], 3);
    Lab_corrected = xyz2lab(corrected_painting_xyz);
    Lab_reference = xyz2lab(painting_ref_xyz_reshaped);
    deltaE = deltaE2000(Lab_reference, Lab_corrected);
    deltaE_img = reshape(deltaE, h, w);
    figure; imagesc(deltaE_img); axis image; colorbar; caxis([0 10]);
    title('\DeltaE_{2000} (RGB LUT correction)');
    fprintf('Mean ΔE2000: %.3f\n', mean(deltaE(:)));
    fprintf('Median ΔE2000: %.3f\n', median(deltaE(:)));
    fprintf('95th percentile ΔE2000: %.3f\n', prctile(deltaE, 95));
end

%% --- Visualize result ---
rgb_corr_img = reshape(rgb_corr, h, w, 3);
figure; imshow(rgb_corr_img); title('Painting Corrected (ProPhoto RGB, RGB LUT)');

%% --- Save as 16-bit TIFF ---
output_path = 'corrected_prophoto_rgbLUT.tif';
% saveProPhotoTIFF(rgb_corr_img, output_path);
fprintf('Saved ProPhoto TIFF: %s\n', output_path);
