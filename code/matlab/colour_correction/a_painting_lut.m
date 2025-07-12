clear; close all;

% painting_path = '/home/oem/eliza/mac-shared/registered_output/yoda_halogen_fuji_exp0_before/yoda_halogen_fuji_exp0_and_yoda_before_registered.mat';
% cc_ref_path = '/home/oem/eliza/mac-shared/colorchecker/reference_cc.mat';
% cc_film_path = '/home/oem/eliza/mac-shared/colorchecker/yoda_halogen_fuji_exp0_colorchecker.mat';


painting_after_path = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
painting_before_path = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
film_path = '/home/oem/eliza/data/xyz_lab_rgb/film/cactus_led_fuji_exp0.mat';
cc_film_path = '/home/oem/eliza/data/xyz_lab_rgb/colorchecker/cactus_led_fuji_exp0_colorchecker.mat';
cc_ref_path = '/home/oem/eliza/data/xyz_lab_rgb/reference/xrite_cc_reference_official.mat';

painting_after = load(painting_after_path);
painting_before = load(painting_before_path);
cc_ref = load(cc_ref_path);
cc_film = load(cc_film_path);
film_data = load(film_path);


% painting_data = load(painting_path);
% cc_ref = load(cc_ref_path);
% cc_film = load(cc_film_path);

xyz_ref = cc_ref.XYZ_ref;
xyz_film = cc_film.patchXYZ;

painting_film_xyz = film_data.XYZ_img;   % h x w x 3
painting_ref_xyz = painting_before.XYZ_img;
[h, w, ~] = size(painting_film_xyz);

%% --- Convert colorchecker patches to ProPhoto RGB (range [0,1]) ---
% Use double precision and ProPhoto RGB conversion
% rgb_film = xyz2rgb(xyz_film ./ 100, 'ColorSpace', 'prophoto-rgb');
% rgb_ref  = xyz2rgb(xyz_ref  ./ 100, 'ColorSpace', 'prophoto-rgb');
rgb_film = cc_film.patchRGB;
rgb_ref = cc_ref.RGB_ref;

%% --- Build scattered interpolant LUT: rgb_film → rgb_ref ---
F_R = scatteredInterpolant(rgb_film(:,1), rgb_film(:,2), rgb_film(:,3), rgb_ref(:,1), 'linear', 'nearest');
F_G = scatteredInterpolant(rgb_film(:,1), rgb_film(:,2), rgb_film(:,3), rgb_ref(:,2), 'linear', 'nearest');
F_B = scatteredInterpolant(rgb_film(:,1), rgb_film(:,2), rgb_film(:,3), rgb_ref(:,3), 'linear', 'nearest');

%% --- Convert painting measured XYZ to ProPhoto RGB ---
% painting_film_rgb = xyz2rgb(reshape(painting_film_xyz, [], 3) ./ 100, 'ColorSpace', 'prophoto-rgb'); % (h*w) x 3
painting_film_rgb = film_data.RGB_img;
%% --- Apply LUT to painting RGB pixelwise ---
painting_film_rgb_flat = reshape(painting_film_rgb, [], 3);  % (H*W) x 3
rgb_corr_flat = zeros(size(painting_film_rgb_flat));
for i = 1:size(painting_film_rgb_flat,1)
    rgb_corr_flat(i,1) = F_R(painting_film_rgb_flat(i,1), painting_film_rgb_flat(i,2), painting_film_rgb_flat(i,3));
    rgb_corr_flat(i,2) = F_G(painting_film_rgb_flat(i,1), painting_film_rgb_flat(i,2), painting_film_rgb_flat(i,3));
    rgb_corr_flat(i,3) = F_B(painting_film_rgb_flat(i,1), painting_film_rgb_flat(i,2), painting_film_rgb_flat(i,3));
end
rgb_corr = reshape(rgb_corr_flat, h, w, 3);


%% --- Optionally convert corrected RGB back to XYZ for DeltaE analysis ---
% Use xyz2rgb's inverse: rgb2xyz (if you have it), else evaluate DeltaE in RGB
% if exist('rgb2xyz', 'file')
%     corrected_painting_xyz = rgb2xyz(rgb_corr, 'ColorSpace', 'prophoto-rgb') * 100; % convert back to XYZ, in same scale as reference
% else
%     warning('rgb2xyz function not found. DeltaE analysis in Lab will not be computed.');
% end

corrected_painting_xyz = prophoto2xyz(rgb_corr, true);

%% --- DeltaE2000 evaluation ---
if exist('rgb2xyz', 'file')
    painting_ref_xyz_reshaped = reshape(painting_ref_xyz, [], 3);
    painting_corr_xyz_reshaped = reshape(corrected_painting_xyz, [], 3);
    Lab_corrected = xyz2lab_custom(painting_corr_xyz_reshaped);
    Lab_reference = xyz2lab_custom(painting_ref_xyz_reshaped);
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
% output_path = 'corrected_prophoto_rgbLUT.tif';
% saveProPhotoTIFF(rgb_corr_img, output_path);
fprintf('Saved ProPhoto TIFF: %s\n', output_path);
