clear; close all;

painting_path = '/home/oem/eliza/mac-shared/registered_output/yoda_halogen_fuji_exp0_before/yoda_halogen_fuji_exp0_and_yoda_before_registered.mat';
cc_ref_path = '/home/oem/eliza/mac-shared/colorchecker/reference_cc.mat';
cc_film_path = '/home/oem/eliza/mac-shared/colorchecker/yoda_halogen_fuji_exp0_colorchecker.mat';

painting_data = load(painting_path);
cc_ref = load(cc_ref_path);
cc_film = load(cc_film_path);

xyz_ref = cc_ref.patchXYZ;
xyz_film = cc_film.patchXYZ;

painting_film_xyz = painting_data.reg_xyz;
painting_ref_xyz = painting_data.ref_cropped;

%% Colour correction in Lab
[h, w, ~] = size(painting_film_xyz);

% ========== Reshape colorchecker data ==========
xyz_film_reshaped = reshape(xyz_film, [], 3);
xyz_ref_reshaped = reshape(xyz_ref, [], 3);

% --- Convert to Lab (D50/XYZ) for patch data
Lab_film = xyz2lab(xyz_film_reshaped);
Lab_ref  = xyz2lab(xyz_ref_reshaped);

% ========== Prepare painting data ==========
painting_film_xyz_reshaped = reshape(painting_film_xyz, [], 3);
painting_ref_xyz_reshaped  = reshape(painting_ref_xyz, [], 3);

% --- Painting to Lab
Lab_painting = xyz2lab(painting_film_xyz_reshaped);

% ========== Polynomial regression in Lab ==========
X_poly = poly3_features(Lab_film);   % N×19
coeffs_lab = pinv(X_poly) * Lab_ref; % 19×3

% ========== Apply correction to painting data ==========
X_poly_painting = poly3_features(Lab_painting);  % (h*w)×19
Lab_painting_corr = X_poly_painting * coeffs_lab; % (h*w)×3

% --- Clamp any crazy values
Lab_painting_corr(~isfinite(Lab_painting_corr)) = 0;

% ========== Convert both corrected and reference to XYZ for display ==========
corrected_painting_xyz = lab2xyz(Lab_painting_corr);       % (h*w)×3
painting_ref_xyz_reshaped = reshape(painting_ref_xyz, [], 3);

% --- Delta E2000 calculation (in Lab)
Lab_reference = xyz2lab(painting_ref_xyz_reshaped);
deltaE = deltaE2000(Lab_reference, Lab_painting_corr);
deltaE_img = reshape(deltaE, h, w);

%% Figure: ΔE map
figure;
imagesc(deltaE_img);
axis image off;
c = colorbar;
title('ΔE_{2000} map (Lab-based color correction)', 'FontSize', 16, 'FontWeight', 'bold');
clim([0 20]);
c.Label.String = 'ΔE_{2000}';
c.Label.FontSize = 14;
c.Label.FontWeight = 'bold';
set(gca, 'FontSize', 13, 'FontWeight', 'bold');
fprintf('Mean ΔE2000: %.3f\n', mean(deltaE(:)));
fprintf('Median ΔE2000: %.3f\n', median(deltaE(:)));
fprintf('95th percentile ΔE2000: %.3f\n', prctile(deltaE, 95));

%% RGB conversions for display
rgb_before = xyz2rgb(painting_film_xyz_reshaped ./ 100, 'ColorSpace', 'prophoto-rgb');
rgb_before_img = reshape(rgb_before, h, w, 3);

rgb_ref_img = xyz2rgb(painting_ref_xyz_reshaped ./ 100, 'ColorSpace', 'prophoto-rgb');
rgb_ref_img = reshape(rgb_ref_img, h, w, 3);

rgb_corr = xyz2prophoto(corrected_painting_xyz ./ 100, true); % gamma 1.8
rgb_corr_img = reshape(rgb_corr, h, w, 3);

figure;
tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
imshow(rgb_before_img);
title('Before Correction', 'FontWeight', 'bold', 'FontSize', 15);
nexttile;
imshow(rgb_corr_img);
title('After Correction (Lab)', 'FontWeight', 'bold', 'FontSize', 15);
nexttile;
imshow(rgb_ref_img);
title('Reference', 'FontWeight', 'bold', 'FontSize', 15);

%% Save
output_path = '/home/oem/eliza/mac-shared/corrected_poly3_Lab.tif';
% saveProPhotoTIFF(rgb_corr_img, output_path);
fprintf('Saved ProPhoto TIFF: %s\n', output_path);

%% Poly3 features (as before)
function X_poly = poly3_features(input_data)
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);
    X_poly = [ ...
        ones(size(a)), ...
        a, b, c, ...
        a.^2, b.^2, c.^2, ...
        a.*b, a.*c, b.*c, ...
        a.^3, b.^3, c.^3, ...
        a.^2.*b, a.^2.*c, ...
        b.^2.*a, b.^2.*c, ...
        c.^2.*a, c.^2.*b, ...
        a.*b.*c];
end
