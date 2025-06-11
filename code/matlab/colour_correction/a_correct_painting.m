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

%% Colour correction
% ========== Reshape colorchecker data ==========
% xyz_film and xyz_ref: N×3 arrays

min_xyz = min(xyz_film,[],1);     % or just [0 0 0]
max_xyz = max(xyz_film,[],1);     % or just [100 100 100]

painting_film_xyz_clamped = painting_film_xyz; % Copy

for c = 1:3
    painting_film_xyz_clamped(:,:,c) = ...
        min( max(painting_film_xyz(:,:,c), min_xyz(c)), max_xyz(c) );
end
painting_film_xyz_reshaped = reshape(painting_film_xyz_clamped, [], 3);

X_poly = poly3_features(xyz_film);   % N×19
coeffs_xyz = pinv(X_poly) * xyz_ref; % 19×3

% ========== Prepare painting data ==========
[h, w, ~] = size(painting_film_xyz);
% painting_film_xyz_reshaped = reshape(painting_film_xyz, [], 3); % (h*w)×3

% ========== Polynomial expansion and correction ==========
X_poly_painting = poly3_features(painting_film_xyz_reshaped);  % (h*w)×19
corrected_painting_xyz = X_poly_painting * coeffs_xyz;         % (h*w)×3

% ========== Convert both corrected and reference to Lab ==========
painting_ref_xyz_reshaped = reshape(painting_ref_xyz, [], 3);  % (h*w)×3

Lab_corrected = xyz2lab(corrected_painting_xyz);       % (h*w)×3
Lab_reference = xyz2lab(painting_ref_xyz_reshaped);    % (h*w)×3

% ========== ompute ΔE2000 ==========
deltaE = deltaE2000(Lab_reference, Lab_corrected);

% ========== Reshape deltaE to image and display or analyze ==========
deltaE_img = reshape(deltaE, h, w);



%% figure;
imagesc(deltaE_img);
axis image off;  
axis image;
c = colorbar;
title('ΔE_{2000} map (XYZ-based color correction)', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('X', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Y', 'FontSize', 14, 'FontWeight', 'bold');
clim([0 20]);
c.Label.String = 'ΔE_{2000}';
c.Label.FontSize = 14;
c.Label.FontWeight = 'bold';
set(gca, 'FontSize', 13, 'FontWeight', 'bold');
% Optional: Print mean/median error
fprintf('Mean ΔE2000: %.3f\n', mean(deltaE(:)));
fprintf('Median ΔE2000: %.3f\n', median(deltaE(:)));
fprintf('95th percentile ΔE2000: %.3f\n', prctile(deltaE, 95));

%%
% Reshape corrected XYZ to image

rgb_before = xyz2rgb(reshape(painting_film_xyz, [], 3) ./ 100, 'ColorSpace', 'prophoto-rgb');
rgb_before_img = reshape(rgb_before, h, w, 3);
rgb_ref_img = xyz2rgb(reshape(painting_ref_xyz, [], 3) ./ 100, 'ColorSpace', 'prophoto-rgb');
rgb_ref_img = reshape(rgb_ref_img, h, w, 3);

corrected_painting_xyz(~isfinite(corrected_painting_xyz)) = 0;
corrected_painting_xyz = max(corrected_painting_xyz, 0);

rgb_corr = xyz2prophoto(corrected_painting_xyz ./100, true); % with gamma 1.8
rgb_corr_img = reshape(rgb_corr, h, w, 3);            % h × w × 3

figure;
imshow(rgb_corr_img,[]); % This assumes values are [0,1]
title('Corrected Image');

%%
figure;
tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
imshow(rgb_before_img);
title('Before Correction', 'FontWeight', 'bold', 'FontSize', 15);

nexttile;
imshow(rgb_corr_img);
title('After Correction', 'FontWeight', 'bold', 'FontSize', 15);

nexttile;
imshow(rgb_ref_img);
title('Reference', 'FontWeight', 'bold', 'FontSize', 15);




%%
output_path = '/home/oem/eliza/mac-shared/corrected_poly3.tif';
saveProPhotoTIFF(rgb_corr_img, output_path);
fprintf('Saved ProPhoto TIFF: %s\n', output_path);
%%


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
