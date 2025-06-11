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

%%
% --------- Parameters ---------
lambda = 0;  % Ridge regularization parameter (try 0.1, 1, or 10)

% --------- Feature expansion for training (quadratic degree 2) ---------
X_poly = poly2_features(xyz_film);   % [N×10]
Y = xyz_ref;                         % [N×3]

nF = size(X_poly,2);                 % number of features, here 10

% --------- Ridge regression training ---------
coeffs_ridge = (X_poly' * X_poly + lambda * eye(nF)) \ (X_poly' * Y); % [10×3]

% --------- Prepare painting data ---------
painting_film_xyz_reshaped = reshape(painting_film_xyz, [], 3);  % [(h*w)×3]

% --------- Feature expansion for painting ---------
X_poly_painting = poly2_features(painting_film_xyz_reshaped);    % [(h*w)×10]

% --------- Apply correction ---------
corrected_painting_xyz = X_poly_painting * coeffs_ridge;         % [(h*w)×3]

% --------- Clamp negative/NaN results ---------
corrected_painting_xyz(~isfinite(corrected_painting_xyz)) = 0;
corrected_painting_xyz = max(corrected_painting_xyz, 0);

% --------- Reshape to image ---------
corrected_painting_xyz_img = reshape(corrected_painting_xyz, size(painting_film_xyz));



%%
% --------- Convert to Lab for ΔE analysis (if needed) ---------
painting_ref_xyz_reshaped = reshape(painting_ref_xyz, [], 3);
Lab_corrected = xyz2lab(corrected_painting_xyz);           % [(h*w)×3]
Lab_reference = xyz2lab(painting_ref_xyz_reshaped);        % [(h*w)×3]
deltaE = deltaE2000(Lab_reference, Lab_corrected);
deltaE_img = reshape(deltaE, size(painting_film_xyz,1), size(painting_film_xyz,2));
figure; imagesc(deltaE_img); axis image; colorbar; clim([0 10]);
title('\DeltaE_{2000} (XYZ, ridge poly2)');

%%
% Optional: Print mean/median error
fprintf('Mean ΔE2000: %.3f\n', mean(deltaE(:)));
fprintf('Median ΔE2000: %.3f\n', median(deltaE(:)));
fprintf('95th percentile ΔE2000: %.3f\n', prctile(deltaE, 95));

%%
% --------- Convert to ProPhoto RGB for visualization ---------
rgb_corr = xyz2prophoto(corrected_painting_xyz ./100, true); % gamma 1.8
rgb_corr_img = reshape(rgb_corr, size(painting_film_xyz,1), size(painting_film_xyz,2), 3);
figure; imshow(rgb_corr_img); title('Corrected Painting (ProPhoto RGB, ridge poly2)');

% --------- Save as 16-bit TIFF ---------
output_path = 'corrected_prophoto_ridge_poly2.tif';
% saveProPhotoTIFF(rgb_corr_img, output_path);

fprintf('Saved ProPhoto TIFF: %s\n', output_path);

% --------- Quadratic polynomial feature function ---------
function X_poly = poly2_features(input_data)
    a = input_data(:,1);
    b = input_data(:,2);
    c = input_data(:,3);

    X_poly = [ ...
        ones(size(a)), ...
        a, b, c, ...
        a.^2, b.^2, c.^2, ...
        a.*b, a.*c, b.*c ];
end
