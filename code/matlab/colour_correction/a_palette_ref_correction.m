clear; close all;

data_path = '/home/oem/eliza/mac-shared/registered_output/yoda_halogen_fuji_exp0_all3/yoda_halogen_fuji_exp0_all3.mat';
data = load(data_path);

after_rgb   = data.after_rgb;
after_xyz   = data.after_xyz;
before_xyz  = data.before_xyz;
before_rgb  = data.before_rgb;
film_xyz    = data.film_xyz;
film_rgb    = data.film_rgb;

% crop image
figure; imshow(after_rgb);
title('Draw a rectangle on AFTER to crop (double-click when done)');
h = drawrectangle('InteractionsAllowed','all');
wait(h);
crop_rect = round(h.Position);

x1 = crop_rect(1);
y1 = crop_rect(2);
w  = crop_rect(3);
h_ = crop_rect(4);

x1 = max(1, x1);   y1 = max(1, y1);
x2 = min(size(after_rgb,2), x1 + w - 1);
y2 = min(size(after_rgb,1), y1 + h_ - 1);

after_rgb_cropped   = after_rgb(y1:y2, x1:x2, :);
before_rgb_cropped  = before_rgb(y1:y2, x1:x2, :);
film_rgb_cropped    = film_rgb(y1:y2, x1:x2, :);

after_xyz_cropped   = after_xyz(y1:y2, x1:x2, :);
before_xyz_cropped  = before_xyz(y1:y2, x1:x2, :);
film_xyz_cropped    = film_xyz(y1:y2, x1:x2, :);

%%
after_xyz_vec = reshape(after_xyz_cropped, [], 3);
after_lab_cropped = xyz2lab(after_xyz_vec);
after_lab_cropped = reshape(after_lab_cropped, size(after_xyz_cropped));

% Before
before_xyz_vec = reshape(before_xyz_cropped, [], 3);
before_lab_cropped = xyz2lab(before_xyz_vec);
before_lab_cropped = reshape(before_lab_cropped, size(before_xyz_cropped));

% Film
film_xyz_vec = reshape(film_xyz_cropped, [], 3);
film_lab_cropped = xyz2lab(film_xyz_vec);
film_lab_cropped = reshape(film_lab_cropped, size(film_xyz_cropped));
%%

thresh = 95;
mask_after  = all(after_xyz_cropped  < thresh, 3);
mask_before = all(before_xyz_cropped < thresh, 3);
mask_film   = all(film_xyz_cropped   < thresh, 3);
mask_all    = mask_after & mask_before & mask_film;

dE_map_ref = deltaE(before_lab_cropped, after_lab_cropped);
dE_masked = dE_map_ref;
dE_masked(~mask_all) = NaN;
%%
figure;
imagesc(dE_masked);
colormap jet; 
% clim([nanmin(dE_masked(:)), 10]);
colorbar;
title('\DeltaE_{00}');
axis image off;

dE_valid = dE_masked(mask_all);
fprintf('Mean ΔE00: %.3f\n', mean(dE_valid));
fprintf('Median ΔE00: %.3f\n', median(dE_valid));
fprintf('95th percentile ΔE00: %.3f\n', prctile(dE_valid, 95));

% Overlay
figure; imshowpair(before_rgb_cropped, after_rgb_cropped);
title('Before vs After Overlay');
