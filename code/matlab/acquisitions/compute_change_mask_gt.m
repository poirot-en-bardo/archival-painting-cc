clear; close all;

% % === Input paths ===
% before_path = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat';
% after_path  = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';
after_path = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
before_path = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
save_dir   = '/home/oem/eliza/data/xyz_lab_rgb/change_gt';  

% === Load data ===
before = load(before_path);
after  = load(after_path);

% === Reshape Lab images to N x 3 ===
lab_before = reshape(before.Lab_img, [], 3);
lab_after  = reshape(after.Lab_img, [], 3);

% === Compute ΔE2000 ===
deltaE = deltaE2000(lab_before, lab_after);
[h, w, ~] = size(before.Lab_img);
deltaE_img = reshape(deltaE, h, w);

% === Apply combined specular mask ===
mask_combined = before.spec_mask | after.spec_mask;
deltaE_img(mask_combined) = 0;

% === Extract prefix from filename (e.g., "yoda" from "yoda_reflectance_before_xyz.mat") ===
[~, before_filename, ~] = fileparts(before_path);
underscore_idx = strfind(before_filename, '_');
prefix = before_filename(1 : underscore_idx(1) - 1);  % Extract prefix before first underscore

% === Construct variable name and file path ===
var_name = [prefix '_deltaE_img'];
save_path = fullfile(save_dir, [var_name '.mat']);

% === Save with dynamic variable name ===
eval([var_name ' = deltaE_img;']);
save(save_path, var_name);

fprintf('ΔE matrix saved to: %s (variable name: %s)\n', save_path, var_name);
%%
% === Reload and display ===
loaded = load(save_path);
deltaE_loaded = loaded.(var_name);  % Access dynamic field

% === Display to confirm ===
figure;
imagesc(deltaE_loaded); axis image off;
title(['\DeltaE_{00} Reloaded: ' var_name], 'FontSize', 16, 'FontWeight', 'bold');
colorbar;
clim([0 20]);

%%
% === Create binary change mask where ΔE >= 6 ===
% === Create binary change mask where ΔE >= threshold ===
threshold = 6;
change_mask = deltaE_img >= threshold;

% === Save the change mask ===
mask_var_name = [prefix '_change_mask'];
mask_save_path = fullfile(save_dir, [mask_var_name '.mat']);
% eval([mask_var_name ' = change_mask;']);
save(mask_save_path, 'change_mask');

fprintf('Change mask saved to: %s (variable name: %s)\n', mask_save_path, mask_var_name);


% === Display the change mask with better layout ===
fig_mask = figure('Units', 'normalized', 'Position', [0.2, 0.3, 0.6, 0.6]);
tl = tiledlayout(1,1, 'Padding','compact', 'TileSpacing','compact');
ax = nexttile;
imshow(change_mask, 'Parent', ax);
title(ax, ['Change Mask (\DeltaE_{00} ≥ ' num2str(threshold) '): ' 'Yoda'], ...
    'FontSize', 18, 'FontWeight', 'bold');

% === Save the mask image with enough space for title ===
mask_img_path = fullfile(save_dir, [prefix '_change_mask.png']);
exportgraphics(fig_mask, mask_img_path, ...
    'Resolution', 300, 'BackgroundColor', 'white', 'ContentType', 'image');

fprintf('Change mask image saved to: %s\n', mask_img_path);
