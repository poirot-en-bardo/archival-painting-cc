clear; close all;

%% === File paths ===
% painting_before_binned_file = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_before_xyz.mat';
% painting_after_binned_file  = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/cactus_reflectance_after_reg_xyz.mat';
% painting_after_ref_file     = '/home/oem/eliza/data/reflectance/registered/cactus_reflectance_after_reg.mat';
% film_file                   = '/home/oem/eliza/data/xyz_lab_rgb/film/cactus_halogen_kodak_exp0.mat';
% change_mask_file            = '/home/oem/eliza/data/xyz_lab_rgb/change_gt/cactus_change_mask.mat';
painting_before_binned_file = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_before_xyz.mat';
painting_after_binned_file  = '/home/oem/eliza/data/xyz_lab_rgb/hyspex/yoda_reflectance_after_reg_xyz.mat';
painting_after_ref_file     = '/home/oem/eliza/data/reflectance/registered/yoda_reflectance_after_reg.mat';
film_file                   = '/home/oem/eliza/data/xyz_lab_rgb/film/yoda_halogen_fuji_exp0.mat';
change_mask_file            = '/home/oem/eliza/data/xyz_lab_rgb/change_gt/yoda_change_mask.mat';


rng(42);

% === Save clustered data ===
[~, film_base_name, ~] = fileparts(film_file);
underscore_idx = strfind(film_base_name, '_');
prefix = film_base_name(1:underscore_idx(1)-1);

save_root = '/home/oem/eliza/data/xyz_lab_rgb/clustered/filtered';
save_dir = fullfile(save_root, film_base_name);
if ~exist(save_dir, 'dir'); mkdir(save_dir); end


%% === Load data ===
painting_after_ref = load(painting_after_ref_file);
painting_after_binned = load(painting_after_binned_file);
painting_before_binned = load(painting_before_binned_file);
film_data = load(film_file);
change_mask = load(change_mask_file).change_mask;

spec_mask_combined = painting_before_binned.spec_mask | painting_after_binned.spec_mask;
data_cube = painting_after_ref.data_cube;

%% === Bin the spectral cube ===
binSize = 2;
[H, W, B] = size(data_cube);
h = floor(H / binSize); w = floor(W / binSize);
data_cube_gpu = gpuArray(data_cube);
binned_cube_gpu = zeros(h, w, B, 'like', data_cube_gpu);
for b = 1:B
    band = data_cube_gpu(:,:,b);
    band = band(1:h*binSize, 1:w*binSize);
    band = reshape(band, binSize, h, binSize, w);
    band = permute(band, [2 4 1 3]);
    binned_cube_gpu(:,:,b) = mean(mean(band, 3), 4);
end
binned_cube_gpu = flip(binned_cube_gpu, 2);

%% === Select valid pixels ===
spec_mask_combined = flip(spec_mask_combined, 2);
change_mask = flip(change_mask, 2);
valid_mask = ~(spec_mask_combined | change_mask);

pixels = reshape(gather(binned_cube_gpu), [], B);
valid_pixels = pixels(valid_mask(:), :);

%% === K-means clustering ===
rng(42); nColors = 80;
[cluster_idx, ~] = kmeans(valid_pixels, nColors, 'MaxIter', 1000, 'Replicates', 3, 'Start', 'sample');
%%
cluster_idx_save = cluster_idx;
%%
cluster_idx_file = fullfile(save_dir, [prefix '_cluster_idx_raw.mat']);
save(cluster_idx_file, 'cluster_idx');


%% ===  CLUSTER FILTERING, SORTING, AND DISPLAY ===
% Reshape color space images to [num_pixels x 3]
Lab_after  = reshape(painting_after_binned.Lab_img, [], 3);
XYZ_after  = reshape(painting_after_binned.XYZ_img, [], 3);
RGB_after  = reshape(painting_after_binned.RGB_img, [], 3);

Lab_before = reshape(painting_before_binned.Lab_img, [], 3);
XYZ_before = reshape(painting_before_binned.XYZ_img, [], 3);

Lab_film   = reshape(film_data.Lab_img, [], 3);
XYZ_film   = reshape(film_data.XYZ_img, [], 3);
RGB_film   = reshape(film_data.RGB_img, [], 3);

% Apply valid_mask
lab_after_valid   = Lab_after(valid_mask(:), :);
xyz_after_valid   = XYZ_after(valid_mask(:), :);
rgb_after_valid   = RGB_after(valid_mask(:), :);

lab_before_valid  = Lab_before(valid_mask(:), :);
xyz_before_valid  = XYZ_before(valid_mask(:), :);

lab_film_valid    = Lab_film(valid_mask(:), :);
xyz_film_valid    = XYZ_film(valid_mask(:), :);
rgb_film_valid    = RGB_film(valid_mask(:), :);

%%
initial_nColors = nColors;
original_cluster_idx = cluster_idx;  

%% cluster filtering
% Use original, unfiltered data
% mean_thresh = 6.0;
% max_thresh  = 20.0;
% 
% cluster_intraDE_mean = zeros(initial_nColors, 1);
% cluster_intraDE_max  = zeros(initial_nColors, 1);
% lab_after_clustered_tmp = zeros(initial_nColors, 3);
% 
% for k = 1:initial_nColors
%     idx_k = (original_cluster_idx == k);
%     if nnz(idx_k) == 0
%         cluster_intraDE_mean(k) = NaN;
%         cluster_intraDE_max(k)  = NaN;
%         continue;
%     end
%     pixel_lab = lab_after_valid(idx_k, :);
%     mean_lab = mean(pixel_lab, 1, 'omitnan');
%     lab_after_clustered_tmp(k,:) = mean_lab;
%     DE = deltaE2000(repmat(mean_lab, size(pixel_lab,1), 1), pixel_lab);
%     cluster_intraDE_mean(k) = mean(DE, 'omitnan');
%     cluster_intraDE_max(k)  = max(DE, [], 'omitnan');
% end
% 
% valid_clusters = find( ...
%     cluster_intraDE_mean < mean_thresh & ...
%     cluster_intraDE_max  < max_thresh);
% 
% fprintf('Kept %d of %d clusters (ΔE₀₀ mean < %.1f, max < %.1f)\n', ...
%     numel(valid_clusters), initial_nColors, mean_thresh, max_thresh);
% 
% % === Update cluster_idx ===
% remap_valid = zeros(1, initial_nColors);
% remap_valid(valid_clusters) = 1:length(valid_clusters);
% 
% cluster_idx_filtered = nan(size(original_cluster_idx));
% for k = 1:numel(valid_clusters)
%     original_id = valid_clusters(k);
%     cluster_idx_filtered(original_cluster_idx == original_id) = k;
% end
% 
% cluster_idx = cluster_idx_filtered;
% nColors = numel(valid_clusters);  % update only for downstream filtered steps

% === Recompute filtered cluster means ===
lab_after_clustered = zeros(nColors, 3);
xyz_after_clustered = zeros(nColors, 3);
lab_before_clustered = zeros(nColors, 3);

for k = 1:nColors
    mask_k = cluster_idx == k;
    lab_after_clustered(k,:) = mean(lab_after_valid(mask_k,:), 1, 'omitnan');
    xyz_after_clustered(k,:) = mean(xyz_after_valid(mask_k,:), 1, 'omitnan');
    lab_before_clustered(k,:) = mean(lab_before_valid(mask_k,:), 1, 'omitnan');
end


%%
% === Sort Lab clusters by chromaticity ===
Lab = lab_after_clustered;
chroma = sqrt(Lab(:,2).^2 + Lab(:,3).^2);
achromaticThresh = 7;
achromaticIdx = find(chroma < achromaticThresh);
chromaticIdx = find(chroma >= achromaticThresh);


nGroups = 3;
rng(123, 'twister');  % Or whatever seed/method you prefer
[grp, ~] = kmeans(Lab(chromaticIdx,2:3), nGroups, 'Replicates', 20, 'Start', 'plus');

finalOrder = [];
[~, achroOrd] = sort(Lab(achromaticIdx,1), 'descend');
finalOrder = [finalOrder; achromaticIdx(achroOrd)];
for g = 1:nGroups
    idx = chromaticIdx(grp == g);
    [~, ord] = sort(Lab(idx,1), 'descend');
    finalOrder = [finalOrder; idx(ord)];
end

% === Remap cluster_idx using sorted order ===
remap_final = zeros(1, nColors);
for i = 1:nColors
    remap_final(finalOrder(i)) = i;
end

valid_idx = ~isnan(cluster_idx);
cluster_idx(valid_idx) = remap_final(cluster_idx(valid_idx));

cluster_idx_sorted_file = fullfile(save_dir, [prefix '_cluster_idx_sorted.mat']);
save(cluster_idx_sorted_file, 'cluster_idx');

% === Reorder cluster means ===
% xyz_after_clustered = xyz_after_clustered(finalOrder, :);
% lab_after_clustered = lab_after_clustered(finalOrder, :);
% === Recompute cluster means after remapping ===
lab_after_clustered = zeros(nColors, 3);
xyz_after_clustered = zeros(nColors, 3);
lab_before_clustered = zeros(nColors, 3);

for k = 1:nColors
    mask_k = cluster_idx == k;
    lab_after_clustered(k,:) = mean(lab_after_valid(mask_k,:), 1, 'omitnan');
    xyz_after_clustered(k,:) = mean(xyz_after_valid(mask_k,:), 1, 'omitnan');
    lab_before_clustered(k,:) = mean(lab_before_valid(mask_k,:), 1, 'omitnan');
end


% === Compute RGB palette ===
RGB = max(0, min(1, xyz2rgb(xyz_after_clustered ./ 100, 'ColorSpace','srgb', 'WhitePoint','d50')));

% === Generate and save RGB palette ===
patchSize = 60;
grid_w = 10; grid_h = ceil(nColors / grid_w);
palette_img = ones(grid_h * patchSize, grid_w * patchSize, 3);

for k = 1:nColors
    row = floor((k-1)/grid_w);
    col = mod((k-1), grid_w);
    r_idx = (row*patchSize+1):((row+1)*patchSize);
    c_idx = (col*patchSize+1):((col+1)*patchSize);
    palette_img(r_idx, c_idx, :) = repmat(reshape(RGB(k,:),1,1,3), patchSize, patchSize, 1);
end
%
figure('Position', [100 100 800 600]);
imshow(palette_img);
exportgraphics(gcf, fullfile(save_dir, [prefix '_cluster_palette_after_rgb.png']), 'Resolution', 300);
%%
cluster_idx_unfiltered = cluster_idx;

%% === Recompute RGB / XYZ / Lab cluster means after final sorted cluster_idx ===
lab_after_clustered = zeros(nColors, 3);
xyz_after_clustered = zeros(nColors, 3);
lab_before_clustered = zeros(nColors, 3);
rgb_film_clustered = zeros(nColors, 3);
rgb_after_clustered = zeros(nColors, 3);
xyz_film_clustered = zeros(nColors, 3);
xyz_before_clustered = zeros(nColors, 3);
lab_film_clustered = zeros(nColors, 3);

for k = 1:nColors
    mask_k = cluster_idx == k;
    lab_after_clustered(k,:) = mean(lab_after_valid(mask_k,:), 1, 'omitnan');
    xyz_after_clustered(k,:) = mean(xyz_after_valid(mask_k,:), 1, 'omitnan');
    lab_before_clustered(k,:) = mean(lab_before_valid(mask_k,:), 1, 'omitnan');
    rgb_film_clustered(k,:) = mean(rgb_film_valid(mask_k,:), 1, 'omitnan');
    rgb_after_clustered(k,:) = mean(rgb_after_valid(mask_k,:), 1, 'omitnan');
    xyz_film_clustered(k,:) = mean(xyz_film_valid(mask_k,:), 1, 'omitnan');
    xyz_before_clustered(k,:) = mean(xyz_before_valid(mask_k,:), 1, 'omitnan');
    lab_film_clustered(k,:) = mean(lab_film_valid(mask_k,:), 1, 'omitnan');
end


srgb_after_clustered = xyz2rgb(xyz_after_clustered ./100, 'WhitePoint','d50');
srgb_before_clustered = xyz2rgb(xyz_before_clustered ./100, 'WhitePoint','d50');

% Palette before vs after
img_diag_before = ones(grid_h * patchSize, grid_w * patchSize, 3);

for k = 1:nColors
    row = floor((k-1) / grid_w);
    col = mod((k-1), grid_w);
    r_idx = (row * patchSize + 1):((row + 1) * patchSize);
    c_idx = (col * patchSize + 1):((col + 1) * patchSize);

    % Prepare film and after patches
    patch_film  = repmat(reshape(srgb_before_clustered(k, :), 1, 1, 3), patchSize, patchSize, 1);
    patch_after = repmat(reshape(srgb_after_clustered(k, :), 1, 1, 3), patchSize, patchSize, 1);

    % Diagonal mask
    [X, Y] = meshgrid(1:patchSize, 1:patchSize);
    diag_mask = Y > X;  % lower triangle gets after

    patch = patch_film; % Start from film values
    for c = 1:3
        temp = patch(:,:,c);
        temp2 = patch_after(:,:,c);
        temp(diag_mask) = temp2(diag_mask); % Correct assignment
        patch(:,:,c) = temp;
    end

    img_diag_before(r_idx, c_idx, :) = patch;
end

% Display and save
figure('Position', [100 100 1200 1100]);
tiledlayout(1,1, 'Padding','compact', 'TileSpacing','compact');

nexttile;
imshow(img_diag_before);
title('Cluster Before-After: Before (top right) vs After (bottom left)', ...
    'FontWeight', 'bold', 'FontSize', 22);

exportgraphics(gcf, fullfile(save_dir, [prefix '_before_after_cluster_rgb_all.png']), 'Resolution', 300);

%%
DE_film_after   = deltaE2000(lab_film_clustered, lab_after_clustered);
DE_before_after = deltaE2000(lab_before_clustered, lab_after_clustered);
DE_film_before = deltaE2000(lab_film_clustered, lab_before_clustered);


% Convert to grid format
dE_grid_film_after   = zeros(grid_h, grid_w);
dE_grid_before_after = zeros(grid_h, grid_w);
dE_grid_film_before = zeros(grid_h, grid_w);
for k = 1:nColors
    row = floor((k-1) / grid_w) + 1;
    col = mod((k-1), grid_w) + 1;
    dE_grid_film_after(row, col)   = DE_film_after(k);
    dE_grid_before_after(row, col) = DE_before_after(k);
    dE_grid_film_before(row, col) = DE_film_before(k);

end

figure('Position', [100 100 800 700]);
imagesc(dE_grid_before_after);
axis image off;
colormap(jet(256));
cb = colorbar;
cb.FontSize = 16;
cb.FontWeight = 'bold';
clim([0 10]);
title('\DeltaE_{00} Before vs After', 'FontSize', 18, 'FontWeight', 'bold');

% Overlay text
for k = 1:nColors
    row = floor((k-1)/grid_w) + 1;
    col = mod((k-1), grid_w) + 1;
    text(col, row, sprintf('%.1f', DE_before_after(k)), ...
        'HorizontalAlignment','center', 'Color','w', ...
        'FontSize', 12, 'FontWeight', 'bold');
end

exportgraphics(gcf, fullfile(save_dir, [prefix '_palette_deltaE_before_vs_after_all.png']), 'Resolution', 300);

%%
% === Filter clusters by ΔE₀₀ Before vs After threshold ===
deltaE_thresh = 6;
valid_clusters = find(DE_before_after <= deltaE_thresh);
nColors = numel(valid_clusters);

fprintf('Keeping %d of %d clusters (ΔE₀₀ ≤ %.1f)\n', ...
    nColors, length(DE_before_after), deltaE_thresh);

% Remap cluster_idx to include only valid clusters
remap_filtered = zeros(1, length(DE_before_after));
remap_filtered(valid_clusters) = 1:nColors;

cluster_idx_filtered = nan(size(cluster_idx));
for i = 1:nColors
    old_k = valid_clusters(i);
    cluster_idx_filtered(cluster_idx == old_k) = i;
end
cluster_idx = cluster_idx_filtered;  % Overwrite
grid_w = 10; grid_h = ceil(nColors / grid_w);

%% === Recompute RGB / XYZ / Lab cluster means after final sorted cluster_idx ===
lab_after_clustered = zeros(nColors, 3);
xyz_after_clustered = zeros(nColors, 3);
lab_before_clustered = zeros(nColors, 3);
rgb_film_clustered = zeros(nColors, 3);
rgb_after_clustered = zeros(nColors, 3);
xyz_film_clustered = zeros(nColors, 3);
xyz_before_clustered = zeros(nColors, 3);
lab_film_clustered = zeros(nColors, 3);

for k = 1:nColors
    mask_k = cluster_idx == k;
    lab_after_clustered(k,:) = mean(lab_after_valid(mask_k,:), 1, 'omitnan');
    xyz_after_clustered(k,:) = mean(xyz_after_valid(mask_k,:), 1, 'omitnan');
    lab_before_clustered(k,:) = mean(lab_before_valid(mask_k,:), 1, 'omitnan');
    rgb_film_clustered(k,:) = mean(rgb_film_valid(mask_k,:), 1, 'omitnan');
    rgb_after_clustered(k,:) = mean(rgb_after_valid(mask_k,:), 1, 'omitnan');
    xyz_film_clustered(k,:) = mean(xyz_film_valid(mask_k,:), 1, 'omitnan');
    xyz_before_clustered(k,:) = mean(xyz_before_valid(mask_k,:), 1, 'omitnan');
    lab_film_clustered(k,:) = mean(lab_film_valid(mask_k,:), 1, 'omitnan');
end

% === Save data ===
save(fullfile(save_dir, [film_base_name '_cluster_rgb.mat']), ...
    'rgb_film_clustered', 'xyz_film_clustered', 'lab_film_clustered');

save(fullfile(save_dir, [prefix '_painting_after_cluster_rgb.mat']), ...
    'rgb_after_clustered', 'xyz_after_clustered', 'xyz_before_clustered', ...
    'lab_after_clustered', 'lab_before_clustered');

save(fullfile(save_dir, [prefix '_spectral_cluster_metadata.mat']), ...
    'cluster_idx', 'valid_mask', 'nColors', '-v7.3');

%% === Visualize comparison === Compute sRGB versions ===

srgb_film_clustered = xyz2rgb(xyz_film_clustered ./100, 'WhitePoint','d50');
srgb_after_clustered = xyz2rgb(xyz_after_clustered ./100, 'WhitePoint','d50');
srgb_before_clustered = xyz2rgb(xyz_before_clustered ./100, 'WhitePoint','d50');


%% Palette film vs after
% === Create diagonal split grid image (Film top-left, After bottom-right) ===
img_diag = ones(grid_h * patchSize, grid_w * patchSize, 3);

for k = 1:nColors
    row = floor((k-1) / grid_w);
    col = mod((k-1), grid_w);
    r_idx = (row * patchSize + 1):((row + 1) * patchSize);
    c_idx = (col * patchSize + 1):((col + 1) * patchSize);

    % Prepare film and after patches
    patch_film  = repmat(reshape(srgb_film_clustered(k, :), 1, 1, 3), patchSize, patchSize, 1);
    patch_after = repmat(reshape(srgb_after_clustered(k, :), 1, 1, 3), patchSize, patchSize, 1);

    % Diagonal mask
    [X, Y] = meshgrid(1:patchSize, 1:patchSize);
    diag_mask = Y > X;  % lower triangle gets after

    patch = patch_film; % Start from film values
    for c = 1:3
        temp = patch(:,:,c);
        temp2 = patch_after(:,:,c);
        temp(diag_mask) = temp2(diag_mask); % Correct assignment
        patch(:,:,c) = temp;
    end

    img_diag(r_idx, c_idx, :) = patch;
end


% Display and save with proper title spacing
figure('Position', [100 100 1200 1100]);
tiledlayout(1,1, 'Padding','compact', 'TileSpacing','compact');

nexttile;
imshow(img_diag);
title('Palette: After (bottom left) vs Film (top right)', ...
    'FontWeight', 'bold', 'FontSize', 22);

exportgraphics(gcf, fullfile(save_dir, [film_base_name 'film_after_cluster_rgb.png']), 'Resolution', 300);


%% Palette film vs before

img_diag = ones(grid_h * patchSize, grid_w * patchSize, 3);

for k = 1:nColors
    row = floor((k-1) / grid_w);
    col = mod((k-1), grid_w);
    r_idx = (row * patchSize + 1):((row + 1) * patchSize);
    c_idx = (col * patchSize + 1):((col + 1) * patchSize);

    % Prepare film and after patches
    patch_film  = repmat(reshape(srgb_film_clustered(k, :), 1, 1, 3), patchSize, patchSize, 1);
    patch_after = repmat(reshape(srgb_before_clustered(k, :), 1, 1, 3), patchSize, patchSize, 1);

    % Diagonal mask
    [X, Y] = meshgrid(1:patchSize, 1:patchSize);
    diag_mask = Y > X;  % lower triangle gets after

    patch = patch_film; % Start from film values
    for c = 1:3
        temp = patch(:,:,c);
        temp2 = patch_after(:,:,c);
        temp(diag_mask) = temp2(diag_mask); % Correct assignment
        patch(:,:,c) = temp;
    end

    img_diag(r_idx, c_idx, :) = patch;
end


% Display and save with proper title spacing
figure('Position', [100 100 1200 1100]);
tiledlayout(1,1, 'Padding','compact', 'TileSpacing','compact');

nexttile;
imshow(img_diag);
title('Palette: Before (bottom left) vs Film (top right)', ...
    'FontWeight', 'bold', 'FontSize', 22);

exportgraphics(gcf, fullfile(save_dir, [film_base_name 'film_before_cluster_rgb.png']), 'Resolution', 300);

%% Palette before vs after
img_diag_before = ones(grid_h * patchSize, grid_w * patchSize, 3);

for k = 1:nColors
    row = floor((k-1) / grid_w);
    col = mod((k-1), grid_w);
    r_idx = (row * patchSize + 1):((row + 1) * patchSize);
    c_idx = (col * patchSize + 1):((col + 1) * patchSize);

    % Prepare film and after patches
    patch_film  = repmat(reshape(srgb_before_clustered(k, :), 1, 1, 3), patchSize, patchSize, 1);
    patch_after = repmat(reshape(srgb_after_clustered(k, :), 1, 1, 3), patchSize, patchSize, 1);

    % Diagonal mask
    [X, Y] = meshgrid(1:patchSize, 1:patchSize);
    diag_mask = Y > X;  % lower triangle gets after

    patch = patch_film; % Start from film values
    for c = 1:3
        temp = patch(:,:,c);
        temp2 = patch_after(:,:,c);
        temp(diag_mask) = temp2(diag_mask); % Correct assignment
        patch(:,:,c) = temp;
    end

    img_diag_before(r_idx, c_idx, :) = patch;
end

% Display and save
figure('Position', [100 100 1200 1100]);
tiledlayout(1,1, 'Padding','compact', 'TileSpacing','compact');

nexttile;
imshow(img_diag_before);
title('Cluster Before-After Stable: Before (top right) vs After (bottom left)', ...
    'FontWeight', 'bold', 'FontSize', 22);

exportgraphics(gcf, fullfile(save_dir, [prefix '_before_after_cluster_rgb.png']), 'Resolution', 300);


%% === Compute Delta E 2000 values ===
DE_film_after   = deltaE2000(lab_film_clustered, lab_after_clustered);
DE_before_after = deltaE2000(lab_before_clustered, lab_after_clustered);
DE_film_before = deltaE2000(lab_film_clustered, lab_before_clustered);


% Convert to grid format
dE_grid_film_after   = zeros(grid_h, grid_w);
dE_grid_before_after = zeros(grid_h, grid_w);
dE_grid_film_before = zeros(grid_h, grid_w);
for k = 1:nColors
    row = floor((k-1) / grid_w) + 1;
    col = mod((k-1), grid_w) + 1;
    dE_grid_film_after(row, col)   = DE_film_after(k);
    dE_grid_before_after(row, col) = DE_before_after(k);
    dE_grid_film_before(row, col) = DE_film_before(k);

end

%% === Display: Film vs After ===
figure('Position', [100 100 800 700]);
imagesc(dE_grid_film_after);
axis image off;
colormap(jet(256));
cb = colorbar;
cb.FontSize = 16;
cb.FontWeight = 'bold';clim([0 10]);  % adjust if needed
title('\DeltaE_{00} Film vs After', 'FontSize', 18, 'FontWeight', 'bold');

% Overlay text
for k = 1:nColors
    row = floor((k-1)/grid_w) + 1;
    col = mod((k-1), grid_w) + 1;
    text(col, row, sprintf('%.1f', DE_film_after(k)), ...
        'HorizontalAlignment','center', 'Color','w', ...
        'FontSize', 12, 'FontWeight', 'bold');
end

exportgraphics(gcf, fullfile(save_dir, [film_base_name 'palette_deltaE_film_vs_after.png']), 'Resolution', 300);
%% delta e film vs before

% Grid format
dE_grid_film_before = zeros(grid_h, grid_w);
for k = 1:nColors
    row = floor((k-1) / grid_w) + 1;
    col = mod((k-1), grid_w) + 1;
    dE_grid_film_before(row, col) = DE_film_before(k);
end

% Display
figure('Position', [100 100 800 700]);
imagesc(dE_grid_film_before);
axis image off;
colormap(jet(256));
cb = colorbar;
cb.FontSize = 16;
cb.FontWeight = 'bold';clim([0 10]);
title('\DeltaE_{00} Film vs Before', 'FontSize', 18, 'FontWeight', 'bold');

% Overlay text
for k = 1:nColors
    row = floor((k-1)/grid_w) + 1;
    col = mod((k-1), grid_w) + 1;
    text(col, row, sprintf('%.1f', DE_film_before(k)), ...
        'HorizontalAlignment','center', 'Color','w', ...
        'FontSize', 12, 'FontWeight', 'bold');
end

exportgraphics(gcf, ...
    fullfile(save_dir, [film_base_name '_palette_deltaE_film_vs_before.png']), ...
    'Resolution', 300);

%% === Display: Before vs After ===
figure('Position', [100 100 800 700]);
imagesc(dE_grid_before_after);
axis image off;
colormap(jet(256));

cb = colorbar;
cb.FontSize = 16;
cb.FontWeight = 'bold';
clim([0 10]);
title('\DeltaE_{00} Before vs After Filtered', 'FontSize', 18, 'FontWeight', 'bold');

% Overlay text
for k = 1:nColors
    row = floor((k-1)/grid_w) + 1;
    col = mod((k-1), grid_w) + 1;
    text(col, row, sprintf('%.1f', DE_before_after(k)), ...
        'HorizontalAlignment','center', 'Color','w', ...
        'FontSize', 12, 'FontWeight', 'bold');
end

exportgraphics(gcf, fullfile(save_dir, [prefix '_palette_deltaE_before_vs_after.png']), 'Resolution', 300);

%% Process new film photo

%load...
% 
% Lab_film_new = reshape(new_film_data.Lab_img, [], 3);  % or RGB/XYZ
% film_valid = Lab_film_new(valid_mask, :);  % Match only valid pixels
% 
% rgb_film_new_clustered = zeros(nColors, 3);
% for k = 1:nColors
%     mask_k = cluster_idx == k;
%     rgb_film_new_clustered(k,:) = mean(film_valid(mask_k,:), 1, 'omitnan');
% end
% 
