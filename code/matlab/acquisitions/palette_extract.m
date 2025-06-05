% --- Step 1: Load cube
cubeFile1 = '/home/oem/eliza/data/processed/reflectance/before/cactus_halogen_reflectance_full.hdr';
cubeFile2 = '/home/oem/eliza/data/processed/reflectance/after/cactus_halogen_reflectance_after_full.hdr';

hcube1 = hypercube(cubeFile1); cube1 = hcube1.DataCube; wl1 = hcube1.Wavelength;
hcube2 = hypercube(cubeFile2); cube2 = hcube2.DataCube; wl2 = hcube2.Wavelength;
valid_idx1 = find(wl1 >= 380 & wl1 <= 780);
cube1 = cube1(:,:,valid_idx1); wl1 = wl1(valid_idx1);
valid_idx2 = find(wl2 >= 380 & wl2 <= 780);
cube2 = cube2(:,:,valid_idx2); wl2 = wl2(valid_idx2);

% Crop 
figure; imagesc(cube1(:,:,10)); axis image; colormap gray; title('Crop FILM');
roi1 = drawrectangle; r1 = round(roi1.Position);
x1a = max(1, r1(1)); y1a = max(1, r1(2));
x1b = min(size(cube1,2), x1a + r1(3) - 1); y1b = min(size(cube1,1), y1a + r1(4) - 1);
cube1_crop = cube1(y1a:y1b, x1a:x1b, :);

figure; imagesc(cube2(:,:,10)); axis image; colormap gray; title('Crop PAINTING');
roi2 = drawrectangle; r2 = round(roi2.Position);
x2a = max(1, r2(1)); y2a = max(1, r2(2));
x2b = min(size(cube2,2), x2a + r2(3) - 1); y2b = min(size(cube2,1), y2a + r2(4) - 1);
cube2_crop = cube2(y2a:y2b, x2a:x2b, :);

%% Step 2: Subsample each
tmp = imresize(cube1_crop(:,:,1), 0.15, 'bilinear');
[new_h1, new_w1] = size(tmp);
cube1_sub = zeros(new_h1, new_w1, size(cube1_crop,3));
for b = 1:size(cube1_crop,3)
    cube1_sub(:,:,b) = imresize(cube1_crop(:,:,b), 0.15, 'bilinear');
end

tmp = imresize(cube2_crop(:,:,1), 0.15, 'bilinear');
[new_h1, new_w1] = size(tmp);
cube2_sub = zeros(new_h1, new_w1, size(cube2_crop,3));
for b = 1:size(cube2_crop,3)
    cube2_sub(:,:,b) = imresize(cube2_crop(:,:,b), 0.15, 'bilinear');
end



%% Step 3: Mask specular/clipped for each (threshold=0.99 or 1)
mask1 = ~all(cube1_sub > 0.99, 3);
mask2 = ~all(cube2_sub > 0.99, 3);

%% Step 4: Make intersection mask of valid pixels
% (If the cropped regions are the same size, proceed. If not, further resize so they match.)
min_h = min(size(cube1_sub,1), size(cube2_sub,1));
min_w = min(size(cube1_sub,2), size(cube2_sub,2));
cube1_sub = cube1_sub(1:min_h,1:min_w,:);
cube2_sub = cube2_sub(1:min_h,1:min_w,:);
mask1 = mask1(1:min_h,1:min_w);
mask2 = mask2(1:min_h,1:min_w);

mask = mask1 & mask2; % Intersection

%% Step 5: Flatten masked pixels for both cubes
pixels1 = reshape(cube1_sub, [], size(cube1_sub,3));
pixels2 = reshape(cube2_sub, [], size(cube2_sub,3));
mask_flat = mask(:);
pixels1_masked = pixels1(mask_flat, :);
pixels2_masked = pixels2(mask_flat, :);

%% Step 6: K-means on film (cube1), apply indices to both
nColors = 60;


% Cluster only on before (film)
% [cluster_idx1, centers] = kmeans(pixels1_masked, nColors, 'Replicates', 5, 'MaxIter', 1000);
% 
% % Assign each after-pixel to its nearest before-cluster
% D2 = pdist2(pixels2_masked, centers); 
% [~, cluster_idx2_proj] = min(D2, [], 2);

% 1. Cluster on "before"
% [cluster_idx, cluster_centers] = kmeans(pixels1_masked, nColors, ...
%     'Replicates', 5, 'MaxIter', 1000);
% 
% % 2. Compute mean spectra for each cluster
% meanSpectra1 = zeros(nColors, size(pixels1_masked,2));
% meanSpectra2 = zeros(nColors, size(pixels2_masked,2));
% for k = 1:nColors
%     idx = (cluster_idx == k);
%     meanSpectra1(k,:) = mean(pixels1_masked(idx,:), 1, 'omitnan');
%     meanSpectra2(k,:) = mean(pixels2_masked(idx,:), 1, 'omitnan'); % <-- exact same indices!
% end

% (Assuming cube1_sub, cube2_sub, mask already cropped to same size)
nBlocksW = 10; nBlocksH = 6;
blockSizeH = floor(size(cube1_sub,1) / nBlocksH);
blockSizeW = floor(size(cube1_sub,2) / nBlocksW);

meanSpectra1 = zeros(nBlocksW*nBlocksH, size(cube1_sub,3));
meanSpectra2 = zeros(nBlocksW*nBlocksH, size(cube2_sub,3));
k = 1;
for row = 1:nBlocksH
    for col = 1:nBlocksW
        rows = (1:blockSizeH) + (row-1)*blockSizeH;
        cols = (1:blockSizeW) + (col-1)*blockSizeW;
        block_mask = mask(rows, cols);
        block_pixels1 = reshape(cube1_sub(rows, cols, :), [], size(cube1_sub,3));
        block_pixels2 = reshape(cube2_sub(rows, cols, :), [], size(cube2_sub,3));
        block_pixels1 = block_pixels1(block_mask(:),:);
        block_pixels2 = block_pixels2(block_mask(:),:);
        meanSpectra1(k,:) = mean(block_pixels1, 1, 'omitnan');
        meanSpectra2(k,:) = mean(block_pixels2, 1, 'omitnan');
        k = k+1;
    end
end





%%
% --- (1) Convert mean spectra to RGB

% Helper: function to turn spectra into RGB using MATLAB's built-in function
% (Assume reflectance, D65, 2° CMFs, wavelengths in nm)
% You can use your own process_rgb if you have it.

cmf = importdata('../../../data/CIE2degCMFs_1931.txt');  % or use your path
ill = importdata('../../../data/CIE_D65.txt');
cmf_wl = cmf(:,1);
cmf_xyz = cmf(:,2:4);
ill_wl = ill(:,1);
ill_val = ill(:,2);

nColors = size(meanSpectra1, 1);

% Interpolate illuminant and cmfs to match your wavelengths
ill_interp = interp1(ill_wl, ill_val, wl1, 'linear', 'extrap');
cmf_interp = interp1(cmf_wl, cmf_xyz, wl1, 'linear', 'extrap');

% Function to convert spectra to xyz then to rgb
spectra2rgb = @(refl) xyz2rgb((refl .* ill_interp') * cmf_interp ./ sum(ill_interp' .* cmf_interp(:,2)), 'ColorSpace','srgb');

% Convert each mean spectrum to RGB
RGBs1 = zeros(nColors, 3);
RGBs2 = zeros(nColors, 3);
for k = 1:nColors
    refl = meanSpectra1(k,:); % 1 x N
    weighted = refl .* ill_interp(:)'; % 1 x N
    xyz1 = weighted * cmf_interp; % 1 x 3
    xyz1 = xyz1 / sum(ill_interp(:)' .* cmf_interp(:,2)');
    xyz1 = max(0, xyz1);
    RGBs1(k,:) = xyz2rgb(xyz1, 'ColorSpace','srgb');

    refl = meanSpectra2(k,:);
    weighted = refl .* ill_interp(:)';
    xyz2 = weighted * cmf_interp;
    xyz2 = xyz2 / sum(ill_interp(:)' .* cmf_interp(:,2)');
    xyz2 = max(0, xyz2);
    RGBs2(k,:) = xyz2rgb(xyz2, 'ColorSpace','srgb');
end


% Clamp RGBs to [0 1]
RGBs1 = max(0, min(1, RGBs1));
RGBs2 = max(0, min(1, RGBs2));

% --- (2) Make a grid image for each palette


nColors = 60;
grid_w = 10;
grid_h = 6;
patchSize = 40;

% --- img1 and img2 initialised as white ---
img1 = ones(grid_h * patchSize, grid_w * patchSize, 3); % Cube 1
img2 = ones(grid_h * patchSize, grid_w * patchSize, 3); % Cube 2 (matched!)

for k = 1:nColors
    row = floor((k-1)/grid_w);
    col = mod((k-1), grid_w);
    r_idx = (row*patchSize+1):((row+1)*patchSize);
    c_idx = (col*patchSize+1):((col+1)*patchSize);
    patch_rgb1 = repmat(reshape(RGBs1(k,:), 1, 1, 3), patchSize, patchSize, 1);
    patch_rgb2 = repmat(reshape(RGBs2(k,:), 1, 1, 3), patchSize, patchSize, 1);
    img1(r_idx, c_idx, :) = patch_rgb1;               % before
    img2(r_idx, c_idx, :) = patch_rgb2;               % after (matched!)
end

%%
% --- (3) Show the palettes
figure; imshow(img1); title('Palette Cactus (Before)');
%%
figure; imshow(img2); title('Palette Cactus (After)');

%% Step 7: Save
save('palette_cactus_before.mat', 'meanSpectra1', 'wl1');
save('palette_cactus_after.mat', 'meanSpectra2', 'wl2');

%%

% --- DeltaE2000 Grid ---
dE = deltaE2000(Lab1, Lab2);
dE_grid = nan(grid_h, grid_w);
for k = 1:nColors
    row = floor((k-1)/grid_w)+1;
    col = mod((k-1), grid_w)+1;
    dE_grid(row, col) = dE(k);
end

figure; imagesc(dE_grid);
axis image off;
colormap(jet(255)); colorbar; clim([0 10]);
title('\DeltaE_{00} between palette patches');
for k = 1:nColors
    row = floor((k-1)/grid_w)+1;
    col = mod((k-1), grid_w)+1;
    if ~isnan(dE_grid(row,col))
        text(col, row, sprintf('%.1f', dE_grid(row,col)), ...
            'HorizontalAlignment','center','Color','w','FontSize',10,'FontWeight','bold');
    end
end

